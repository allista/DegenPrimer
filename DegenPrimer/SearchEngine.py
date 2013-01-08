# Copyright (C) 2012 Allis Tauri <allista@gmail.com>
# 
# degen_primer is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# indicator_gddccontrol is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
Created on Jan 1, 2013

@author: Allis Tauri <allista@gmail.com>
'''

import os
import errno
import numpy as np
import multiprocessing as mp
from Queue import Empty
from time import time, sleep
from array import array
from scipy.fftpack import fft, ifft
from SecStructures import Duplex
from StringTools import print_exception


class SearchEngine(object):
    '''Fast search of a pattern sequence in a given set of template sequences 
    with parallelization using multiprocessing.'''
    
    _w_3_0 = 1                        #trivial 3d root of unity 
    _w_3_1 = (-1/2.0+np.sqrt(3)/2.0j) #3d root of unity in a power of 1
    _w_3_2 = (-1/2.0-np.sqrt(3)/2.0j) #3d root of unity in a power of 2
    
    _unambiguous = array('b','ATGC')  #ATGC characters as byte array
    
    #alphabet mappings to the 3d roots of unity for templates and patterns
    _T_AT_mapping = dict(zip(_unambiguous, (_w_3_1,_w_3_2,0,0))) 
    
    _T_GC_mapping = dict(zip(_unambiguous, (0,0,_w_3_1,_w_3_2)))

    _P_AT_mapping = {'A':_w_3_2,
                     'T':_w_3_1,
                     'G':_w_3_0,
                     'C':_w_3_0,
                     
                     'R':_w_3_2,
                     'Y':_w_3_1,
                     'S':_w_3_0,
                     'W':_w_3_2+_w_3_1,
                     'K':_w_3_1,
                     'M':_w_3_2,
                     'B':_w_3_1,
                     'D':_w_3_2+_w_3_1,
                     'H':_w_3_2+_w_3_1,
                     'V':_w_3_2,
                     'N':_w_3_2+_w_3_1}
    
    _P_GC_mapping = {'A':_w_3_0,
                     'T':_w_3_0,
                     'G':_w_3_2,
                     'C':_w_3_1,
                     
                     'R':_w_3_2,
                     'Y':_w_3_1,
                     'S':_w_3_2+_w_3_1,
                     'W':_w_3_0,
                     'K':_w_3_2,
                     'M':_w_3_1,
                     'B':_w_3_2+_w_3_1,
                     'D':_w_3_2,
                     'H':_w_3_1,
                     'V':_w_3_2+_w_3_1,
                     'N':_w_3_2+_w_3_1}
    
    #template sequences grater than this value will be split into slices of 
    #approximately that size 
    _max_chunk_size = 2**12
    
    #maximum jobs to launch within a single search
    _cpu_count   = mp.cpu_count()
    
    ###########################################################################


    def __init__(self):
        self._all_jobs        = []
        self._abort_lock      = mp.Lock()
        self._abort_events    = []
    #end def
    
    
    @classmethod
    def set_max_chunk(cls, chunk_size):
        cls._max_chunk_size = chunk_size
    

    @classmethod
    def _map_letter(cls, letter, _map):
        try: return _map[letter]
        except KeyError: return 0
    #end def
    
    
    @classmethod
    def _map_pattern(cls, pattern, map_len):
        '''Map pattern sequence to an alphabet of 3d roots of unity.'''
        #this naive algorithm works ~5 times faster 
        #than the one used for template mapping due to short length of patterns
        AT_map = np.zeros(map_len, dtype=complex)
        GC_map = np.zeros(map_len, dtype=complex)
        for i,letter in enumerate(pattern):
            AT_map[i] = cls._map_letter(letter, cls._P_AT_mapping)
            GC_map[i] = cls._map_letter(letter, cls._P_GC_mapping)
        return (AT_map, GC_map)
    #end def


    @classmethod
    def _compile_duplexes(cls, template, primer, matches):
        '''Given a template strand, a primer and a list of locations where the 
        primer matches the template, return a list of Duplexes formed by 
        unambiguous components of the primer at each match location.'''
        p_len = len(primer)
        results = []
        for i in matches:
            duplexes = []
            for var in primer.sequences:
                duplexes.append(Duplex(var, template[i:i+p_len].reverse_complement()))
            results.append((i+p_len-1, duplexes))
        return results
    #end def
    
    
    @classmethod
    def _find_in_chunk(cls, t_chunk, p_fft, correction, c_size, c_stride):
        '''Find number of matches of pattern at each position in a given 
        chunk of a template.
        Pattern is given as a polynomial evaluated at n-th roots of unity 
        using fft.
        map_len is a length of a map of template to the alphabet of 3-d roots 
        of unity chunk to build. It's a power of 2 integer for fft to work fast.
        c_stride -- a part of chunk for which matches are calculated 
        (it is less than map_len, so chunks overlap each other)'''
        t_AT_map = np.fromiter(array('b',str(t_chunk)), dtype=complex) 
        t_AT_map.resize(c_size)
        t_GC_map = t_AT_map.copy()
        for k,v in cls._T_AT_mapping.iteritems():
            t_AT_map[t_AT_map == k] = v
        for k,v in cls._T_GC_mapping.iteritems():
            t_GC_map[t_GC_map == k] = v
        AT_score = ifft(fft(t_AT_map[::-1])*p_fft[0])[::-1][:c_stride]
        GC_score = ifft(fft(t_GC_map[::-1])*p_fft[1])[::-1][:c_stride]
        score    = AT_score.real + GC_score.real
        score    = (score + correction - score/3.0)
        return score
    #end def
    
    
    @classmethod
    def _calculate_chunk_size(cls, t_len, p_len):
        rem = lambda(c): ((t_len/(c-p_len))*(c-p_len)+c-t_len)
        if t_len <= cls._max_chunk_size:
            chunk = 2**int(np.ceil(np.log2(t_len)))
            if chunk % t_len == 0: return chunk
        else: chunk = cls._max_chunk_size
        r = rem(chunk)
        min_chunk = 2**int(np.ceil(np.log2(2*p_len)))
        max_rem   = chunk/2+1
        while r > max_rem \
        and chunk > min_chunk:
            chunk /= 2
            r = rem(chunk)
        return chunk
    #end def
    
    
    @classmethod
    def _check_length_inequality(cls, t_len, p_len):
        if t_len < p_len or p_len == 0:
            raise ValueError('SearchEngine.find: template sequence should be '
                             'longer or equal to primer sequence and both '
                             'should be grater than zero.')

    
    @classmethod
    def mp_better(cls, t_len, p_len):
        return cls._cpu_count > 1 and t_len > 25000


    @classmethod
    def _join_jobs(cls, abort_e, jobs, parser, *args):
        while jobs and not abort_e.is_set():
            finished_job = None
            for i,job in enumerate(jobs):
                try: 
                    out = job[1].get(True,1e-3)
                    parser(out, *args)
                    job[0].join()
                    finished_job = i
                    break
                except Empty: 
                    continue
                except IOError, e:
                    if e.errno == errno.EINTR:
                        continue
                    else:
                        print 'Unhandled IOError:', e.message
                        raise
                except Exception, e:
                    print 'Unhandled Exception:'
                    print_exception(e)
                    raise
            if finished_job is not None:
                del jobs[finished_job]
    #end def
    

    @classmethod
    def _start_find_worker(cls, jobs, queue, abort_e, start, fwd_chunk, rev_chunk, 
                           p_fft, correction, t_len, p_len, c_size, c_stride):
        def worker(queue, abort_e, start, fwd_chunk, rev_chunk, p_fft, 
                   correction, t_len, p_len, c_size, c_stride):
            fwd_score = np.ndarray(0,dtype=float)
            rev_score = np.ndarray(0,dtype=float)
            pos = 0
            while pos < t_len and not abort_e.is_set():
                front = min(t_len, pos+c_size)
                score = cls._find_in_chunk(fwd_chunk[pos:front], p_fft, correction,
                                           c_size, c_stride)
                fwd_score = np.concatenate([fwd_score, score])
                score = cls._find_in_chunk(rev_chunk[pos:front], p_fft, correction,
                                           c_size, c_stride)
                rev_score = np.concatenate([rev_score, score])
                pos += c_stride
            if abort_e.is_set():
                queue.cancel_join_thread()
            else: queue.put((start, 
                             fwd_score[:t_len-p_len], 
                             rev_score[:t_len-p_len]))
            queue.close()
            return
        #end def
        job = mp.Process(target=worker, args=(queue, abort_e, start, fwd_chunk, rev_chunk, p_fft, 
                                              correction, t_len, p_len, c_size, c_stride))
        job.start()
        jobs.append((job,queue))
    #end def
    
    
    def find_mp_test(self, template, primer, mismatches, mult):
        '''Find all occurrences of a primer sequence in both strands of a 
        template sequence with at most k mismatches. Multiprocessing version.'''
        p_len,t_len  = len(primer),len(template)
        self._check_length_inequality(t_len, p_len)
        slice_size   = t_len/mult+p_len+1
        slice_stride = slice_size-p_len
        chunk_size   = self._calculate_chunk_size(slice_size, p_len)
        chunk_stride = chunk_size-p_len
        p_maps       = self._map_pattern(str(primer.master_sequence.seq), chunk_size)
        p_fft        = (fft(p_maps[0]),fft(p_maps[1]))
        fwd_seq      = template.seq
        rev_seq      = template.seq.reverse_complement()
        correction   = np.ndarray(chunk_stride); correction.fill(p_len/3.0)
        fwd_score    = []
        rev_score    = []
        jobs         = []
        i = 0; abort_event = mp.Event()
        self._abort_events.append(abort_event)
        #start find_in_chunk jobs
        while i < t_len and not abort_event.is_set():
            front = min(t_len, i+slice_size)
            self._start_find_worker(jobs, mp.Queue(), abort_event, i, 
                                    fwd_seq[i:front], rev_seq[i:front], p_fft, 
                                    correction, front-i, p_len, chunk_size, chunk_stride)
            self._all_jobs.append(jobs[-1])
            i += slice_stride
        #join all jobs
        def parse_out(out, fwd_list, rev_list):
            fwd_list.append((out[0],out[1]))
            rev_list.append((out[0],out[2]))
        self._join_jobs(abort_event, list(jobs), parse_out, fwd_score, rev_score)
        #if search was aborted, return empty results
        if abort_event.is_set():
            return [],[]
        #else, cleanup
        with self._abort_lock:
            if abort_event in self._abort_events:
                self._abort_events.remove(abort_event)
            for job in jobs: 
                if job in self._all_jobs:
                    self._all_jobs.remove(job)
        #sort, correct and concatenate scores
        fwd_score.sort(key=lambda(x): x[0])
        fwd_score = [result[1] for result in fwd_score]
        rev_score.sort(key=lambda(x): x[0])
        rev_score = [result[1] for result in rev_score]
        fwd_score = np.concatenate(fwd_score)
        rev_score = np.concatenate(rev_score)
        #match indices
        matches     = max(1, p_len - mismatches)-0.5
        fwd_matches = np.arange(t_len-p_len)[fwd_score[:t_len-p_len] >= matches]
        rev_matches = np.arange(t_len-p_len)[rev_score[:t_len-p_len] >= matches]
        del fwd_score, rev_score
        #construct duplexes
        fwd_results = self._compile_duplexes(fwd_seq, primer, fwd_matches)
        rev_results = self._compile_duplexes(rev_seq, primer, rev_matches)
        return fwd_results,rev_results
    #end def
    
    
    @classmethod
    def _optimal_slices(cls, t_len):
        linear = max(cls._cpu_count, int(t_len*1.75e-5 + 1.75))
        return min(60, linear)
    #end def
    
    
    def find_mp(self, template, primer, mismatches, abort_event=None):
        '''Find all occurrences of a primer sequence in both strands of a 
        template sequence with at most k mismatches. Multiprocessing version.'''
        p_len,t_len  = len(primer),len(template)
        self._check_length_inequality(t_len, p_len)
        slice_size   = t_len/self._optimal_slices(t_len)+p_len+1
        slice_stride = slice_size-p_len
        chunk_size   = self._calculate_chunk_size(slice_size, p_len)
        chunk_stride = chunk_size-p_len
        p_maps       = self._map_pattern(str(primer.master_sequence.seq), chunk_size)
        p_fft        = (fft(p_maps[0]),fft(p_maps[1]))
        fwd_seq      = template.seq
        rev_seq      = template.seq.reverse_complement()
        correction   = np.ndarray(chunk_stride); correction.fill(p_len/3.0)
        fwd_score    = []
        rev_score    = []
        jobs         = []
        if abort_event is None:
            abort_event = mp.Event()
            self._abort_events.append(abort_event)
        #start find_in_chunk jobs
        i = 0 
        while i < t_len and not abort_event.is_set():
            front = min(t_len, i+slice_size)
            self._start_find_worker(jobs, mp.Queue(), abort_event, i, 
                                    fwd_seq[i:front], rev_seq[i:front], p_fft, 
                                    correction, front-i, p_len, chunk_size, chunk_stride)
            self._all_jobs.append(jobs[-1])
            i += slice_stride
        #join all jobs
        def parse_out(out, fwd_list, rev_list):
            fwd_list.append((out[0],out[1]))
            rev_list.append((out[0],out[2]))
        self._join_jobs(abort_event, list(jobs), parse_out, fwd_score, rev_score)
        #if search was aborted, return empty results
        if abort_event.is_set():
            return [],[]
        #else, cleanup
        with self._abort_lock:
            if abort_event in self._abort_events:
                self._abort_events.remove(abort_event)
            for job in jobs: 
                if job in self._all_jobs:
                    self._all_jobs.remove(job)
        #sort, correct and concatenate scores
        fwd_score.sort(key=lambda(x): x[0])
        fwd_score = [result[1] for result in fwd_score]
        rev_score.sort(key=lambda(x): x[0])
        rev_score = [result[1] for result in rev_score]
        fwd_score = np.concatenate(fwd_score)
        rev_score = np.concatenate(rev_score)
        #match indices
        matches     = max(1, p_len - mismatches)-0.5
        fwd_matches = np.arange(t_len-p_len)[fwd_score[:t_len-p_len] >= matches]
        rev_matches = np.arange(t_len-p_len)[rev_score[:t_len-p_len] >= matches]
        del fwd_score, rev_score
        #construct duplexes
        fwd_results = self._compile_duplexes(fwd_seq, primer, fwd_matches)
        rev_results = self._compile_duplexes(rev_seq, primer, rev_matches)
        return fwd_results,rev_results
    #end def

    
    def find(self, template, primer, mismatches, abort_event=None):
        '''Find all occurrences of a primer sequence in both strands of a 
        template sequence with at most k mismatches.'''
        p_len,t_len  = len(primer),len(template)
        self._check_length_inequality(t_len, p_len)
        chunk_size   = self._calculate_chunk_size(t_len, p_len)
        chunk_stride = chunk_size-p_len
        p_maps       = self._map_pattern(str(primer.master_sequence.seq), chunk_size)
        p_fft        = (fft(p_maps[0]),fft(p_maps[1]))
        fwd_seq      = template.seq
        rev_seq      = template.seq.reverse_complement()
        fwd_score    = []
        rev_score    = []
        correction   = np.ndarray(chunk_stride); correction.fill(p_len/3.0)
        if abort_event is None:
            abort_event = mp.Event()
            self._abort_events.append(abort_event)
        #find in chunks of a template, which is faster due to lower cost of memory allocation
        i = 0
        while i < t_len and not abort_event.is_set():
            front = min(t_len, i+chunk_size)
            fwd_score.append(self._find_in_chunk(fwd_seq[i:front], p_fft, correction, 
                                                chunk_size, chunk_stride))
            
            rev_score.append(self._find_in_chunk(rev_seq[i:front], p_fft, correction, 
                                                chunk_size, chunk_stride))
            i += chunk_stride
        #if search was aborted, return empty results
        if abort_event.is_set(): return [],[]
        self._abort_events.remove(abort_event)
        #concatenate scores
        fwd_score = np.concatenate(fwd_score)
        rev_score = np.concatenate(rev_score)
        #match indices
        matches     = max(1, p_len - mismatches)-0.5
        fwd_matches = np.arange(t_len-p_len)[fwd_score[:t_len-p_len] >= matches]
        rev_matches = np.arange(t_len-p_len)[rev_score[:t_len-p_len] >= matches]
        del fwd_score, rev_score
        #construct duplexes
        fwd_results = self._compile_duplexes(fwd_seq, primer, fwd_matches)
        rev_results = self._compile_duplexes(rev_seq, primer, rev_matches)
        return fwd_results,rev_results
    #end def
    
    
    def _start_short_find_job(self, jobs, queue, abort_event, t_id, template, primer, mismatches):
        def worker(searcher, queue, abort_event, t_id, template, primer, mismatches):
            result = searcher.find(template, primer, mismatches, abort_event)
            if abort_event.is_set():
                queue.cancel_join_thread()
            else: queue.put((t_id, result))
            queue.close()
            return
        #end def
        job = mp.Process(target=worker, args=(self, queue, abort_event, 
                                              t_id, template, primer, mismatches))
        job.start()
        jobs.append((job,queue))
    #end def
    
    
    def batch_find(self, templates, primer, mismatches):
        results = dict()
        p_len   = len(primer)
        #sort templates into short, suitable for plain search and long -- for mp 
        short_templates = []
        long_templates  = []
        for t in templates:
            if self.mp_better(t[0], p_len):
                long_templates.append(t)
            else: short_templates.append(t)
        #add aboty event
        abort_event = mp.Event()
        self._abort_events.append(abort_event)
        #first, launch short jobs
        jobs = []
        for t_len, t_id, template in short_templates:
            if abort_event.is_set(): break
            self._start_short_find_job(jobs, mp.Queue(), abort_event, t_id, 
                                       template, primer, mismatches)
        #join all launched jobs and gather the results
        def parse_out(out, results):
            results[out[0]] = out[1]
        self._join_jobs(abort_event, list(jobs), parse_out, results)
        #when they are finished, launch long searches one by one
        for t_len, t_id, template in long_templates:
            if abort_event.is_set(): break
            results[t_id] = self.find_mp(template, primer, mismatches, abort_event)
        #if abort event is set, return empty dict
        if abort_event.is_set(): 
            return dict().fromkeys([t[1] for t in templates])
        #else, cleanup
        with self._abort_lock:
            if abort_event in self._abort_events:
                self._abort_events.remove(abort_event)
            for job in jobs: 
                if job in self._all_jobs:
                    self._all_jobs.remove(job)
        return results
    #end def
    
    
    def abort(self):
        if self._abort_events:
            with self._abort_lock:
                for event in self._abort_events:
                    event.set()
            sleep(1)
            with self._abort_lock:
                while self._all_jobs:
                    finished_job = None
                    for i,job in enumerate(self._all_jobs):
                        try: 
                            while True: job[1].get_nowait()
                        except Empty: pass
                        if not job[0].is_alive():
                            job[0].join()
                            finished_job = i
                            break
                    if finished_job is not None:
                        del self._all_jobs[finished_job]
                self._abort_events = []
    #end def
#end class



#tests
import signal
import cProfile
import sys, csv, timeit
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Primer import Primer

searcher = None
ppid     = -1
data1    = []
data2    = []

def sig_handler(signal, frame):
    if ppid != os.getpid():
        return
    print 'Aborting main process %d' % os.getpid()
    if searcher:
        searcher.abort()
    if data1:
        print 'Write out gathered data1...'
        out_file = open('gather_data1-%d.csv' % time(), 'wb')
        csv_writer = csv.writer(out_file, delimiter='\t', quotechar='"')
        csv_writer.writerows(data1)
        out_file.close()
        print 'Done.'
    if data2:
        print 'Write out gathered data2...'
        out_file = open('gather_data2-%d.csv' % time(), 'wb')
        csv_writer = csv.writer(out_file, delimiter='\t', quotechar='"')
        csv_writer.writerows(data2)
        out_file.close()
        print 'Done.'
    sys.exit(1)
#end def


if __name__ == '__main__':
    #setup signal handler
    signal.signal(signal.SIGINT,  sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)
    signal.signal(signal.SIGQUIT, sig_handler)

    os.chdir('../')
    try:
        record_file = open('Ch5_gnm.fa', 'r')
    except IOError, e:
        print 'Unable to open Ch5_gnm.fa'
        print_exception(e)
        sys.exit(1)
    record = SeqIO.read(record_file, 'fasta', IUPAC.unambiguous_dna)
    record_file.close()
    
    query = Seq('ATATTCTACRACGGCTATCC', IUPAC.ambiguous_dna)
    
    template = record+record+record+record
    query    = query+query+query+query+query
    searcher = SearchEngine()

    def print_out(out, name):
        print name
        for results in out:
            if not results: 
                print 'No results.'
            else:
                for res in results:
                    print res[0]
                    for dup in res[1]:
                        print dup
                    print '\n'
                print '\n'
    #end def
    
    
    T,P,mult = 0,0,0
    def gather_statistics(start_len, end_len, delta):
        global ppid, data1, data2, T,P,mult
        ppid  = os.getpid()
        mults = range(searcher._cpu_count, max(searcher._cpu_count, end_len/50000)+10)
        patts = range(20,21,5)
        data1  = [['multiplier']]
        data1 += [[mult] for mult in range(1,max(mults)+1)]
        data2  = [('t_len', 'p_len', 'min multiplier', 'min time (sec)')]
        for t_len in range(start_len, end_len+1, delta):
            print 't_len:',t_len
            l_mults = range(searcher._cpu_count, max(searcher._cpu_count,t_len/50000)+10)
            for p_len in patts:
                print 'p_len:',p_len
                data1[0].append('%d-%d' % (t_len, p_len))
                times = dict([(m,[]) for m in l_mults])
                for i in range(5):
                    print 'iteration:', i
                    for mult in l_mults:
                        print 'mult:',mult
                        T = template[:t_len]
                        P = Primer(SeqRecord(query[:p_len], id='test'), 0.1e-6)
                        et = timeit.timeit('searcher.find_mp_test(T, P, 3, mult)', 
                                           setup='from __main__ import T,P,mult,searcher', 
                                           number=1)
                        times[mult].append(et)
                        print ''
                for mult in mults:
                    if mult in l_mults:
                        data1[mult].append(min(times[mult]))
                    else:
                        data1[mult].append(None)
                min_mult = min(l_mults)
                min_time = data1[min_mult][-1]
                for mult in l_mults:
                    if  data1[mult][-1] < min_time:
                        min_time = data1[mult][-1]
                        min_mult = mult
                data2.append((t_len, p_len, min_mult, min_time))
        filename1 = 'find_mp_data1-%d-%d-%d.csv' % (start_len, end_len, time())
        filename2 = 'find_mp_data2-%d-%d-%d.csv' % (start_len, end_len, time())
        out_file = open(filename1, 'wb')
        csv_writer = csv.writer(out_file, delimiter='\t', quotechar='"')
        csv_writer.writerows(data1)
        out_file.close()
        out_file = open(filename2, 'wb')
        csv_writer = csv.writer(out_file, delimiter='\t', quotechar='"')
        csv_writer.writerows(data2)
        out_file.close()
        data1 = []
        data2 = []
        print 'Gather statistics: data was written to %s and %s' % (filename1, filename2)
    #end def
    
    
    def gather_statistics1(start_len, end_len, delta, title=None):
        global ppid, data1, data2, T,P,mult
        ppid  = os.getpid()
        p_len = 20
        lens  = range(start_len, end_len+1, delta)
        if title:
            data2 = [('t_len', 'p_len', 'find time %s'%title, 
                      'find_mp time %s'%title)]
        else: data2 = [('t_len', 'p_len', 'find time', 'find_mp time')]
        find_times    = dict([(l,[]) for l in lens])
        find_mp_times = dict([(l,[]) for l in lens])
        for i in range(20):
            print 'iteration', i
            for t_len in lens:
                print 't_len:',t_len
                T = template[:t_len]
                P = Primer(SeqRecord(query[:p_len], id='test'), 0.1e-6)
                et = timeit.timeit('searcher.find(T, P, 3)', 
                                   setup='from __main__ import T,P,searcher', 
                                   number=1)
                find_times[t_len].append(et)
                et = timeit.timeit('searcher.find_mp(T, P, 3)', 
                                   setup='from __main__ import T,P,searcher', 
                                   number=1)
                find_mp_times[t_len].append(et)
                print ''
        for t_len in lens:
            data2.append((t_len, p_len, 
                          min(find_times[t_len]), min(find_mp_times[t_len])))
        if title:
            filename = 'find_mp_vs_find-%s-%d.csv' % (title, time())
        else: filename = 'find_mp_vs_find-%d.csv' % time()
        out_file = open(filename, 'wb')
        csv_writer = csv.writer(out_file, delimiter='\t', quotechar='"')
        csv_writer.writerows(data2)
        out_file.close()
        data2 = []
        print 'Gather statistics 1: data was written to %s' % filename 
    #end def

    ppid = os.getpid()
    
    
    gather_statistics1(5000, 500000, 5000)
    
#    for l in [2**10, 2**11, 2**12, 2**13, 2**14, 2**15, 2**16]:
#        searcher.set_max_chunk(l)
#        gather_statistics1(205000, 400000, 5000, str(l))
        
#    while True:
#        gather_statistics(50000, 500000, 50000)
#        gather_statistics(600000, 1000000, 100000)
#        gather_statistics(1100000, 1500000, 100000)
#        gather_statistics(1600000, 2000000, 100000)
#        gather_statistics(1800000, 2600000, 200000)
#        gather_statistics(2800000, 4000000, 400000)
        #gather_statistics(5000000, 10000000, 1000000)
        
    #searcher.find(template[:215000], Primer(SeqRecord(query[:20], id='test'), 0.1e-6), 3)
    #searcher.find(template[:220000], Primer(SeqRecord(query[:20], id='test'), 0.1e-6), 3)

    #primer = Primer(SeqRecord(query, id='test'), 0.1e-6)
    #out0,out1,out2,out3 = [],[],[],[]
    #cProfile.run("for i in range(10): searcher.find_mp(template, primer, 3);", 'find_mp2.profile')
    #sleep(1)
    #cProfile.run("for i in range(10): searcher.find(template[:33000], Primer(SeqRecord(query[:15], id='test'), 0.1e-6), 3)", 'find_33.profile')
    #cProfile.run("for i in range(10): searcher.find(template[:34000], Primer(SeqRecord(query[:15], id='test'), 0.1e-6), 3)", 'find_34.profile')

    
#    out  = None
#    cProfile.run("for i in range(10):\
#        out = searcher.find(template[:215000], Primer(SeqRecord(Seq('ATATTCTACRACGGCTATCC', IUPAC.ambiguous_dna), id='test'), 0.1e-6), 3)",
#        'find_215000.profile')
#    print_out(out, 'find')
#    cProfile.run("for i in range(10):\
#        out = searcher.find(template[:220000], Primer(SeqRecord(Seq('ATATTCTACRACGGCTATCC', IUPAC.ambiguous_dna), id='test'), 0.1e-6), 3)",
#        'find_220000.profile')
#    print_out(out, 'find')
#    ar1, ar2 = None,None
#    for l in [2**10, 2**11, 2**12, 2**13, 2**14, 2**15, 2**16]:
#        ar1 = np.random.random_sample(l).astype(complex)
#        ar2 = np.random.random_sample(l).astype(complex)
#        times = []
#        for i in range(10):
#            et = timeit.timeit('ifft(fft(ar1)*fft(ar2))', 
#                               setup='from __main__ import ar1, ar2\n'
#                                     'from scipy.fftpack import fft, ifft', 
#                               number=10)
#            times.append(et)
#        print '%d %f' % (l, min(times))
#        cProfile.run("for i in range(10):\
#            ifft(fft(ar1)*fft(ar2))",
#            'ifft_fft_%d.profile'%l)
    
    
    #print_out(out0, 'find')
    #print_out(out1, 'find_mp')
    print 'Done'
