# Copyright (C) 2012 Allis Tauri <allista@gmail.com>
# 
# degen_primer is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# degen_primer is distributed in the hope that it will be useful, but
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

import numpy as np
from array import array
from scipy.fftpack import fft, ifft
from SecStructures import Duplex
from StringTools import print_exception
from MultiprocessingBase import MultiprocessingBase
from WorkCounter import WorkCounter


class SearchEngine(MultiprocessingBase):
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
    ###########################################################################


    @classmethod
    def set_max_chunk(cls, chunk_size):
        cls._max_chunk_size = chunk_size
    

    @staticmethod
    def _map_letter(letter, _map):
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
    

    @staticmethod
    def _compile_duplexes_for_position(position, template, primer, 
                                           t_len, p_len, reverse):
        '''Given a template strand, a primer and a location where the 
        primer matches the template, return a list of Duplexes formed by 
        unambiguous components of the primer'''
        duplexes = []
        for var in primer.seq_records:
            dup = Duplex(str(var.seq), str(template[position:position+p_len].complement())[::-1])
            if not dup: continue
            duplexes.append((dup, var.id))
        if not reverse:
            return position+p_len, tuple(duplexes)
        else:
            return t_len+1-(position+p_len), tuple(duplexes)
    #end def
    _compile_duplexes_mapper = staticmethod(MultiprocessingBase.data_mapper(_compile_duplexes_for_position.__func__)) 
    
    
    @staticmethod
    @MultiprocessingBase.results_assembler
    def _duplexes_assembler(index, result, output):
        if result[1]: output.append(result)
    #end def
    

    def compile_duplexes_mp(self,  
                              fwd_seq, rev_seq, primer, 
                              t_len, p_len,
                              fwd_matches, rev_matches, counter=None):
        '''Compile duplexes for both strands of a template in parallel'''
        if counter is None: counter = WorkCounter()
        counter.set_subwork(2, map(len, (fwd_matches, rev_matches)))
        counter[0].set_work(len(fwd_matches))
        counter[1].set_work(len(rev_matches))
        #prepare and start two sets of jobs
        fwd_jobs = self.prepare_jobs(self._compile_duplexes_mapper, 
                                     fwd_matches, None, 
                                     fwd_seq, primer, t_len, p_len, False)
        rev_jobs = self.prepare_jobs(self._compile_duplexes_mapper, 
                                     rev_matches, None,
                                     rev_seq, primer, t_len, p_len, True)
        self.start_jobs(fwd_jobs,rev_jobs)
        #allocate containers for results
        fwd_results = []
        rev_results = []
        #assemble results
        self.get_result(fwd_jobs, 0.1, 
                       self._duplexes_assembler, fwd_results, counter=counter[0])
        self.get_result(rev_jobs, 0.1, 
                       self._duplexes_assembler, rev_results, counter=counter[1])
        #if aborted, return None
        if self._abort_event.is_set(): return None
        #else, cleanup
        self.clean_jobs(fwd_jobs+rev_jobs)
        del fwd_jobs; del rev_jobs 
        #sort duplexes by position and return them
        fwd_results.sort(key=lambda x: x[0])
        rev_results.sort(key=lambda x: x[0])
        return fwd_results,rev_results
    #end def
    
    
    def compile_duplexes(self,  
                           fwd_seq, rev_seq, primer, 
                           t_len, p_len,
                           fwd_matches, rev_matches, counter=None):
        '''Compile duplexes for both strands of a template'''
        if counter is None: counter = WorkCounter()
        counter.set_work(len(fwd_matches)+len(rev_matches))
        fwd_results = []
        for pos in fwd_matches:
            if self._abort_event.is_set(): break
            duplexes = self._compile_duplexes_for_position(pos, fwd_seq, primer, 
                                                           t_len, p_len, reverse=False)
            if duplexes[1]: fwd_results.append(duplexes)
            counter.count()
        rev_results = []
        for pos in rev_matches:
            if self._abort_event.is_set(): break
            duplexes = self._compile_duplexes_for_position(pos, rev_seq, primer, 
                                                           t_len, p_len, reverse=True)
            if duplexes[1]: rev_results.append(duplexes)
            counter.count()
        #if aborted, return None
        if self._abort_event.is_set(): return None
        return fwd_results,rev_results
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
        return max(chunk, min_chunk)
    #end def
    
    
    @staticmethod
    def _check_length_inequality(t_len, p_len):
        if t_len < p_len or p_len == 0:
            raise ValueError('SearchEngine._find: template sequence should be '
                             'longer or equal to primer sequence and both '
                             'should be grater than zero.')
    #end def
    
    
    @classmethod
    def mp_better(cls, t_len):
        #based on computation time statistics
        return cls.cpu_count > 1 and t_len > 25000
    #end def
    
    
    @classmethod
    def _optimal_slices(cls, t_len, p_len):
        #linear regression of measured computing time with respect to 
        #number of slices and template length
        linear = max(cls.cpu_count, int(t_len*1.75e-5 + 1.75)) 
        return min(60, linear, t_len/p_len)
    #end def
    
    
    @staticmethod
    @MultiprocessingBase.worker
    def _find_in_slice(abort_e, _id, start, end, fwd_seq, rev_seq, 
                         p_fft, correction, s_stride, c_size, c_stride):
        fwd_score = []
        rev_score = []
        pos = start
        while pos < end and not abort_e.is_set():
            front = min(end, pos+c_size)
            fwd_score.append(SearchEngine._find_in_chunk(fwd_seq[pos:front], 
                                                         p_fft, correction,
                                                         c_size, c_stride))
            rev_score.append(SearchEngine._find_in_chunk(rev_seq[pos:front], 
                                                         p_fft, correction,
                                                         c_size, c_stride))
            pos += c_stride
        fwd_score = np.concatenate(fwd_score)
        rev_score = np.concatenate(rev_score)
        return (start, 
                fwd_score[:s_stride], 
                rev_score[:s_stride])
    #end def
    

    def _find_mp(self, template, primer, t_len, p_len, mismatches):
        '''Find all occurrences of a primer sequence in both strands of a 
        template sequence with at most k mismatches. Multiprocessing version.'''
        slice_size   = t_len/self._optimal_slices(t_len, p_len)+p_len+1
        slice_stride = slice_size-p_len
        chunk_size   = self._calculate_chunk_size(slice_size, p_len)
        chunk_stride = chunk_size-p_len
        p_maps       = self._map_pattern(str(primer.master_sequence.seq), chunk_size)
        p_fft        = (fft(p_maps[0]),fft(p_maps[1]))
        fwd_seq      = template.seq
        rev_seq      = template.seq.reverse_complement()
        correction   = np.ndarray(chunk_stride); correction.fill(p_len/3.0)
        jobs         = []
        #start find_in_slice jobs
        pos = 0
        while pos < t_len and not self._abort_event.is_set():
            front = min(t_len, pos+slice_size)
            queue = self._Queue()
            job = self._Process(target=self._find_in_slice, 
                                args=(queue, self._abort_event, 
                                      len(jobs), 
                                      pos, front, 
                                      fwd_seq, rev_seq, 
                                      p_fft, correction,  
                                      slice_stride, chunk_size, chunk_stride))
            job.daemon = 1
            job.start()
            jobs.append((job,queue))
            self.register_job(jobs[-1])
            pos  += slice_stride
        #join all jobs; if scores arrays are allocated beforehand, the memory
        #is returned upon deletion
        scores_len = slice_stride*len(jobs)
        scores = [np.zeros(scores_len), np.zeros(scores_len)]
        def assemble_scores(out, scores):
            scores[0][out[0]:out[0]+slice_stride] = out[1]
            scores[1][out[0]:out[0]+slice_stride] = out[2]
        self.get_result(jobs, 1, assemble_scores, scores)
        #if search was aborted, return empty results
        if self._abort_event.is_set(): return None
        #else, cleanup
        self.clean_jobs(jobs); del jobs
        #compute match indices
        matches     = max(1, p_len - mismatches)-0.5
        fwd_matches = np.where(scores[0][:t_len-p_len+1] >= matches)[0]
        rev_matches = np.where(scores[1][:t_len-p_len+1] >= matches)[0] 
        return fwd_seq, rev_seq, primer, t_len, p_len, fwd_matches, rev_matches
    #end def

    
    def _find(self, template, primer, t_len, p_len, mismatches):
        '''Find all occurrences of a primer sequence in both strands of a 
        template sequence with at most k mismatches.'''
        chunk_size   = self._calculate_chunk_size(t_len, p_len)
        chunk_stride = chunk_size-p_len
        p_maps       = self._map_pattern(str(primer.master_sequence.seq), chunk_size)
        p_fft        = (fft(p_maps[0]),fft(p_maps[1]))
        fwd_seq      = template.seq
        rev_seq      = template.seq.reverse_complement()
        correction   = np.ndarray(chunk_stride); correction.fill(p_len/3.0)
        fwd_score    = np.zeros(t_len+chunk_stride)
        rev_score    = np.zeros(t_len+chunk_stride)
        #_find in chunks of a template, which is faster due to lower cost of memory allocation
        pos = 0
        while pos < t_len and not self._abort_event.is_set():
            front = min(t_len, pos+chunk_size)
            fwd_score[pos:pos+chunk_stride] = self._find_in_chunk(fwd_seq[pos:front], 
                                                                  p_fft, correction, 
                                                                  chunk_size, chunk_stride) 
            rev_score[pos:pos+chunk_stride] = self._find_in_chunk(rev_seq[pos:front], 
                                                                  p_fft, correction, 
                                                                  chunk_size, chunk_stride)
            pos += chunk_stride
        #if search was aborted, return empty results
        if self._abort_event.is_set(): return None
        #match indexes
        matches     = max(1, p_len - mismatches)-0.5
        fwd_matches = np.where(fwd_score[:t_len-p_len+1] >= matches)[0]; del fwd_score
        rev_matches = np.where(rev_score[:t_len-p_len+1] >= matches)[0]; del rev_score
        return fwd_seq, rev_seq, primer, t_len, p_len, fwd_matches, rev_matches
    #end def
    

    def find_matches(self, template, primer, mismatches):
        '''Find all occurrences of a primer sequence in both strands of a 
        template sequence with at most k mismatches. This method uses 
        multiprocessing to speedup the search process. Use it to perform 
        search in a long sequence.''' 
        p_len,t_len = len(primer),len(template)
        self._check_length_inequality(t_len, p_len)
        if self.mp_better(t_len): find = self._find_mp
        else: find = self._find
        return find(template, primer, t_len, p_len, mismatches)
    #end def
    

    def find(self, counter, template, primer, mismatches):
        '''Find occurrences of a degenerate primer in a template sequence.
        Return positions and Duplexes formed. This method uses 
        multiprocessing to speedup the search process. Use it to perform 
        search in a long sequence.''' 
        matches = self.find_matches(template, primer, mismatches)
        if matches is None: return None
        return self.compile_duplexes_mp(*matches, counter=counter)
    #end def
    
    
    def _batch_find(self, template, primer, p_len, mismatches):
        t_id, t_len, _template = template
        matches = self._find(_template, primer, t_len, p_len, mismatches)
        if matches is None: return t_id, None
        return t_id, self.compile_duplexes(*matches)
    #end def
    
    def batch_find(self, counter, templates, primer, mismatches):
        '''Find occurrences of a degenerate primer in each of the provided 
        templates. Return dictionary of results using template IDs as keys.
        It uses multiprocessing to parallelize searches, each of which does 
        not use parallelization. Use it to search in many short sequences.'''
        results = self.parallelize_work(0.1, self._batch_find,
                                        templates, primer, len(primer), 
                                        mismatches, counter=counter)
        if results is None: return None
        return dict(results)
    #end def
#end class

#shortcut functions
mp_better = SearchEngine.mp_better

def find(counter, template, primer, mismatches, abort_event):
    return SearchEngine(abort_event).find(counter, template, primer, mismatches)

def batch_find(counter, templates, primer, mismatches, abort_event):
    return SearchEngine(abort_event).batch_find(counter, templates, primer, mismatches)


#tests
if __name__ == '__main__':
    #tests
    import signal
    import cProfile
    import gc
    from time import time
    import sys, os, csv, timeit
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC
    from Primer import Primer
    from multiprocessing import Event
    
    searcher = None
    ppid     = -1
    data1    = []
    data2    = []
    abort_event = Event()
    
    def sig_handler(signal, frame):
        if ppid != os.getpid():
            return
        print 'Aborting main process %d' % os.getpid()
        abort_event.set()
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

    #setup signal handler
    signal.signal(signal.SIGINT,  sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)
    signal.signal(signal.SIGQUIT, sig_handler)

    seq_file = '../ThGa.fa'
    try:
        record_file = open(seq_file, 'r')
    except IOError, e:
        print 'Unable to open %s' % seq_file
        print_exception(e)
        sys.exit(1)
    template = SeqIO.read(record_file, 'fasta', IUPAC.unambiguous_dna)
    record_file.close()
    
    ftgam = Seq('ATATTCTACRACGGCTATCC', IUPAC.ambiguous_dna)
    rtgam = Seq('GAASGCRAAKATYGGGAAC', IUPAC.ambiguous_dna)

#    primer   = Primer(SeqRecord(ftgam, id='ftgam'), 0.43e-6, True)
    primer   = Primer(SeqRecord(rtgam, id='rtgam'), 0.43e-6, True) #48C, 9mism, results: 564.85Mb, 410.37Mb, 240.07Mb
    
    def print_out(out, name):
        print name
        for i, results in enumerate(out):
            if not results: 
                print 'No results.'
            else:
                print 'Results %d' % i
                for res in results:
                    print res[0]
                    for dup, _id in res[1]:
                        print _id
                        print dup
                        #print 'mismatches:', dup.mismatches
                    print '\n'
                print '\n'
    #end def

    ppid = os.getpid()
    import TD_Functions
    TD_Functions.PCR_P.PCR_T = 60
    
    searcher = SearchEngine(abort_event)
    
    def mem_test(num):
        for _n in xrange(num):
            t0 = time()
            c = WorkCounter()
            searcher = SearchEngine(abort_event)
            results = searcher.find(c, template, primer, 4)
            t1 = (time()-t0)
#            print 'results use: %fMb' % (asizeof(results)/1024.0/1024.0)
            print results
#            print_out(results, '')
            print '[%6.2f%%] elapsed %f seconds\n' % (c.percent(), t1)
    #end def
    
    mem_test(1)
    sys.exit()
        
    
    #compile duplexes serial vs parallel statistics
    fwd_seq = template.seq
    rev_seq = template.seq.reverse_complement()
    t_len = len(template); p_len = len(primer)
    x = []; s = []; p = []
    for n in xrange(5,105,5):
        for _i in xrange(10):
            x.append(n)
            fwd_matches = np.random.randint(0, t_len, n/2)
            rev_matches = np.random.randint(0, t_len, n-n/2)
            s.append(timeit.timeit('searcher.compile_duplexes(fwd_seq, rev_seq, primer, t_len, p_len, fwd_matches, rev_matches)', 
                                   'from __main__ import searcher, fwd_seq, rev_seq, t_len, p_len, primer, fwd_matches, rev_matches', 
                                   number=1))
            p.append(timeit.timeit('searcher.compile_duplexes_mp(fwd_seq, rev_seq, primer, t_len, p_len, fwd_matches, rev_matches)', 
                                   'from __main__ import searcher, fwd_seq, rev_seq, t_len, p_len, primer, fwd_matches, rev_matches', 
                                   number=1))
    
    from matplotlib import pyplot as plt
    from scipy.stats import linregress
    x = np.fromiter(x, dtype=int)
    slm = linregress(x, s)
    plm = linregress(x, p)
    plt.plot(x, s, 'b.', x, x*slm[0]+slm[1], 'b-', 
             x, p, 'r.', x, x*plm[0]+plm[1], 'r-')
    plt.show()
