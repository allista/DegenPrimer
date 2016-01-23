# coding=utf-8
#
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
Created on 2016-01-14

@author: Allis Tauri <allista@gmail.com>
'''

def test():
    #tests
    import signal
    from time import time
    import sys, os, csv
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC
    from DegenPrimer.Primer import Primer
    from multiprocessing import Event
    import cProfile
    
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

    seq_file = '../data/ThGa.fa'
    try:
        record_file = open(seq_file, 'r')
    except IOError, e:
        print 'Unable to open %s' % seq_file
        print e
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
    import DegenPrimer.TD_Functions as tdf
    from DegenPrimer.SearchEngine import SearchEngine
    from DegenPrimer.WorkCounter import WorkCounter
    tdf.PCR_P.PCR_T = 60
    
    searcher = SearchEngine(abort_event)
    
    cProfile.runctx('for x in xrange(10): searcher._find(WorkCounter(), template, primer, len(template), len(primer), 6)', 
                    globals(), locals(), 'SearchEngine._find.profile')
    
    matches = searcher._find_mp(WorkCounter(), template[:2000], primer, len(template), len(primer), 6)
    cProfile.runctx('for x in xrange(10): print searcher.compile_duplexes(WorkCounter(), *matches)', 
                    globals(), locals(), 'SearchEngine.compile_duplexes.profile')
    
    def mem_test(num):
        for _n in xrange(num):
            t0 = time()
            c = WorkCounter()
            searcher = SearchEngine(abort_event)
            results = searcher.find(c, template, primer, 6)
            t1 = (time()-t0)
#            print 'results use: %fMb' % (asizeof(results)/1024.0/1024.0)
            print results
#            print_out(results, '')
            print '[%6.2f%%] elapsed %f seconds\n' % (c.percent(), t1)
    #end def
    
#    mem_test(1)
#    sys.exit()
        
    
    #compile duplexes serial vs parallel statistics
#    fwd_seq = template.seq
#    rev_seq = template.seq.reverse_complement()
#    t_len = len(template); p_len = len(primer)
#    n = 1000; j = []; s = []; p = []
#    jobs = (4,32,4); tries = 1
#    c = WorkCounter(len(xrange(*jobs))*tries)
#    np.random.seed(42)
#    fwd_matches = np.random.randint(0, t_len, n/2)
#    rev_matches = np.random.randint(0, t_len, n-n/2)
#    dups = None
##    dups = searcher.compile_duplexes_mp(WorkCounter(), fwd_seq, rev_seq, primer, t_len, p_len, fwd_matches, rev_matches)
#    cProfile.run('dups = searcher.compile_duplexes(WorkCounter(), fwd_seq, rev_seq, primer, t_len, p_len, fwd_matches, rev_matches)', 'compile_duplexes.profile')
#    print dups
#    for ji in xrange(*jobs):
#        for _i in xrange(tries):
#            j.append(ji)
#            p.append(timeit.timeit('searcher.compile_duplexes_mp(WorkCounter(), fwd_seq, rev_seq, primer, t_len, p_len, fwd_matches, rev_matches, ji)', 
#                                   'from __main__ import searcher, fwd_seq, rev_seq, t_len, p_len, primer, fwd_matches, rev_matches, ji\n'
#                                   'from WorkCounter import WorkCounter', 
#                                   number=1))
#            c.count()
#            print c.percent()
#    
#    from matplotlib import pyplot as plt
#    from scipy.stats import linregress
#    j = np.fromiter(j, dtype=int)
#    lm = linregress(j, p)
#    plt.plot(j, p, 'b.', j, j*lm[0]+lm[1], 'r-') 
#    plt.show()

    #find vs regexp
    from DegenPrimer.SeqUtils import IUPAC_ambiguous, IUPAC_unambiguous
    def primer2re(primer):
        re_str = ''
        for l in primer.master_sequence:
            if l in IUPAC_unambiguous:
                re_str += l
            else:
                re_str += '[%s]' % (''.join(IUPAC_ambiguous[l]))
        return re.compile(re_str)
    
    import re
    def find_vs_regexp(primer, template):
        find_times = []
        re_times = []
        for _i in xrange(10): 
            #find
            t0 = time()
            searcher = SearchEngine(abort_event)
            p_len,t_len = len(primer),len(template)
            res = searcher._find(WorkCounter(), template, primer, t_len, p_len, 10)
            find_times.append(time()-t0)
            #regexp
            t0 = time()
            fwd_seq = str(template.seq)
            rev_seq = str(template.seq.reverse_complement())
            primer_re = primer2re(primer)
            fwd_res = primer_re.findall(fwd_seq)
            rev_res = primer_re.findall(rev_seq)
            re_times.append(time()-t0)
        print 'find results: f %d, r %d' % (len(res[5]), len(res[6]))
        print 're results:   f %d, r %d' % (len(fwd_res), len(rev_res))
#        print template.seq[res[5][0]:res[5][0]+p_len]
#        print fwd_res[0]
        return find_times, re_times
            
#    print 'template length:   %d' % len(template)
#    print 'primer length:     %d' % len(primer)
#    print 'primer components: %d' % primer.num_components
#    print
#    ft, rt = find_vs_regexp(primer, template)
#    print
#    print 'find mean: %f sec' % (sum(ft)/float(len(ft)))
#    print 're mean: %f sec' % (sum(rt)/float(len(rt)))
    
    print 'Done.'
    
if __name__ == '__main__':
    test()