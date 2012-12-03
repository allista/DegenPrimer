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
Created on Jul 25, 2012

@author: Allis Tauri <allista@gmail.com>
'''

#imports
import sys
import errno
from time import time
from datetime import timedelta
from copy import deepcopy
from Queue import Empty
from multiprocessing import Process, Queue
from multiprocessing.managers import BaseManager
from contextlib import contextmanager
from DegenPrimer.OligoFunctions import generate_unambiguous, self_complement
from DegenPrimer import TD_Functions    
from DegenPrimer.TD_Functions import calculate_Tm, source_feature, add_PCR_conditions, format_PCR_conditions
from DegenPrimer.Blast import Blast
from DegenPrimer.SecStructures import AllSecStructures
from DegenPrimer.StringTools import hr, wrap_text, print_exception, time_hr
from DegenPrimer.iPCR import iPCR
try:
    from Bio import SeqIO
except Exception, e:
    print_exception(e)
    raise ImportError('The BioPython must be installed in your system.')
################################################################################ 


#managers for subprocessed classes
class AllSecStructures_Manager(BaseManager): pass
AllSecStructures_Manager.register('AllSecStructures', AllSecStructures)

class iPCR_Manager(BaseManager): pass
iPCR_Manager.register('iPCR', iPCR)

class Blast_Manager(BaseManager): pass
Blast_Manager.register('Blast', Blast)


class OutQueue(object):
    '''A file-like object which puts text written into it 
    in a cross-process queue'''
    def __init__(self, queue=None):
        if queue == None:
            self._queue = Queue()
        else: self._queue = queue
        
    @property
    def queue(self):
        return self._queue
        
    def write(self, text):
        self._queue.put(text)
    
    def flush(self): pass
#end class


#a context manager to capture output of print into OutQueue object
@contextmanager
def capture_to_queue(out_queue=None):
    oldout,olderr = sys.stdout, sys.stderr
    try:
        out = OutQueue(out_queue)
        sys.stdout = sys.stderr = out
        yield out
    except Exception, e:
        print e.message
        raise
    finally:
        sys.stdout,sys.stderr = oldout, olderr
#end def


#subprocess worker factory and launcher
def _subprocess(func, args, queue, p_entry, p_name=None):
    '''
    Run a func in a subrocess with args provided as dictionary.
    Output will be put into a queue q.
    Return Process object. 
    '''
    process_name = p_name if p_name else func.__name__
     
    def worker(_args, _q):
        #capture stdout into a queue
        with capture_to_queue(_q):
            #remember start time
            time0 = time()
            #run the function
            if _args: func(*_args)
            else: func()
            print '\nTask has finished:\n   %s\nElapsed time: %s' % (process_name, timedelta(seconds=time()-time0))
    #end def
    
    #start subprocess
    p = Process(target=worker, args=[args, queue], name=process_name)
    p_entry['process'] = p
    print '\nStarting CPU-intensive task:\n   %s\nThis may take awhile...' % process_name
    p.start()
#end def


class DegenPrimerPipeline(object):
    def __init__(self):
        self._processes  = [] #list of launched processes with their managers and output queues
        self._out_queue  = Queue()
        self._terminated = False #flag which is set when pipeline is terminated
    #end def


    def generate_subprocess_entry(self, manager, queue):
        p_entry = {'process': None,
                   'manager': manager,
                   'queue'  : queue,
                   'output' : []}
        self._processes.append(p_entry)
        return p_entry
    #end def

    def terminate(self):
        self._terminated = True
        for p_entry in self._processes:
            del p_entry['manager']
            del p_entry['queue']
            if p_entry['process'] != None:
                p_entry['process'].terminate()
        self._processes = []
    #end def


    def run(self, args):
        #reset global state of the module
        self._terminated = False
        self._processes  = []
        
        #set job ID and primers list
        primers = args.primers
        job_id  = args.job_id
        
        #set concentrations
        TD_Functions.C_Na   = args.Na
        TD_Functions.C_Mg   = args.Mg
        TD_Functions.C_dNTP = args.dNTP
        TD_Functions.C_DNA  = args.DNA
        TD_Functions.C_Prim = args.Primer
        TD_Functions.C_DMSO = args.DMSO
        TD_Functions.PCR_T  = args.PCR_T
        
        #check if at least one primer is provided
        if not args.sense_primer and not args.antisense_primer:
            print 'At least one primer (sense or antisense) should be provided.'
            return 1
        #test for self-complementarity
        for primer in primers:
            if primer and self_complement(primer[0]):
                print 'Error: %s primer "%s" [%s] is self-complementary.' % \
                        (primer[0].description, primer[0].id, primer[0].seq)
                print 'You should not use self-complementary oligonucleotides as primers.\n'
                return 1
        #save the configuration only after preliminary checks
        args.save_configuration()
        print ''
        #--------------------------------------------------------------------------#
        
        #zero time, used to calculate elapsed time in the end
        time0 = time()
        #generate unambiguous primer sets
        for primer in primers:
            try:
                if primer: primer[1] = generate_unambiguous(primer[0])
            except Exception, e:
                print_exception(e)
                return 1
        #--------------------------------------------------------------------------#
        
        
        #primer lists
        def pfam_primers(pfam):
            '''Return list of unambiguous primers:
            if pfam=0, list of forward primers
            if pfam=1, list of reverse primers'''
            if not primers[pfam]: return [] 
            if not primers[pfam][1]:
                return [primers[pfam][0]]
            else: return primers[pfam][1]
        #end def
        
        def all_primers():
            return pfam_primers(0) + pfam_primers(1)
        #--------------------------------------------------------------------------#
        
        
        #calculate melting temperatures for primers
        for primer in primers:
            if not primer: continue
            #if original primer is not a degenerate
            if not primer[1]: 
                primer_Tm = calculate_Tm(primer[0])
                primer[2]['Tm'] = primer_Tm
                continue
            #else...
            n_primers = len(primer[1])
            mean_Tm, min_Tm, max_Tm = 0, 10000, 0 
            for p in range(n_primers):
                primer_Tm = calculate_Tm(primer[1][p])
                if primer_Tm:
                    primer_Tm = primer_Tm
                    mean_Tm  += primer_Tm
                    if min_Tm >= primer_Tm: min_Tm = primer_Tm
                    if max_Tm <  primer_Tm: max_Tm = primer_Tm
                else: n_primers -= 1
            feature = source_feature(primer[0])
            add_PCR_conditions(feature)
            feature.qualifiers['Tm_min']  = str(min_Tm)
            feature.qualifiers['Tm_max']  = str(max_Tm)
            feature.qualifiers['Tm_mean'] = str(mean_Tm/n_primers)
            primer[2]['Tm_min']  = min_Tm
            primer[2]['Tm_max']  = max_Tm
            primer[2]['Tm_mean'] = mean_Tm/n_primers
        #--------------------------------------------------------------------------#
        
        
        #fasta file with primers
        fmt = 'fasta'
        fmt_filename = job_id+'.'+fmt
        primers_list = list()
        for primer in primers:
            if not primer: continue
            primers_list.append(primer[0])
            primers_list += primer[1]
        SeqIO.write(primers_list, fmt_filename, fmt)
        print '\nUnambiguous primers were written to:\n   ',fmt_filename
        #--------------------------------------------------------------------------#
        
        
        #--following computations are CPU intensive, so they need parallelization--#
        #check primers for hairpins, dimers and cross-dimers
        with capture_to_queue() as out:
            p_entry = self.generate_subprocess_entry(AllSecStructures_Manager(), out.queue)
            p_entry['manager'].start()
            all_sec_structures = p_entry['manager'].AllSecStructures(pfam_primers(0), pfam_primers(1))
        side_reactions       = deepcopy(all_sec_structures.reactions())
        side_concenctrations = deepcopy(all_sec_structures.concentrations())
        _subprocess(all_sec_structures.calculate_equilibrium, None, out.queue, p_entry,
                    'Calculate quantities of secondary structures.')
        #--------------------------------------------------------------------------#
        
        
        #ipcress program file and test for primers specificity by iPCR
        #this is only available for pairs of primers
        ipcr = None
        if args.sense_primer and args.antisense_primer:
            fwd_primers  = pfam_primers(0)
            rev_primers  = pfam_primers(1)
            with capture_to_queue() as out:
                p_entry = self.generate_subprocess_entry(iPCR_Manager(), out.queue)
                p_entry['manager'].start()
                ipcr = p_entry['manager'].iPCR(job_id, 
                        fwd_primers, rev_primers, 
                        args.min_amplicon, args.max_amplicon, 
                        args.with_exonuclease, 
                        side_reactions, 
                        side_concenctrations)
            ipcr.writeProgram()
            #if target sequences are provided and the run_icress flag is set, run iPCR...
            if args.run_ipcress and args.fasta_files:
                #change name of the method to name the subprocess
                _subprocess(ipcr.executeProgram, 
                            (args.fasta_files, 
                             args.max_mismatches), out.queue, p_entry,
                            'Execute iPCRess program and calculate quantities of possible products.')
            else:
                #load previously saved results
                print '\nLoading raw iPCR results from the previous ipcress run.'
                _subprocess(ipcr.load_results, None, out.queue, p_entry,
                            'Calculate quantities of possible PCR products found by iPCRess.')
        #--------------------------------------------------------------------------#
        
        
        #test for primers specificity by BLAST
        with capture_to_queue() as out:
            p_entry = self.generate_subprocess_entry(Blast_Manager(), out.queue)
            p_entry['manager'].start()
            blast = p_entry['manager'].Blast(job_id,
                                            args.min_amplicon, 
                                            args.max_amplicon, 
                                            args.with_exonuclease)
        #if --do-blast command was provided, make an actual query
        if args.do_blast:
            #construct Entrez query
            entrez_query = ''
            if args.organisms:
                for organism in args.organisms:
                    if entrez_query: entrez_query += ' OR '
                    entrez_query += organism+'[organism]'
            #do the blast and analyze results
            def blast_and_analyze():
                blast.blast_primers(all_primers(), entrez_query)
                if blast.have_results():
                    blast.find_PCR_products(side_reactions,
                                            side_concenctrations)
            _subprocess(blast_and_analyze, None, out.queue, p_entry,
                        'Make BLAST query and search for possible PCR products')
        #otherwise, reparse previously saved results with current parameters
        else:
            #load blast results
            blast.load_results()
            if blast.have_results():
                _subprocess(blast.find_PCR_products, 
                            (side_reactions,
                             side_concenctrations), 
                            out.queue, p_entry,
                            'Search for possible PCR products in BLAST results.')
        #--------------------------------------------------------------------------#
        
        
        #collect processes output, write it to stdout, then join all processes
        while not self._terminated:
            processes_alive = False
            for p_entry in self._processes:
                if p_entry['process'] == None:
                    continue
                try:
                    p_entry['output'].append(p_entry['queue'].get(True, 1))
                except Empty:
                    if not p_entry['process'].is_alive():
                        p_entry['process'].join()
                        p_entry['process'] = None
                        print ''.join(message for message in p_entry['output'])
                    continue
                except IOError, e:
                    if e.errno == errno.EINTR:
                        processes_alive = True 
                        continue
                    else:
                        self.terminate()
                        raise
                processes_alive = True
            if not processes_alive: break
        
        
        #if processes have been terminated, terminate the whole pipeline
        if self._terminated: return 1
        #--------------------------------------------------------------------------#
        
        
        #-----------------------now write all the reports--------------------------#
        #write full and short reports
        structures_full_report_filename  = job_id+'-full-report.txt'
        structures_short_report_filename = job_id+'-short-report.txt'
        full_structures_file  = open(structures_full_report_filename, 'w')
        short_structures_file = open(structures_short_report_filename, 'w')
        #write header
        for f in (full_structures_file, short_structures_file):
            f.write(time_hr())
            f.write(wrap_text('For each degenerate primer provided, a set of unambiguous primers is generated. '
                              'For each such set the minimum, maximum and mean melting temperatures are calculated. '
                              'For each primer in each set stable self-dimers and hairpins are predicted. '
                              'For every possible combination of two unambiguous primers cross-dimers are also predicted. '
                              'If an unambiguous primer is provided, it is treated as a set with a single element.\n\n'))
            f.write(hr(' PCR conditions '))
            f.write(format_PCR_conditions()+'\n')
            f.write(hr(' primers and their melting temperatures '))
            #temperatures
            temps  = []
            
            if primers[0] and primers[1]:
                max_id = max(len(primers[0][0].id), len(primers[1][0].id))
            else: max_id = 0
            for primer in primers:
                if not primer: continue
                spacer = max_id-len(primer[0].id)
                if not primer[1]: #not degenerate
                    f.write('%s:%s  %02db  %s\n' % (primer[0].id, 
                                                    ' '*spacer, 
                                                    len(primer[0].seq), 
                                                    str(primer[0].seq)) +\
                            '   Tm:      %.1f C\n' % primer[2]['Tm'])
                    temps.append(primer[2]['Tm'])
                else:
                    f.write('%s:%s  %02db  %s\n' % (primer[0].id, 
                                                    ' '*spacer, 
                                                    len(primer[0].seq), 
                                                    str(primer[0].seq)) +\
                            '   Tm max:  %.1f C\n' % primer[2]['Tm_max'] +\
                            '   Tm mean: %.1f C\n' % primer[2]['Tm_mean'] +\
                            '   Tm min:  %.1f C\n' % primer[2]['Tm_min'])
                    temps.append(primer[2]['Tm_min'])
            #warning
            if len(temps) == 2:
                if abs(temps[0]-temps[1]) >= 5:
                    f.write('\nWarning: lowest melting temperatures of sense and antisense primes \n'
                            '         differ more then by 5C\n')
            f.write('\n\n')
        #write secondary structures information
        full_structures_file.write(all_sec_structures.print_structures())
        short_structures_file.write(all_sec_structures.print_structures_short())
        full_structures_file.close()
        short_structures_file.close()
        print '\nFull report with all secondary structures was written to:\n   ',structures_full_report_filename
        print '\nShort report with a summary of secondary structures was written to:\n   ',structures_short_report_filename
        args.register_report('Tm and secondary structures', structures_short_report_filename)
        #--------------------------------------------------------------------------#
    
    
        #write iPCR reports if they are available
        if ipcr != None and ipcr.have_results():
            ipcr.write_PCR_report()
            for report in ipcr.reports():
                if report:
                    args.register_report(**report)
        #--------------------------------------------------------------------------#
        
        
        #write BLAST reports
        if blast.have_results():
            blast.write_hits_report()
            blast.write_PCR_report()
            for report in blast.reports():
                if report:
                    args.register_report(**report)
        #--------------------------------------------------------------------------#
        
        #terminate managers and delete queues
        self.terminate()
        
        print '\nDone. Total elapsed time: %s' % timedelta(seconds=time()-time0)
        return 0
    #end def