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
import os
import sys
import errno
import signal
from time import time, sleep
from datetime import timedelta
from Queue import Empty
from multiprocessing import Queue
from multiprocessing.managers import BaseManager
from threading import Thread
from contextlib import contextmanager
from DegenPrimer.OligoFunctions import generate_unambiguous, self_complement
from DegenPrimer import TD_Functions    
from DegenPrimer.TD_Functions import calculate_Tm, source_feature, add_PCR_conditions, format_PCR_conditions
from DegenPrimer.Blast import Blast
from DegenPrimer.SecStructures import AllSecStructures
from DegenPrimer.StringTools import hr, wrap_text, time_hr, print_exception
from DegenPrimer.iPCR import iPCR
try:
    from Bio import SeqIO
except ImportError:
    print'The BioPython must be installed in your system.'
    raise
################################################################################ 


#managers for subprocessed classes
class ClassManager(BaseManager):
    def getpid(self):
        if '_process' in self.__dict__:
            return self._process.pid
        else: return None
#end class
ClassManager.register('AllSecStructures', AllSecStructures)
ClassManager.register('iPCR', iPCR)
ClassManager.register('Blast', Blast)


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


class WaitThread(Thread):
    '''Thread that specially handle EOFError and IOError 4, interpreting them 
    as a signal that the target activity (which is supposed to run in another 
    process, e.g. in a Manager) has been terminated.'''
    def __init__(self, target=None, name=None, args=(), kwargs=dict()):
        Thread.__init__(self, name=name)
        self._target = target
        self._args   = args
        self._kwargs = kwargs
        
    def run(self):
        try:
            print 'start thread'
            self._target(*self._args, **self._kwargs)
            print 'thread finished'
        #Ctrl-C
        except KeyboardInterrupt: return
        #EOF means that target activity was terminated in the Manager process
        except EOFError: return
        #IO code=4 means the same
        except IOError, e:
            if e.errno == errno.EINTR:
                return
            else: raise
    #end def
#end class


#a context manager to capture output of print into OutQueue object
@contextmanager
def capture_to_queue(out_queue=None):
    oldout,olderr = sys.stdout, sys.stderr
    try:
        out = OutQueue(out_queue)
        #need to left stderr out due to the http://bugs.python.org/issue14308
        sys.stdout = out #sys.stderr = out
        yield out
    except Exception, e:
        print_exception(e)
        raise
    finally:
        sys.stdout,sys.stderr = oldout, olderr
#end def


#subprocess worker factory and launcher
def _subroutine(func, args, queue, p_entry, p_name=None):
    '''
    Run a func in a subrocess with args provided as dictionary.
    Output will be put into a queue q.
    Return Process object. 
    '''
    routine_name = p_name if p_name else func.__name__
     
    def worker(_args, _q):
        #remember start time
        time0 = time()
        #run the function
        if _args is None: func()
        else: func(*_args)
        _q.put('\nTask has finished:\n   %s\nElapsed time: %s' % (routine_name, timedelta(seconds=time()-time0)))
    #end def
    
    #start subprocess
    subroutine = WaitThread(target=worker, args=[args, queue], name=routine_name)
    subroutine.daemon = True
    p_entry['process'] = subroutine
    print '\nStarting CPU-intensive task:\n   %s\nThis may take awhile...' % routine_name
    subroutine.start()
#end def


class DegenPrimerPipeline(object):
    '''This class gathers all calculations in a single pipeline with 
    parallelized subroutines for CPU-intensive computations.'''
    
    #file format for saving primers
    _fmt = 'fasta'
    

    def __init__(self):
        self._subroutines = [] #list of launched processes with their managers and output queues
        self._terminated  = False #flag which is set when pipeline is terminated
    #end def

    def __del__(self):
        self.terminate()


    def terminate(self):
        if self._terminated: return
        self._terminated = True
        #try to shutdown managers
        for p_entry in self._subroutines:
            try: 
                p_entry['manager'].shutdown()
                if p_entry['children']:
                    for child in p_entry['children']:
                        os.kill(child, signal.SIGTERM)
                        sleep(0.1) #control shot
                        os.kill(child, signal.SIGKILL)
            except: pass
        self._subroutines = []
    #end def


    def _generate_subroutine_entry(self, queue):
        p_entry = {'process'  : None,
                   'children' : [],
                   'manager'  : ClassManager(),
                   'queue'    : queue,
                   'output'   : []}
        self._subroutines.append(p_entry)
        return p_entry
    #end def


    def run(self, args):
        #reset global state of the module
        self._terminated  = False
        self._subroutines = []
        
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
            return False
        #test for self-complementarity
        for primer in primers:
            if primer and self_complement(primer[0]):
                print 'Error: %s primer "%s" [%s] is self-complementary.' % \
                        (primer[0].description, primer[0].id, primer[0].seq)
                print 'You should not use self-complementary oligonucleotides as primers.\n'
                return False
        #save the configuration only after preliminary checks
        args.save_configuration()
        print ''
        #----------------------------------------------------------------------#
        
        #zero time, used to calculate elapsed time in the end
        time0 = time()
        
        #generate unambiguous primer sets and calculate melting temperatures
        if not self._process_degenerate_primers(primers): return False
        #save primers 
        if not self._save_primers(primers, job_id): return False
        
        #primer lists 
        fwd_primers = self._pfam_primers(primers, 0)
        rev_primers = self._pfam_primers(primers, 1)
        all_primers = fwd_primers + rev_primers
        #----------------------------------------------------------------------#
    
        #following computations are CPU intensive, so they need parallelization#
        #check primers for hairpins, dimers and cross-dimers
        with capture_to_queue() as out:
            p_entry = self._generate_subroutine_entry(out.queue)
            p_entry['manager'].start()
        all_sec_structures  = p_entry['manager'].AllSecStructures(fwd_primers, rev_primers)
        side_reactions      = all_sec_structures.reactions()
        side_concentrations = all_sec_structures.concentrations()
        _subroutine(all_sec_structures.calculate_equilibrium, None, out.queue, p_entry,
                    'Calculate quantities of secondary structures.')
        #----------------------------------------------------------------------#
        
        #ipcress program file and test for primers specificity by iPCR
        #this is only available for pairs of primers
        ipcr = None
        if args.sense_primer and args.antisense_primer:
            with capture_to_queue() as out:
                p_entry = self._generate_subroutine_entry(out.queue)
                p_entry['manager'].start()
            ipcr = p_entry['manager'].iPCR(job_id, 
                                           fwd_primers, rev_primers, 
                                           args.min_amplicon, args.max_amplicon, 
                                           args.with_exonuclease, 
                                           side_reactions, 
                                           side_concentrations)
            #if target sequences are provided and the run_icress flag is set, run iPCR...
            if args.run_ipcress and args.fasta_files:
                ipcr.write_program()
                _subroutine(ipcr.execute_and_analyze, 
                            (args.fasta_files, 
                             args.max_mismatches), out.queue, p_entry,
                            'Execute iPCRess program and calculate quantities of possible products.')
                #get ipcress pid
                sleep(0.1) #wait some time for process to fork
                ipcress_pid = ipcr.ipcress_pid()
                if ipcress_pid != None:
                    p_entry['children'].append(ipcress_pid)
            #else, try to load previously saved results and analyze them with current parameters
            elif ipcr.load_results():
                print '\nFound raw iPCR results from the previous ipcress run.'
                _subroutine(ipcr.calculate_quantities, None, out.queue, p_entry,
                            'Calculate quantities of possible PCR products found by iPCRess.')
            else: self._print_queue(p_entry['queue'])
        #----------------------------------------------------------------------#
        
        #test for primers specificity by BLAST
        with capture_to_queue() as out:
            p_entry = self._generate_subroutine_entry(out.queue)
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
            _subroutine(blast.blast_and_analyze_primers, 
                        (all_primers(),
                         entrez_query,
                         side_reactions,
                         side_concentrations),
                        out.queue, p_entry,
                        'Make BLAST query and search for possible PCR products.')
        #else, try to load previously saved results and analyze them with current parameters
        elif blast.load_results():
            print '\nFound saved BLAST results.'
            _subroutine(blast.find_PCR_products, 
                        (side_reactions,
                         side_concentrations), 
                        out.queue, p_entry,
                        'Search for possible PCR products in BLAST results.')
        else: self._print_queue(p_entry['queue'])
        #----------------------------------------------------------------------#
        
        #collect subroutines output, write it to stdout, wait for them to finish
        try: self._listen_subroutines()
        except: self.terminate()
        #if subroutines have been terminated, abort pipeline
        if self._terminated:
            print '\nAborted.' 
            return False
        #----------------------------------------------------------------------#
        
        
        #---------------------now write all the reports------------------------#
        #write full and short reports
        structures_full_report_filename  = job_id+'-full-report.txt'
        structures_short_report_filename = job_id+'-short-report.txt'
        full_structures_file  = open(structures_full_report_filename, 'w')
        short_structures_file = open(structures_short_report_filename, 'w')
        #write header
        for f in (full_structures_file, short_structures_file):
            f.write(self._format_primers_report_header(primers))
        #write secondary structures information
        full_structures_file.write(all_sec_structures.print_structures())
        short_structures_file.write(all_sec_structures.print_structures_short())
        full_structures_file.close()
        short_structures_file.close()
        print '\nFull report with all secondary structures was written to:\n   ',structures_full_report_filename
        print '\nShort report with a summary of secondary structures was written to:\n   ',structures_short_report_filename
        args.register_report('Tm and secondary structures', structures_short_report_filename)
    
        #write iPCR reports if they are available
        if ipcr != None and ipcr.have_results():
            ipcr.write_PCR_report()
            for report in ipcr.reports():
                if report:
                    args.register_report(**report)
        
        #write BLAST reports
        if blast.have_results():
            blast.write_hits_report()
            blast.write_PCR_report()
            for report in blast.reports():
                if report:
                    args.register_report(**report)
                    
        #print last queued messages
        sleep(0.1)
        for p_entry in self._subroutines:
            self._print_queue(p_entry['queue'])
        #----------------------------------------------------------------------#
        
        #terminate managers and delete queues
        self.terminate()
        
        print '\nDone. Total elapsed time: %s' % timedelta(seconds=time()-time0)
        return True
    #end def


    @classmethod
    def _pfam_primers(cls, primers, pfam):
        '''Return a list of unambiguous primers:
        if pfam=0, list of forward primers
        if pfam=1, list of reverse primers'''
        if not primers[pfam]: return [] 
        if not primers[pfam][1]:
            return [primers[pfam][0]]
        else: return primers[pfam][1]
    #end def
        
    
    @classmethod
    def _process_degenerate_primers(cls, primers):
        '''Disambiguate primers and calculate their melting temperatures'''
        for primer in primers:
            if not primer: continue
            #generate unambiguous primers 
            try:
                primer[1] = generate_unambiguous(primer[0])
            except ValueError, e:
                print e.message
                return False
            #calculate melting temperatures
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
        return True
    #end def
    
    
    @classmethod
    def _save_primers(cls, primers, job_id):
        fmt_filename = job_id+'.'+cls._fmt
        primers_list = list()
        for primer in primers:
            if not primer: continue
            primers_list.append(primer[0])
            primers_list += primer[1]
        try:
            SeqIO.write(primers_list, fmt_filename, cls._fmt)
        except: return False 
        print '\nUnambiguous primers were written to:\n   ',fmt_filename
        return True
    #end def
    
    
    @classmethod
    def _print_queue(cls, queue):
        output = []
        while True:
            try: output.append(queue.get(False))
            except Empty:
                print ''.join(message for message in output)
                return
    #end def
                
    
    def _listen_subroutines(self):
        '''While waiting for their termination, get output from subroutines. 
        When some subroutine has finished, print out it's output.'''
        while not self._terminated:
            processes_alive = False
            for p_entry in self._subroutines:
                if p_entry['process'] == None:
                    continue
                try:
                    p_entry['output'].append(p_entry['queue'].get(False))
                except Empty:
                    if p_entry['process'].is_alive():
                        processes_alive = True
                    else:
                        p_entry['process'].join()
                        p_entry['process'] = None
                        print ''.join(message for message in p_entry['output'])
                    continue
                except IOError, e:
                    if e.errno == errno.EINTR:
                        processes_alive = True 
                        continue
                    else:
                        print 'Unhandled IOError:', e.message
                        raise
                except Exception, e:
                    print 'Unhandled Exception:'
                    print_exception(e)
                    raise
                processes_alive = True
            if not processes_alive: break
            sleep(0.1)
    #end def
    
    
    @classmethod
    def _format_primers_report_header(cls, primers):
        header_string  = ''
        header_string += time_hr()
        header_string += wrap_text('For each degenerate primer provided, a set '
                                   'of unambiguous primers is generated. '
                                   'For each such set the minimum, maximum and '
                                   'mean melting temperatures are calculated. '
                                   'For each primer in each set stable self-'
                                   'dimers and hairpins are predicted. '
                                   'For every possible combination of two '
                                   'unambiguous primers cross-dimers are also '
                                   'predicted. If an unambiguous primer is '
                                   'provided, it is treated as a set with a '
                                   'single element.\n\n')
        header_string += hr(' PCR conditions ')
        header_string += format_PCR_conditions()+'\n'
        header_string += hr(' primers and their melting temperatures ')
        #print melting temperatures
        temps  = []
        if primers[0] and primers[1]:
            max_id = max(len(primers[0][0].id), len(primers[1][0].id))
        else: max_id = 0
        for primer in primers:
            if not primer: continue
            spacer = max_id-len(primer[0].id)
            if not primer[1]: #not degenerate
                header_string += '%s:%s  %02db  %s\n' % (primer[0].id, 
                                                         ' '*spacer, 
                                                         len(primer[0].seq), 
                                                         str(primer[0].seq))
                header_string += '   Tm:      %.1f C\n' % primer[2]['Tm']
                temps.append(primer[2]['Tm'])
            else:
                header_string += '%s:%s  %02db  %s\n' % (primer[0].id, 
                                                         ' '*spacer, 
                                                         len(primer[0].seq), 
                                                         str(primer[0].seq))
                header_string += '   Tm max:  %.1f C\n' % primer[2]['Tm_max']
                header_string += '   Tm mean: %.1f C\n' % primer[2]['Tm_mean']
                header_string += '   Tm min:  %.1f C\n' % primer[2]['Tm_min']
                temps.append(primer[2]['Tm_min'])
        #warning
        if len(temps) == 2:
            if abs(temps[0]-temps[1]) >= 5:
                header_string += '\nWarning: lowest melting temperatures of sense and antisense primes \n'
                header_string += '         differ more then by 5C\n'
        header_string += '\n\n'
        return header_string
#end class