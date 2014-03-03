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
Created on Jul 25, 2012

@author: Allis Tauri <allista@gmail.com>
'''

#imports
import sys
import errno
import traceback as tb
from time import time, sleep
from datetime import timedelta
from Queue import Empty
from multiprocessing import Queue
from threading import Thread, Lock
from contextlib import contextmanager
import TD_Functions
from PrimerTaskBase import PrimerTaskBase
from BlastPrimers import BlastPrimers
from AllSecStructures import AllSecStructures
from MultiprocessingBase import UManager
from StringTools import print_exception
from SeqDB import SeqDB
from iPCR import iPCR
try:
    from Bio import SeqIO
except ImportError:
    print'The BioPython must be installed in your system.'
    raise
################################################################################ 


#managers for subprocessed classes
class SubroutineManager(UManager): pass
SubroutineManager.register('AllSecStructures', AllSecStructures)
SubroutineManager.register('SeqDB', SeqDB)
SubroutineManager.register('iPCR', iPCR)
SubroutineManager.register('BlastPrimers', BlastPrimers)


class OutQueue(object):
    '''A file-like object which puts text written into it 
    in a cross-process queue'''
    def __init__(self, queue=None):
        if queue == None:
            self._queue = Queue()
        else: self._queue = queue
    #end def
    
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
    def __init__(self, lock, t_id, target=None, name=None, args=(), kwargs=dict()):
        Thread.__init__(self, name=name)
        self._lock   = lock
        self._id     = t_id
        self._target = target
        self._args   = args
        self._kwargs = kwargs
        self.daemon  = True
    #end def
        
    def _print_exception(self):
        with self._lock:
            print '\nError in thread: %s' % self.name
            tb.print_exc()
    #end def
        
    def run(self):
        elapsed_time = 0
        try: elapsed_time = self._target(*self._args, **self._kwargs)
        #Ctrl-C
        #EOF means that target activity was terminated in the Manager process
        except (KeyboardInterrupt, EOFError): return
        #IO code=4 means the same
        except IOError, e:
            if e.errno == errno.EINTR: return
            elif e.errno == errno.EBADMSG:
                with self._lock:
                    print '\nError in thread: %s' % self.name
                    print e.message
                    print '\n*** It seems that an old degen_primer ' \
                    'subprocess is running in the system. Kill it and try again. ***\n'
            else: self._print_exception()
        except Exception: self._print_exception()
        finally:
            #send back a message with elapsed time
            sleep(0.1)
            with self._lock:
                print('\nTask #%d has finished:\n   %s\nElapsed time: %s' % \
                      (self._id, self.name, timedelta(seconds=elapsed_time)))
    #end def
#end class


class AnalysisTask(PrimerTaskBase):
    '''This class gathers all calculations in a single pipeline with 
    parallelized subroutines for CPU-intensive computations.'''
    
    #file format for saving primers
    _fmt = 'fasta'
    
    #a context manager to capture output of print into OutQueue object
    @staticmethod
    @contextmanager
    def capture_to_queue(out_queue=None):
        oldout,olderr = sys.stdout, sys.stderr
        try:
            out = OutQueue(out_queue)
            #need to leave stderr due to the http://bugs.python.org/issue14308
            sys.stdout = out #sys.stderr = out
            yield out
        except Exception, e:
            print_exception(e)
            raise
        finally:
            sys.stdout,sys.stderr = oldout, olderr
    #end def
    

    def __init__(self, abort_event):
        PrimerTaskBase.__init__(self, abort_event)
        self._print_lock  = Lock()
        self._subroutines = [] #list of launched processes with their managers and output queues
    #end def
    
    
    def __del__(self): self._clean_subroutines()


    def _clean_subroutines(self):
        for p_entry in self._subroutines:
            try: p_entry['manager'].shutdown()
            except: pass
        self._subroutines = []
    #end def


    def _generate_subroutine_entry(self):
        with self.capture_to_queue() as out:
            manager = SubroutineManager()
            manager.start()
        e_id     = (len(self._subroutines)+1)
        p_entry  = {'process'  : None,
                    'manager'  : manager,
                    'queue'    : out.queue,
                    'output'   : [],
                    'id'       : e_id}
        self._subroutines.append(p_entry)
        return p_entry
    #end def
    
    
    #subprocess worker factory and launcher
    def _run_subroutine(self, func, args, p_entry, p_name=None):
        '''Run a func in a thread with args provided as dictionary'''
        def worker(_args):
            #remember start time
            time0 = time()
            #run the function
            if _args: func(*_args)
            else: func()
            #return elapsed time
            return (time()-time0)
        #start subprocess
        routine_name = p_name if p_name else func.__name__
        routine_id   = p_entry['id']
        subroutine   = WaitThread(self._print_lock, routine_id, target=worker, 
                                  name=routine_name, args=[args,])
        p_entry['process'] = subroutine
        with self._print_lock: 
            print('\nStarting CPU-intensive task #%d:\n   %s\nThis may take awhile...' % \
                  (routine_id, routine_name))
        subroutine.start()
    #end def
    
    
    def _run_analysis(self, args):
        analysis_routines = []
        #-------------generate components and calculate primers' Tm------------#
        print ('\nCalculating unambiguous components and melting temperatures '
               'of the primers...')
        for primer in args.primers:
            if self._abort_event.is_set(): return None
            primer.generate_components()
            primer.calculate_Tms(self._abort_event)
        print 'Done.'
        #save primers 
        if not self._save_primers(args.primers): return None
        #----------------------------------------------------------------------#
        #following computations are CPU intensive, so they need parallelization#
        #check primers for hairpins, dimers and cross-dimers                   #
        print ('\nSearching for stable secondary structures of the primers...')
        p_entry = self._generate_subroutine_entry()
        all_sec_structures  = p_entry['manager'].AllSecStructures(self._abort_event, args.job_id, args.primers)
        all_sec_structures.find_structures()
        if self._abort_event.is_set(): return None
        analysis_routines.append(all_sec_structures)
        side_reactions      = all_sec_structures.reactions()
        side_concentrations = all_sec_structures.concentrations()
        with self._print_lock: print 'Done.'
        self._run_subroutine(all_sec_structures.calculate_equilibrium, None, p_entry,
                             'Calculate relative concentrations of secondary structures at Tm.')
        #----------------------------------------------------------------------#
        #in silica PCR simulation. This is only available if sequence database #
        #is provided in some form                                              #
        if args.fasta_files or args.sequence_db:
            p_entry = self._generate_subroutine_entry()
            ipcr = p_entry['manager'].iPCR(self._abort_event,
                                           args.max_mismatches,
                                           args.job_id, 
                                           args.primers, 
                                           args.min_amplicon, 
                                           args.max_amplicon, 
                                           args.polymerase, 
                                           args.with_exonuclease, 
                                           args.cycles,
                                           side_reactions, 
                                           side_concentrations,
                                           args.analyse_all_annealings)
            analysis_routines.append(ipcr)
            #connect to sequence database
            seq_files = []
            if args.sequence_db: seq_files.append(args.sequence_db)
            else: seq_files = args.fasta_files
            self._run_subroutine(ipcr.simulate_PCR, 
                                 (seq_files, 
                                  args.use_sequences), p_entry,
                                 'Simulate PCR using provided sequences as DNA templates.')
        #-----------------test for primers specificity by BLAST----------------#
        p_entry = self._generate_subroutine_entry()
        blast_primers = p_entry['manager'].BlastPrimers(self._abort_event,
                                                        args.job_id, 
                                                        args.primers, 
                                                        args.min_amplicon, 
                                                        args.max_amplicon, 
                                                        args.polymerase, 
                                                        args.with_exonuclease, 
                                                        args.cycles,
                                                        side_reactions, 
                                                        side_concentrations,
                                                        include_side_annealings=args.analyse_all_annealings)
        analysis_routines.append(blast_primers)
        #if --do-blast flag was provided, make an actual query
        if args.do_blast:
            #construct Entrez query
            entrez_query = ''
            if args.organisms:
                for organism in args.organisms:
                    if entrez_query: entrez_query += ' OR '
                    entrez_query += organism+'[organism]'
            #do the blast and analyze the results
            self._run_subroutine(blast_primers.blast_and_analyze, 
                                 (entrez_query,),
                                 p_entry,
                                 'Make BLAST query, then simulate PCR using returned alignments.')
        #else, try to load previously saved results and analyze them with current parameters
        elif blast_primers.have_saved_results():
            self._run_subroutine(blast_primers.load_and_analyze, 
                                 None, 
                                 p_entry,
                                 'Simulate PCR using alignments in BLAST results.')
        else: self._print_queue(p_entry['queue'])
        #----------------------------------------------------------------------#
        return analysis_routines
    #end def


    def _save_primers(self, primers):
        for primer in primers:
            if self._abort_event.is_set(): return False
            if not primer: continue
            filename = primer.id+'.'+self._fmt
            try:
                SeqIO.write(primer.all_seq_records, filename, self._fmt)
            except Exception, e:
                print '\nFailed to write %s primer and it\'s unambiguous components to:\n   %s' % (primer.id, filename)
                print_exception(e)
                return False 
            print '\n%s primer and it\'s unambiguous components were written to:\n   %s' % (primer.id, filename)
        return True
    #end def
    
    
    def _print_queue(self, queue):
        output = []
        while not self._abort_event.is_set():
            try: output.append(queue.get(False))
            except Empty:
                if output: print ''.join(message for message in output)
                return
            except Exception, e:
                print_exception(e)
    #end def
                
    
    def _listen_subroutines(self):
        '''While waiting for their termination, get output from subroutines. 
        When some subroutine has finished, print out it's output.'''
        while True: #not self._abort_event.is_set():
            processes_alive = False
            for p_entry in self._subroutines:
                if p_entry['process'] == None: continue
                try:
                    while True: p_entry['output'].append(p_entry['queue'].get_nowait())
                except Empty:
                    if p_entry['process'].is_alive():
                        processes_alive = True
                    else:
                        p_entry['process'].join()
                        p_entry['process'] = None
                    if p_entry['output']:
                        lines = ''.join(unicode(message) for message in p_entry['output']).split('\n')
                        with self._print_lock:
                            for line in lines:
                                if line: print('[%d] %s' % (p_entry['id'], line))
                                else: print ''
                        p_entry['output'] = []
                    continue
                except IOError, e:
                    if e.errno == errno.EINTR:
                        processes_alive = True
                        continue
                    else:
                        with self._print_lock:
                            print '\nAnalysisTask._listen_subroutines:'
                            print_exception(e)
                except Exception, e:
                    with self._print_lock:
                        print '\nAnalysisTask._listen_subroutines:'
                        print_exception(e)
                processes_alive = True
            if not processes_alive: break
            sleep(0.1)
    #end def

    
    def run(self, args):
        #reset global state of the module
        self._subroutines = []
        #zero time, used to calculate elapsed time
        time0 = time()
        #set PCR parameters
        TD_Functions.PCR_P.set(args.options)
        #run analysis routines and wait for them to finish
        analysis_routines = self._run_analysis(args)
        self._listen_subroutines()
        #if subroutines have been terminated, abort pipeline
        if self._abort_event.is_set():
            self._clean_subroutines()
            with self._print_lock: 
                print '\nAnalysisTask aborted.'
                print '\nTotal elapsed time: %s' % timedelta(seconds=time()-time0)
            return False
        #write reports
        with self._print_lock: print '\nWriting analysis reports...'
        for routine in analysis_routines:
            if routine.have_results():
                routine.write_reports()
                args.register_reports(routine.reports())
        #print last messages
        for p_entry in self._subroutines: self._print_queue(p_entry['queue'])
        #terminate managers
        self._clean_subroutines()
        print '\nDone. Total elapsed time: %s' % timedelta(seconds=time()-time0)
        return 1
    #end def
#end class