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

import errno
from time import time
from datetime import timedelta, datetime
from Queue import Empty
from threading import Lock
from Bio import SeqIO
from BioUtils.Tools.Output import OutQueue
from BioUtils.Tools import EchoLogger
from BioUtils.Tools.UMP import UManager
from BioUtils.Tools import WaitingThread

from .PrimerTaskBase import PrimerTaskBase
from .BlastPrimers import BlastPrimers
from .AllSecStructures import AllSecStructures

from .WorkCounter import WorkCounterManager
from .SeqDB import SeqDB
from .iPCR import iPCR
from . import TD_Functions as tdf


class SubroutineManager(UManager): pass
SubroutineManager.register('AllSecStructures', AllSecStructures)
SubroutineManager.register('SeqDB', SeqDB)
SubroutineManager.register('iPCR', iPCR)
SubroutineManager.register('BlastPrimers', BlastPrimers)


class AnalysisTask(PrimerTaskBase):
    '''This class gathers all calculations in a single pipeline with 
    parallelized subroutines for CPU-intensive computations.'''
    
    #file format for saving primers
    _fmt = 'fasta'
    _counter_threshold = 60

    def __init__(self, abort_event):
        PrimerTaskBase.__init__(self, abort_event)
        self._print_lock   = Lock()
        self._out_queue    = OutQueue()
        self._counter_mgr  = WorkCounterManager()
        self._counter_mgr.start()
        self._subroutines  = []
        self._time0        = time()
        self._newline_last = False
    #end def
    
    def __del__(self):
        self._counter_mgr.shutdown() 
        self._clean_subroutines()
    #end def

    def _clean_subroutines(self):
        for p_entry in self._subroutines:
            try: p_entry['manager'].shutdown()
            except: pass
        self._subroutines = []
    #end def


    def _generate_subroutine_entry(self):
        _id = (len(self._subroutines)+1)
        with self._out_queue.set_id(_id):
            manager = SubroutineManager()
            manager.start()
        p_entry  = {'id'       : _id,
                    'manager'  : manager,
                    'process'  : None,
                    'counter'  : self._counter_mgr.WorkCounter()}
        self._subroutines.append(p_entry)
        return p_entry
    #end def
    
    
    def _run_subroutine(self, func, args, p_entry, p_name=None):
        '''Run a func in a thread with args'''
        routine_name = p_name if p_name else func.__name__
        routine_id   = p_entry['id']
        args = [p_entry['counter']]+(list(args) if args else [])
        subroutine   = WaitingThread(self._print_lock, routine_id, target=func, 
                                     name=routine_name, args=args)
        p_entry['process'] = subroutine
        with self._print_lock: 
            self._print(('\nStarting CPU-intensive task #%d:\n  '
                         ' %s\nThis may take awhile...') %
                        (routine_id, routine_name))
        subroutine.start()
    #end def
    
    
    def _run_analysis(self, args):
        analysis_routines = []
        #--------generate components calculate Tms and save the primers--------#
        if not self._prepare_primers(args): return None
        if not self._save_primers(): return None
        #----------------------------------------------------------------------#
        #following computations are CPU intensive, so they need parallelization#
        #check primers for hairpins, dimers and cross-dimers                   #
        self._print('\nSearching for stable secondary structures of the primers...')
        p_entry = self._generate_subroutine_entry()
        all_sec_structures  = p_entry['manager'].AllSecStructures(self._abort_event, args.job_id, self._primers)
        all_sec_structures.find_structures()
        if self.aborted(): return None
        analysis_routines.append(all_sec_structures)
        side_reactions      = all_sec_structures.reactions()
        side_concentrations = all_sec_structures.concentrations()
        with self._print_lock: self._print('Done.')
        self._run_subroutine(all_sec_structures.calculate_equilibrium, None, p_entry,
                             'Calculate relative concentrations of secondary structures at Tm.')
        #----------------------------------------------------------------------#
        #in-silico PCR simulation. This is only available if sequence database #
        #is provided in some form                                              #
        if args.fasta_files or args.sequence_db:
            p_entry = self._generate_subroutine_entry()
            ipcr = p_entry['manager'].iPCR(self._abort_event,
                                           args.max_mismatches,
                                           args.job_id, 
                                           self._primers, 
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
                                                        self._primers, 
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
        #----------------------------------------------------------------------#
        return analysis_routines
    #end def


    def _save_primers(self):
        for primer in self._primers:
            if self.aborted(): return False
            if not primer: continue
            filename = primer.id+'.'+self._fmt
            try:
                SeqIO.write(primer.all_seq_records, filename, self._fmt)
            except Exception, e:
                self._print('\nFailed to write %s primer and it\'s unambiguous components to:\n   %s' % (primer.id, filename))
                print e
                return False 
            self._print('\n%s primer and it\'s unambiguous components were written to:\n   %s' % (primer.id, filename))
        return True
    #end def
    
    
    def _elapsed(self):
        return str(timedelta(seconds=time()-self._time0))
    
    def _join_finished(self):
        processes_alive = False
        for p_entry in self._subroutines:
            if p_entry['process'] == None: continue
            #check if process is alive
            if p_entry['process'].is_alive():
                processes_alive = True
            else:
                p_entry['process'].join()
                p_entry['process'] = None
        return processes_alive
    #end def
    
    def _print_counters(self):
        with self._print_lock:
            counters = ''
            for p_entry in self._subroutines:
                if p_entry['process'] == None: continue
                if p_entry['counter'].last_checked() < self._counter_threshold: continue
                percent = p_entry['counter'].changed_percent()
                if percent is None: continue
                counters += self._fmt_msg(p_entry['id'], percent, 'elapsed: %s\n' % self._elapsed()) 
            if counters: self._print('\n'+counters)
    #end def
    
    def _fmt_msg(self, p_id, percent, msg):
        return '[%d] [%6.2f%%] [%s] %s' % \
                (p_id, percent, 
                 datetime.now().strftime('%d.%m.%Y %H:%M:%S'), 
                 msg)
    
    def _listen_subroutines(self):
        processes_alive = True
        prev_id         = None
        while processes_alive:
            messages = []
            try:
                while True:
                    message = self._out_queue.get(True, 0.1)
                    if message[1] and message[1] != '\n': 
                        messages.append(message)
            except Empty:
                if messages:
                    with self._print_lock:
                        for message in messages:
                            p_id    = message[0]
                            percent = self._subroutines[p_id-1]['counter'].percent()
                            lines   = message[1].strip('\n').split('\n')
                            if len(lines) > 1 or prev_id != p_id: self._print() 
                            for line in lines:
                                if line: self._print(self._fmt_msg(p_id, percent, unicode(line)))
                                else: self._print()
                            if len(lines) > 1: self._print() 
                            prev_id = p_id 
                else: self._print_counters()
            #errors
            except IOError, e:
                if e.errno != errno.EINTR:
                    with self._print_lock:
                        print '\nAnalysisTask._listen_subroutines:'
                        print e
            except Exception, e:
                with self._print_lock:
                    print '\nAnalysisTask._listen_subroutines:'
                    print e
            finally: processes_alive = self._join_finished()
    #end def
    
    def _print_out_queue(self):
        output = []
        while not self.aborted():
            try: output.append(self._out_queue.get_nowait())
            except Empty:
                if output: 
                    self._print('\n'+''.join(message[1] for message in output))
                return
            except Exception, e:
                print e
    #end def
   
    def run(self, args):
        with EchoLogger(args.job_id):
            #reset global state of the module
            self._subroutines = []
            #zero time, used to calculate elapsed time
            self._time0 = time()
            #set PCR parameters
            tdf.PCR_P.set(args.options)
            #run analysis routines and wait for them to finish
            analysis_routines = self._run_analysis(args)
            self._listen_subroutines()
            #if subroutines have been terminated, abort pipeline
            if not analysis_routines or self.aborted():
                self._clean_subroutines()
                with self._print_lock: 
                    self._print('\nAnalysis Task aborted')
                    self._print('\nTotal elapsed time: %s' % self._elapsed())
                return 1 if analysis_routines else -1
            #write reports
            with self._print_lock: self._print('\nWriting analysis reports...')
            for routine in analysis_routines:
                if routine.have_results():
                    routine.write_reports()
                    args.register_reports(routine.reports())
            #print last messages
            self._print_out_queue()
            #terminate managers
            self._clean_subroutines()
            self._print('\nDone. Total elapsed time: %s' % self._elapsed())
        return 1
    #end def
#end class