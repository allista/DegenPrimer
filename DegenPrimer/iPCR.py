#!/usr/bin/python
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
Created on Feb 27, 2013

@author: Allis Tauri <allista@gmail.com>
'''

from StringTools import wrap_text, hr
from iPCR_Base import iPCR_Base
from SearchEngine import mp_better
from UMP import FuncManager, at_manager
from PCR_ProductsFinder import PPFManager


#manager to isolate forking point in a low-memory process
SearchManager = FuncManager('SearchManager', 
                            (at_manager(PPFManager, 'find'), 
                             at_manager(PPFManager, 'batch_find')))

class iPCR(iPCR_Base):
    '''
    Using PCR_Simulation and SeqDB classes runs PCR simulation with given 
    primers and report results in human readable form in a text file.
    '''
    
    def __init__(self, *args, **kwargs):
        iPCR_Base.__init__(self, *args, **kwargs)
        self._PCR_Simulation = self._new_PCR_Simulation()
        self._ProductsFinder = self._new_PCR_ProductsFinder()
        self._searcher       = SearchManager()
        self._searcher.start()
    #end def
    
    
    def _find_products_in_templates(self, counter, t_ids, products_finder):
        templates   = self._seq_db.get_seqs(t_ids)
        if len(t_ids) < 2:
            return self._searcher.find(counter, self._seq_names[t_ids[0]],
                                       templates[0][2], 
                                       self._max_mismatches, 
                                       products_finder)
        else:
            return self._searcher.batch_find(counter, self._seq_names, 
                                             templates,  
                                             self._max_mismatches, 
                                             products_finder)
    #end def
    
    
    def _find_products_in_db(self, counter, PCR_Sim, P_Finder, seqs_info):
        print '\nPSR Simulation: searching for annealing sites in provided sequences...'
        #sort templates into short, suitable for plain search and long -- for mp
        short_templates = []
        long_templates  = []
        for t_id, t_len in seqs_info.iteritems():
            if mp_better(t_len):
                long_templates.append(t_id)
            else: short_templates.append(t_id)
        #setup work counters
        if short_templates and long_templates: 
            counter.set_subwork(2, (sum(seqs_info[t_id] for t_id in short_templates),
                                    sum(seqs_info[t_id] for t_id in long_templates)))
            short_counter  = counter[0]
            long_counter   = counter[1]
        else: long_counter = short_counter = counter
        if long_templates: long_counter.set_subwork(len(long_templates), [seqs_info[t_id] for t_id in long_templates])
        #search short templates in batch
        if short_templates:
            short_counter.set_subwork(len(short_templates))
            results = self._find_products_in_templates(short_counter, 
                                                       short_templates, 
                                                       P_Finder)
            if not results or self.aborted(): return False
            for t_name, m_path in results.iteritems():
                PCR_Sim.add_mixture(t_name, m_path)
        #if there're long templates, search sequentially
        for i, t_id in enumerate(long_templates):
            result = self._find_products_in_templates(long_counter[i], (t_id,), P_Finder)
            if result is None:
                if self.aborted(): return False
                else: continue
            PCR_Sim.add_mixture(self._seq_names[t_id], result)
        return PCR_Sim.not_empty()
    #end def
    
    
    def _find_products(self, counter, PCR_Sim, P_Finder, seq_files, seq_ids=None):
        #try to connect to a database
        if not seq_files or not self._try_connect_db(seq_files): return False
        #get names and legths of sequences
        self._seq_names = self._get_names(seq_ids)
        if not self._seq_names:
            self._seq_db.close() 
            return False
        seq_lengths = self._seq_db.get_lengths(seq_ids)
        #find primer annealing sites
        result = self._find_products_in_db(counter, PCR_Sim, P_Finder, seq_lengths)
        self._seq_db.close()
        return result
    #end def
    
    
    def simulate_PCR(self, counter, seq_files, seq_ids=None):
        counter.set_subwork(2, (10, 1+5*self._include_side_annealings))
        counter[0]
        if self._find_products(counter[0], 
                               self._PCR_Simulation, 
                               self._ProductsFinder, 
                               seq_files, seq_ids):
            self._PCR_Simulation.run(counter[1])
            self._have_results = bool(self._PCR_Simulation)
            return self._have_results
        return False
    #end def
    
    
    def _format_header(self):
        header =  wrap_text('All possible PCR products are ' 
                            'filtered by amplicon size.\n'
                            'If --no-exonuclease option was ' 
                            'provided, products which may be formed by '
                            'primers with ' 
                            "mismatches on 3'-end are ignored.\n"
                            'Quantities of the remaining products ' 
                            'are estimated using equilibrium equations '
                            'and current PCR parameters.\n\n')
        header += iPCR_Base._format_header(self)
        header += self._PCR_Simulation.format_report_header()
        return header
    #end def
    
    
    def _format_report_body(self):
        body  = ''
        body += self._PCR_Simulation.format_quantity_explanation()
        #PCR products by hit 
        if len(self._PCR_Simulation.hits()) == 1:
            #all products histogram
            hit = self._PCR_Simulation.hits()[0]
            body += hr(' histogram of all possible PCR products ', symbol='=')
            body += self._PCR_Simulation.per_hit_header(hit)
            body += self._PCR_Simulation.per_hit_histogram(hit)
            body += '\n'
            body += hr(' electrophorogram of all possible PCR products ', symbol='=')
            body += self._PCR_Simulation.per_hit_electrophoresis(hit)
        else:
            #all products histogram
            body += hr(' histogram of all possible PCR products ', symbol='=')
            body += self._PCR_Simulation.all_products_histogram()
            body += '\n\n\n'
            #per hit histogram and phoresis
            body += hr(' histograms and electrophorograms of PCR products of each hit ', symbol='=')
            body += self._PCR_Simulation.all_graphs_grouped_by_hit()     
        return body
    #end def
#end class



#tests
if __name__ == '__main__':
    import sys, time
    from multiprocessing import Manager
    from Primer import Primer
    from SeqUtils import load_sequence
    from WorkCounter import WorkCounterManager
    from WaitingThread import WaitingThread
    from threading import Lock
    import TD_Functions as tdf
    
    mgr = Manager()
    abort_event = mgr.Event()
    
    tdf.PCR_P.PCR_T = 53
    tdf.PCR_P.Mg    = 3e-3
    tdf.PCR_P.dNTP  = 300e-6
    tdf.PCR_P.DNA   = 1e-10
    fwd_primer = Primer(load_sequence('ATATTCTACRACGGCTATCC', 'fwd_test', 'fwd_test'), 0.43e-6, True)
    rev_primer = Primer(load_sequence('GAASGCRAAKATYGGGAAC', 'rev_test', 'rev_test'), 0.43e-6, True)
    ipcr = iPCR(abort_event,
                max_mismatches=6,
                job_id='test-job', 
                primers=[fwd_primer,
                         rev_primer], 
                min_amplicon=50, 
                max_amplicon=2000, 
                polymerase=40000, 
                with_exonuclease=False, 
                num_cycles=30,
                side_reactions=None, 
                side_concentrations=None,
                include_side_annealings=True)
    cmgr = WorkCounterManager()
    cmgr.start()
    counter = cmgr.WorkCounter()
    
    plock = Lock()
        
    job = WaitingThread(plock, 1, target=ipcr.simulate_PCR, 
                        name='simulate_PCR',
                        args=(counter, 
                              ('../ThGa.fa', #single sequence
#                               '../Ch5_gnm.fa', '../ThDS1.fa', '../ThES1.fa', #long sequences 
                               '../ThDS1-FC.fa', '../ThDS1-850b-product.fa', #short sequences
                              ),))
    job.start(); print ''
    while job.is_alive():
        if counter.changed_percent():
            with plock: print counter
        time.sleep(1)
    job.join()
    with plock: print counter
    
#    ipcr.write_products_report()
    ipcr.write_report()
    
############################## statistics ######################################
#TD_Functions.PCR_P.PCR_T = 53
#TD_Functions.PCR_P.Mg    = 3e-3
#TD_Functions.PCR_P.dNTP  = 300e-6
#TD_Functions.PCR_P.DNA   = 1e-10
#(abort_event,
#max_mismatches=7,
#job_id='test-job', 
#primers=[#fwd_primer,
#         rev_primer], 
#min_amplicon=50, 
#max_amplicon=2000, 
#polymerase=40000, 
#with_exonuclease=False, 
#num_cycles=30,
#side_reactions=None, 
#side_concentrations=None,
#include_side_annealings=False)

#Results use: 44.114922Mb
#PCR Simulation: searching for possible PCR products in Thermococcus gammatolerans EJ3...
#Results use: 29.034309Mb
#PCR Simulation: searching for possible PCR products in TBCH5...
#Results use: 41.029671Mb
#PCR Simulation: searching for possible PCR products in DS1 complete temp...
#Results use: 21.963631Mb
#PCR Simulation: searching for possible PCR products in Thermococcus sp. ES1...