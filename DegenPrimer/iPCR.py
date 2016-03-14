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

from BioUtils.Tools.Output import simple_timeit
from BioUtils.Tools.UMP import FuncManager, at_manager
from BioUtils.Tools.Multiprocessing import raise_tb_on_error
from BioUtils.Tools.Text import wrap_text, hr

from .PCR_ProductsFinder import PPFManager
from .SearchEngine import mp_better
from .iPCR_Base import iPCR_Base


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
    
    def _find_products_in_templates(self, counter, templates, P_Finder):
        if not templates: return None
        if len(templates) == 1:
            return self._searcher.find(counter,
                                       templates[0], 
                                       self._max_mismatches, 
                                       P_Finder)
        else:
            return self._searcher.batch_find(counter, 
                                             templates,
                                             self._max_mismatches, 
                                             P_Finder, copy_data=True)
    #end def
    
    
    def _find_products_in_db(self, counter, seq_ids, PCR_Sim, P_Finder):
        #sort templates into short, suitable for plain search and long -- for mp
        seq_lengths = dict()
        short_templates = []
        long_templates  = []
        with simple_timeit('PCR Simulation: sorting templates'):
            for t_id in seq_ids:
                template = self._seq_db[t_id]
                if template is None:
                    print 'Sequence %s not found. Something is not right with the sequence database.' % t_id 
                    continue
                t_len = len(template)
                seq_lengths[t_id] = t_len
                if mp_better(t_len): long_templates.append(t_id)
                else: short_templates.append(t_id)
            if self.aborted(): return False
        if len(short_templates) == 1: 
            long_templates += short_templates
            short_templates = []
        #setup work counters
        with simple_timeit('PCR Simulation: counting work to be done'):
            if short_templates and long_templates: 
                counter.set_subwork(2, (sum(seq_lengths[t_id] for t_id in short_templates),
                                        sum(seq_lengths[t_id] for t_id in long_templates)))
                short_counter  = counter[0]
                long_counter   = counter[1]
            else: long_counter = short_counter = counter
            if long_templates: 
                long_counter.set_subwork(len(long_templates), 
                                         [seq_lengths[t_id] for t_id in long_templates])
            if short_templates:
                short_counter.set_subwork(len(short_templates))
        print '\nPCR Simulation: searching for annealing sites in provided sequences...'
        #search short templates in batch
        if short_templates:
            results = self._find_products_in_templates(short_counter, 
                                                       self._seq_db.subview(short_templates), 
                                                       P_Finder)
            if not results or self.aborted(): return False
            for t_id, m_path in results.iteritems():
                PCR_Sim.add_mixture(t_id, m_path)
        #if there're long templates, search sequentially
        for i, t_id in enumerate(long_templates):
            result = self._find_products_in_templates(long_counter[i], [self._seq_db[t_id]], P_Finder)
            if result is None:
                if self.aborted(): return False
                else: continue
            PCR_Sim.add_mixture(*result.items()[0])
        return PCR_Sim.not_empty()
    #end def
    
    
    def _find_products(self, counter, PCR_Sim, P_Finder, seq_files, seq_ids):
        with simple_timeit('PCR Simulation: loading templates'):
            if not seq_files or not self._load_db(seq_files):
                print 'No templates were loaded from: %s' % str(seq_files) 
                return False
            if self.aborted(): return False
        if not seq_ids: seq_ids = self._seq_db.keys()
        else: seq_ids = [str(sid) for sid in seq_ids]
        print 'Number of templates to process: %d' % len(seq_ids)
        return self._find_products_in_db(counter, seq_ids, PCR_Sim, P_Finder)
    #end def
    
    @raise_tb_on_error
    def simulate_PCR(self, counter, seq_files, seq_ids=None):
        counter.set_subwork(2, (10, 1+5*self._include_side_annealings))
        if self._find_products(counter[0], 
                               self._PCR_Simulation, 
                               self._ProductsFinder, 
                               seq_files, seq_ids):
            if self.aborted(): return False
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
        body += hr(' histogram of all possible PCR products ', symbol='=') 
        if len(self._PCR_Simulation.hits()) == 1:
            #all products histogram
            hit = self._PCR_Simulation.hits()[0]
            body += self._PCR_Simulation.per_hit_header(hit)
            body += self._PCR_Simulation.per_hit_histogram(hit)
            body += '\n'
            body += hr(' electrophorogram of all possible PCR products ', symbol='=')
            body += self._PCR_Simulation.per_hit_electrophoresis(hit)
        else:
            products = self._PCR_Simulation.products().items()
            products.sort(key=lambda(d): max(p.quantity for p in d[1].values()), reverse=True)
            #all products histogram
            body += self._PCR_Simulation.all_products_histogram(products)
            body += '\n\n\n'
            #per hit histogram and phoresis
            body += hr(' histograms and electrophorograms of PCR products of each hit ', symbol='=')
            body += self._PCR_Simulation.all_graphs_grouped_by_hit(products)     
        return body
    #end def
#end class