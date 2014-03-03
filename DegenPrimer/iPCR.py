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

from StringTools import wrap_text, time_hr, hr
from iPCR_Base import iPCR_Base

class iPCR(iPCR_Base):
    '''
    Using PCR_Simulation and SeqDB classes runs PCR simulation with given 
    primers and report results in human readable form in a text file.
    '''
    
    def __init__(self, *args, **kwargs):
        iPCR_Base.__init__(self, *args, **kwargs)
        self._PCR_Simulation = self._new_PCR_Simulation()
        self._ProductsFinder = self._new_PCR_ProductsFinder()
    #end def
    
    
    def simulate_PCR(self, seq_files, seq_ids=None):
        if self._find_products(self._PCR_Simulation, self._ProductsFinder, seq_files, seq_ids):
            self._PCR_Simulation.run()
            self._have_results = bool(self._PCR_Simulation)
            return self._have_results
        return False
    #end def
    
    
    def write_products_report(self):
        if not self._have_results: return
        #open report file
        ipcr_products = self._open_report('iPCR products', self._PCR_products_filename)
        ipcr_products.write(time_hr())
        if self._PCR_Simulation:
            ipcr_products.write(self._PCR_Simulation.format_products_report())
        else: ipcr_products.write(hr(' No PCR products have been found ', symbol='!')) 
        ipcr_products.close()
        print '\nThe list of PCR products was written to:\n   %s' % self._PCR_products_filename
        self._add_report('iPCR products', self._PCR_products_filename)
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
        body = ''
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
    from multiprocessing import Manager
    from Primer import Primer, load_sequence
    import TD_Functions
    
    mgr = Manager()
    abort_event = mgr.Event()
    
    TD_Functions.PCR_P.PCR_T = 53
    TD_Functions.PCR_P.Mg    = 3e-3
    TD_Functions.PCR_P.dNTP  = 300e-6
    TD_Functions.PCR_P.DNA   = 1e-10
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
    
    ipcr.simulate_PCR(('../ThGa.fa', #single sequence
                    '../Ch5_gnm.fa', '../ThDS1.fa', '../ThES1.fa', #long sequences 
                    '../ThDS1-FC.fa', '../ThDS1-850b-product.fa', #short sequences
                       ))
    
    ipcr.write_products_report()
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