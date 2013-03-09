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
# indicator_gddccontrol is distributed in the hope that it will be useful, but
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

from StringTools import wrap_text, time_hr, hr, print_exception
from iPCR_Interface import iPCR_Interface
from SeqDB import SeqDB

class iPCR(iPCR_Interface):
    '''
    Using PCR_Simulation and SeqDB classes runs PCR simulation with given 
    primers and report results in human readable form in a text file.
    '''

    def __init__(self, job_id, *args, **kwargs):
        iPCR_Interface.__init__(self, *args, **kwargs)
        #ipcr parameters
        self._seq_db         = SeqDB()
        self._db_filename    = None
        self._max_mismatches = None
        #report files
        self._PCR_report_filename   = job_id+'-iPCR-report.txt'
        self._PCR_products_filename = job_id+'-iPCR-products.txt'
        #results
        self._PCR_Simulation = self._PCR_Simulation_factory()
    #end def
    
    
    def __del__(self):
        self._seq_db.abort_search()
        self._seq_db.close()
    #end def
    
    
    def find_possible_products(self, seq_files, max_mismatches, seq_ids=None):
        if not seq_files: return False
        self._max_mismatches = max_mismatches
        #create database from given files
        try:
            if len(seq_files) == 1 and seq_files[0].endswith('.db'):
                if os.path.isfile(seq_files[0]):
                    self._seq_db.connect(seq_files[0])
                else: print '\niPCR: no such file: %s' % seq_files[0]
            elif not self._seq_db.create_db_from_files(':memory:', seq_files): 
                print '\niPCR: unable to create sequence database from given files.'
                return False
        except Exception, e:
            print_exception(e)
            return False
        #find annealing sites
        seq_names      = self._seq_db.get_names(seq_ids)
        fwd_primer_ann = self._seq_db.find_in_db(self._fwd_primer, 
                                                 max_mismatches, seq_ids)
        rev_primer_ann = self._seq_db.find_in_db(self._rev_primer, 
                                                 max_mismatches, seq_ids)
        self._seq_db.close()
        #find possible products
        if not seq_ids: seq_ids = seq_names.keys()
        self._possible_products = False
        for seq_id in seq_ids:
            fwd_annealings = []
            rev_annealings = []
            if fwd_primer_ann and fwd_primer_ann[seq_id]:
                fwd_annealings.extend(fwd_primer_ann[seq_id][0])
                rev_annealings.extend(fwd_primer_ann[seq_id][1])
            if rev_primer_ann and rev_primer_ann[seq_id]:
                fwd_annealings.extend(rev_primer_ann[seq_id][0])
                rev_annealings.extend(rev_primer_ann[seq_id][1])
            self._possible_products |= self._PCR_Simulation.add_annealings(seq_names[seq_id], 
                                                                           fwd_annealings, 
                                                                           rev_annealings)
        if not self._possible_products:
            print '\nNone of the possible PCR products satisfy given reaction ' \
                  'parameters.'
        return self._possible_products
    #end def
    
    
    def find_and_analyze(self, seq_files, max_mismatches, seq_ids=None):
        return self.find_possible_products(seq_files, max_mismatches) \
            and self.simulate_PCR()
    #end def
    
    
    def simulate_PCR(self):
        if not self._possible_products: return False
        #compute PCR products quantities
        self._PCR_Simulation.run()
        if self._PCR_Simulation:
            self._have_results = True
        return self._have_results
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
        print '\nThe list of PCR products was written to:\n   ',self._PCR_products_filename
        self._add_report('iPCR products', self._PCR_products_filename)
    #end def
    
    
    def write_PCR_report(self):
        if not self._have_results: return
        #open report file
        ipcr_report = self._open_report('iPCR', self._PCR_report_filename)
        #header
        ipcr_report.write(time_hr())
        ipcr_report.write(wrap_text('All possible PCR products are ' 
                                    'filtered by amplicon size.\n'
                                    'If --no-exonuclease option was ' 
                                    'provided, products which may be formed by '
                                    'primers with ' 
                                    "mismatches on 3'-end are ignored.\n"
                                    'Quantities of the remaining products ' 
                                    'are estimated using equilibrium equations '
                                    'and current PCR parameters.\n\n'))
        #filter parameters
        ipcr_report.write(self._PCR_Simulation.format_report_header())
        if self._max_mismatches != None:
            ipcr_report.write('Number of mismatches allowed: %d\n\n' % self._max_mismatches)
        #if no PCR products have been found
        if not self._PCR_Simulation:
            ipcr_report.write(hr(' No PCR products have been found ', symbol='!'))
            ipcr_report.close()
            return
        #else...
        ipcr_report.write(self._PCR_Simulation.format_quantity_explanation())
        #PCR products by hit 
        if len(self._PCR_Simulation.hits()) == 1:
            #all products histogram
            hit = self._PCR_Simulation.hits()[0]
            ipcr_report.write(hr(' histogram of all possible PCR products ', symbol='='))
            ipcr_report.write(self._PCR_Simulation.per_hit_header(hit))
            ipcr_report.write(self._PCR_Simulation.per_hit_histogram(hit))
            ipcr_report.write('\n')
            ipcr_report.write(hr(' electrophorogram of all possible PCR products ', symbol='='))
            ipcr_report.write(self._PCR_Simulation.per_hit_electrophoresis(hit))
        else:
            #all products histogram
            ipcr_report.write(hr(' histogram of all possible PCR products ', symbol='='))
            ipcr_report.write(self._PCR_Simulation.all_products_histogram())
            ipcr_report.write('\n\n\n')
            #per hit histogram and phoresis
            ipcr_report.write(hr(' histograms and electrophorograms of PCR products of each hit ', symbol='='))
            ipcr_report.write(self._PCR_Simulation.all_graphs_grouped_by_hit())            
        ipcr_report.close()
        print '\niPCR report was written to:\n   ',self._PCR_report_filename
        self._add_report('iPCR report', self._PCR_report_filename)
    #end def
#end class



#tests
if __name__ == '__main__':
    from Primer import Primer, load_sequence
    import TD_Functions
    import cProfile
    import os, sys
    os.chdir('../')
    fwd_primer = Primer(load_sequence('ATATTCTACRACGGCTATCC', 'fwd_test', 'fwd_test'), 0.43e-6)
    rev_primer = Primer(load_sequence('GAASGCRAAKATYGGGAAC', 'rev_test', 'rev_test'), 0.43e-6)
    ipcr = iPCR('test-job', 
               fwd_primer, 
               rev_primer, 
               50, 
               1500, 
               40000,
               False, 
               33,
               include_side_annealings=True)
    TD_Functions.PCR_T = 53
    TD_Functions.C_Mg  = 3e-3
    TD_Functions.C_dNTP = 300e-6
    TD_Functions.C_DNA = 1e-10
    
#    out_file = open('iPCR.out', 'w')
#    stdout = sys.stdout 
#    sys.stdout = out_file
    
#    l1 = [1,]*1000000000
#    l2 = [2,]*1000000000
#    l3 = [1,]*1000000000
#    def iadd(x, y): 
#        x += y
#        return x
#    cProfile.run('''
#iadd(l3,l2)
#l1.extend(l2)
#''',
#'iadd-extend.profile')
    
    cProfile.run('''
ipcr.find_possible_products(('ThGa.fa',), 5)
ipcr.simulate_PCR()''',
    'iPCR.profile')
    
#    sys.stdout = stdout
#    out_file.close()
    
    
    ipcr.write_products_report()
    ipcr.write_PCR_report()