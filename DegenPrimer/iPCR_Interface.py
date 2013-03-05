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
Created on Dec 15, 2012

@author: Allis Tauri <allista@gmail.com>
'''

from ReporterInterface import ReporterInterface
from PCR_Simulation import PCR_Simulation


class iPCR_Interface(ReporterInterface):
    '''Interface class to create facilities which use PCR_Simulation'''

    def __init__(self, 
                 fwd_primer, 
                 rev_primer, 
                 min_amplicon, 
                 max_amplicon, 
                 polymerase, 
                 with_exonuclease,
                 num_cycles,
                 side_reactions=None,
                 side_concentrations=None):
        ReporterInterface.__init__(self)
        #parameters
        self._fwd_primer          = fwd_primer
        self._rev_primer          = rev_primer
        self._primers             = []
        if fwd_primer: self._primers.append(fwd_primer)
        if rev_primer: self._primers.append(rev_primer)
         
        self._min_amplicon        = min_amplicon
        self._max_amplicon        = max_amplicon
        
        self._polymerase          = polymerase
        self._with_exonuclease    = with_exonuclease
        
        self._num_cycles          = num_cycles
        
        self._side_reactions      = side_reactions
        self._side_concentrations = side_concentrations
        
        self._PCR_report_filename = None
        self._possible_products   = None
        self._have_results        = False
    #end def
    
    
    #factory for PCR_Simulation objects
    def _PCR_Simulation_factory(self):
        _PCR_simulation = PCR_Simulation(self._primers,
                                         self._min_amplicon,
                                         self._max_amplicon,
                                         self._polymerase,
                                         self._with_exonuclease,
                                         self._num_cycles)
        if self._side_reactions:
            _PCR_simulation.add_side_reactions(self._side_reactions)
        if self._side_concentrations:
            _PCR_simulation.add_side_concentrations(self._side_concentrations)
        return _PCR_simulation
    #end def
    
    
    def _add_products(self, simulation, hit, fwd_annealings, rev_annealings):
        '''
        Add possible products to the simulation:
        simulation - a PCR_Simulation object to add products to
        hit - name of a target sequence
        fwd_annealings - a list of primer annealing sites on a direct strand 
        of the target sequence: [(position, [(duplex, primer_id), ...]), ...]
        rev_annealings - a list of primer annealing sites on a reverse strand 
        of the target sequence (the structure is the same).
        Return True if some products were added. False otherwise.
        '''
        products_added = False
        #sort lists
        fwd_annealings.sort(key=lambda x: x[0])
        rev_annealings.sort(key=lambda x: x[0])
        #check only those annealing that are in range [min_amplicon, max_amplicon]
        rev_start = 0
        for fwd_p in fwd_annealings:
            start = fwd_p[0]+1
            for rev_pi, rev_p in enumerate(rev_annealings[rev_start:]):
                end = rev_p[0]-1
                if start >= end:
                    rev_start = rev_pi
                    continue
                if not (self._min_amplicon <= end-start+1 <= self._max_amplicon): 
                    break
                if simulation.add_product(hit, start, end, fwd_p[1], rev_p[1]):
                    if not products_added: products_added = True
        return products_added
    #end def
    
    
    #property functions
    def have_results(self): return self._have_results
    
    
    #pure virtual
    def simulate_PCR(self): pass
#end class