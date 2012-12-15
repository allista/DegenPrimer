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

from PCR_Simulation import PCR_Simulation


class iPCR_Interface(object):
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
        #parameters
        self._fwd_primer          = fwd_primer
        self._rev_primer          = rev_primer
        self._min_amplicon        = min_amplicon
        self._max_amplicon        = max_amplicon
        self._polymerase          = polymerase
        self._with_exonuclease    = with_exonuclease
        self._num_cycles          = num_cycles
        self._side_reactions      = side_reactions
        self._side_concentrations = side_concentrations
    #end def
    
    
    def _PCR_Simulation_factory(self):
        _PCR_simulation = PCR_Simulation([self._fwd_primer, 
                                          self._rev_primer],
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
    
    
    def simulate_PCR(self): pass
#end class