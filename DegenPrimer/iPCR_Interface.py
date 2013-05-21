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
Created on Dec 15, 2012

@author: Allis Tauri <allista@gmail.com>
'''

from ReporterInterface import ReporterInterface
from PCR_Simulation import PCR_Simulation
from StringTools import hr, time_hr


class iPCR_Interface(ReporterInterface):
    '''Interface class to create facilities which use PCR_Simulation'''

    _PCR_report_suffix = 'iPCR'

    def __init__(self,
                 job_id,  
                 primers, 
                 min_amplicon, 
                 max_amplicon, 
                 polymerase, 
                 with_exonuclease,
                 num_cycles,
                 side_reactions=None,
                 side_concentrations=None,
                 include_side_annealings=False):
        ReporterInterface.__init__(self)
        #parameters
        try:
            if len(primers) == 0: 
                raise ValueError('iPCR_Interface: no primers given.')
        except TypeError:
            raise TypeError(('iPCR_Interface: primers should be an iterable. '
                            'Given %s instead') % str(primers))
        
        self._job_id              = job_id
        self._primers             = primers
         
        self._min_amplicon        = min_amplicon
        self._max_amplicon        = max_amplicon
        
        self._polymerase          = polymerase
        self._with_exonuclease    = with_exonuclease
        
        self._num_cycles          = num_cycles
        
        self._side_reactions      = side_reactions
        self._side_concentrations = side_concentrations
        self._include_side_annealings = include_side_annealings
        
        #report files
        self._PCR_report_filename = '%s-%s-report.txt' % (self._job_id, 
                                                          self._PCR_report_suffix)
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
                                         self._num_cycles,
                                         self._include_side_annealings)
        if self._side_reactions:
            _PCR_simulation.add_side_reactions(self._side_reactions)
        if self._side_concentrations:
            _PCR_simulation.add_side_concentrations(self._side_concentrations)
        return _PCR_simulation
    #end def
    
    
    #property functions
    def have_results(self): return self._have_results
    
    
    #virtual
    def simulate_PCR(self): pass
    
    def _format_header(self): return ''
    
    def _format_report_body(self): return ''
    
    
    def write_report(self):
        if not self._have_results: return
        #open report file
        report = self._open_report(self._PCR_report_suffix, self._PCR_report_filename)
        if not report: return
        #header
        report.write(time_hr())
        report.write(self._format_header())
        #if no PCR products have been found
        if not self._have_results:
            report.write(hr(' No PCR products have been found ', symbol='!'))
            report.close()
            return
        #report body
        report.write(self._format_report_body())
        report.close()
        print '\n%s report was written to:\n   %s' % (self._PCR_report_suffix,
                                                      self._PCR_report_filename)
        self._add_report(self._PCR_report_suffix.replace('_', ' '), self._PCR_report_filename)
    #end def
#end class