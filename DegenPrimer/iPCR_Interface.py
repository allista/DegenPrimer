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

from abc import ABCMeta, abstractmethod

from BioUtils.Tools.Text import hr, time_hr
from BioUtils.Tools.Debug import raise_tb_on_error
from BioUtils.Tools.Output import simple_timeit

from .PCR_ProductsFinder import PCR_ProductsFinder
from .PCR_Simulation import PCR_Simulation_Interface, PCR_Simulation
from .ReporterInterface import ReporterInterface
from . import TD_Functions as tdf

class iPCR_Interface(PCR_Simulation_Interface, ReporterInterface):
    '''Interface class to create facilities which use PCR_Simulation'''
    __metaclass__ = ABCMeta

    _PCR_report_suffix = 'iPCR'

    def __init__(self, abort_event, job_id, *args, **kwargs):
        ReporterInterface.__init__(self)
        PCR_Simulation_Interface.__init__(self, abort_event, *args, **kwargs)
        self._job_id              = job_id
        self._PCR_report_filename = '%s-%s-report.txt' % (self._job_id, 
                                                          self._PCR_report_suffix)
        self._PCR_products_filename = '%s-%s-products.txt' % (self._job_id, 
                                                              self._PCR_report_suffix)
    #end def
    
    
    #factory for PCR_ProductsFinder objects
    def _new_PCR_ProductsFinder(self):
        return PCR_ProductsFinder(self._abort_event,
                                  self._primers,
                                  self._min_amplicon, 
                                  self._max_amplicon,
                                  self._with_exonuclease,
                                  self._include_side_annealings)
    #end def
    
    #factory for PCR_Simulation objects
    def _new_PCR_Simulation(self):
        return PCR_Simulation(self._abort_event,
                              self._primers,
                              self._min_amplicon,
                              self._max_amplicon,
                              self._polymerase,
                              self._with_exonuclease,
                              self._num_cycles,
                              self._side_reactions,
                              self._side_concentrations,
                              self._include_side_annealings)
    #end def
    
    
    @abstractmethod
    def simulate_PCR(self): pass
    
    @abstractmethod
    def _format_header(self): return ''
    
    @abstractmethod
    def _format_report_body(self): return ''
    
    @abstractmethod
    def write_products_report(self): pass
    
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
        report_name = self._PCR_report_suffix.replace('_', ' ').replace('-', ' ')
        self._add_report(report_name, self._PCR_report_filename)
    #end def
    
    @raise_tb_on_error
    def write_reports(self):
        with tdf.AcquireParameters():
            self.write_products_report()
            self.write_report()
    #end def
#end class