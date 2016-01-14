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

from abc import ABCMeta, abstractmethod


class ReporterInterface(object):
    '''Provides common interface to register reports if DegenPrimerConfig'''
    __metaclass__ = ABCMeta
    
    def __init__(self):
        self._reports = []
        self._have_results = False
    #end def
    
    def have_results(self): return self._have_results
    
    @abstractmethod
    def write_reports(self): pass
         
    #property function to use with a Proxy
    def reports(self): return self._reports
    
    
    def _add_report(self, name, filename):
        self._reports.append((name, filename))
    
    
    def _open_report(self, name, filename):
        try:
            report_handle = open(filename, 'w')
        except IOError, e:
            print '\nFailed to open %s report file for writing:\n   %s' \
                    % (name, filename)
            print e.message
            return
        return report_handle
    #end def
#end class