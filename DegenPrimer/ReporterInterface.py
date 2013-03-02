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

class ReporterInterface(object):
    '''Provides common interface to register reports if DegenPrimerConfig'''
    
    def __init__(self):
        self._reports = []
    
         
    #property function to use with a Proxy
    def reports(self): return self._reports
    
    
    def _add_report(self, name, filename):
        self._reports.append({'report_name': name, 'report_file': filename})
    
    
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