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

import os

from BioUtils.Tools import mktmp_name, tmpStorage

from .StringTools import print_exception, time_hr, hr
from .iPCR_Interface import iPCR_Interface
from .SeqDB import SeqDB


class iPCR_Base(iPCR_Interface):
    '''
    Using PCR_Simulation and SeqDB classes runs PCR simulation with given 
    primers and report results in human readable form in a text file.
    '''
    
    def __init__(self, abort_event, max_mismatches, *args, **kwargs):
        iPCR_Interface.__init__(self, abort_event, *args, **kwargs)
        self._max_mismatches = max_mismatches
        self._seq_db         = SeqDB(self._abort_event)
        self._seq_names      = None
        self._PCR_Simulation = None
    #end def
    
    
    def __del__(self):
        self._seq_db.close()
        try: self._searcher.shutdown()
        except: pass
    #end def
    
    
    def _try_connect_db(self, seq_files):
        try: #to connect to a database
            if len(seq_files) == 1 and seq_files[0].endswith('.db'):
                if os.path.isfile(seq_files[0]):
                    if not self._seq_db.connect(seq_files[0]):
                        print '\niPCR: unable to connect to %s' % seq_files[0]
                        return False
                else: 
                    print '\niPCR: no such file: %s' % seq_files[0]
                    return False
            #or make one in temp file
            else:
                dbfile = mktmp_name('.db')
                tmpStorage.register_tmp_file(dbfile)
                if not self._seq_db.create_db_from_files(dbfile, seq_files): 
                    print '\niPCR: unable to create sequence database from given files.'
                    tmpStorage.cleanup_file(dbfile)
                    return False
        except Exception, e:
            print_exception(e)
            return False
        return True
    #end def
    
    
    def _get_names(self, seq_ids):
        seq_names = self._seq_db.get_names(seq_ids)
        if not seq_names:
            if seq_ids:
                print ('\nPCR Simulation: there\'re no sequences with ' 
                       'ids in %s in the database.') % str(seq_ids)
            else: print '\nPCR Simulation: the database is empty.'
            return None
        return seq_names
    #end def
    
    
    def _format_header(self):
        header = iPCR_Interface._format_header(self)
        if self._max_mismatches != None:
            header += 'Number of mismatches allowed: %d\n\n' % self._max_mismatches
        return header
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
#end class