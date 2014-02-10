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
from StringTools import print_exception, hr, time_hr
from iPCR_Interface import iPCR_Interface
from SeqDB import SeqDB

class iPCR_Base(iPCR_Interface):
    '''
    Using PCR_Simulation and SeqDB classes runs PCR simulation with given 
    primers and report results in human readable form in a text file.
    '''
    
    def __init__(self, max_mismatches, *args, **kwargs):
        iPCR_Interface.__init__(self, *args, **kwargs)
        self._max_mismatches    = max_mismatches
        #ipcr parameters
        self._seq_db            = SeqDB()
        self._seq_names         = None
        self._seq_annealings    = None
        #simulation
        self._PCR_Simulation    = None
    #end def
    
    
    def __del__(self):
        self._seq_db.abort_search()
        self._seq_db.close()
        if self._PCR_Simulation is not None:
            self._PCR_Simulation.abort()
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
            #or make one in the memory from provided sequences
            elif not self._seq_db.create_db_from_files(':memory:', seq_files): 
                print '\niPCR: unable to create sequence database from given files.'
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
                print ('\niPCR: there\'re no sequences with ' 
                       'ids in %s in the database.') % str(seq_ids)
            else: print '\niPCR: the database is empty.'
            return None
        return seq_names
    #end def
    
    
    def _find_annealings(self, seq_files, seq_ids=None):
        #try to connect to a database
        if not seq_files or not self._try_connect_db(seq_files): return False
        #setup parameters
        self._seq_names = self._get_names(seq_ids)
        if not self._seq_names: return False
        #find primer annealing sites
        primer_annealings = []
        for primer in self._primers:
            print '\n\nStarting search for annealing sites of %s\n' % str(primer)
            results = self._seq_db.find_in_db(primer, 
                                              self._max_mismatches, 
                                              self._seq_names.keys())
            print '\nResults obtained for %s\n' % str(primer)
            primer_annealings.append(results)
        self._seq_db.close()
        self._seq_annealings = dict()
        for seq_id in self._seq_names:
            fwd_annealings = []
            rev_annealings = []
            for annealings in primer_annealings:
                if annealings[seq_id]:
                    fwd_annealings.extend(annealings[seq_id][0])
                    rev_annealings.extend(annealings[seq_id][1])
            self._seq_annealings[seq_id] = (fwd_annealings, rev_annealings)
        return True
    #end def
    
    
    def _add_annealings_for_seqs(self, _PCR_Sim):
        products_added = False
        for seq_id, seq_name in self._seq_names.items():
            if _PCR_Sim.add_annealings(seq_name, 
                                       self._seq_annealings[seq_id][0], 
                                       self._seq_annealings[seq_id][1]):
                products_added |= True
        return products_added
    #end def
    
    
    def _format_header(self):
        header = iPCR_Interface._format_header(self)
        if self._max_mismatches != None:
            header += 'Number of mismatches allowed: %d\n\n' % self._max_mismatches
        return header
#end class