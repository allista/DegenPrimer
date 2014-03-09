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
from StringTools import print_exception
from iPCR_Interface import iPCR_Interface
from SeqDB import SeqDB
from MultiprocessingBase import FuncManager, cpu_count, even_chunks
from SearchEngine import mp_better
from PCR_ProductsFinder import PPFManager


#manager to isolate forking point in a low-memory process
def _find(counter, t_name, template, mismatches, products_finder):
    finder = PPFManager()
    finder.start()
    result_path = finder.find(counter, t_name, template, mismatches, products_finder)
    result_path = result_path._getvalue()
    finder.shutdown()
    return result_path
#end def

def _batch_find(counter, t_names, templates, mismatches, products_finder):
    finder = PPFManager()
    finder.start()
    result_paths = finder.batch_find(counter, t_names, templates, mismatches, products_finder)
    result_paths = result_paths._getvalue()
    finder.shutdown()
    return result_paths
#end def

SearchManager = FuncManager((_find, _batch_find))


class iPCR_Base(iPCR_Interface):
    '''
    Using PCR_Simulation and SeqDB classes runs PCR simulation with given 
    primers and report results in human readable form in a text file.
    '''
    
    def __init__(self, abort_event, max_mismatches, *args, **kwargs):
        iPCR_Interface.__init__(self, abort_event, *args, **kwargs)
        self._max_mismatches = max_mismatches
        self._seq_db         = SeqDB(self._abort_event)
        self._searcher       = SearchManager()
        self._seq_names      = None
        self._PCR_Simulation = None
        self._searcher.start()
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
                print ('\nPCR Simulation: there\'re no sequences with ' 
                       'ids in %s in the database.') % str(seq_ids)
            else: print '\nPCR Simulation: the database is empty.'
            return None
        return seq_names
    #end def
    
    
    def _find_products_in_templates(self, counter, t_ids, products_finder):
        templates   = self._seq_db.get_seqs(t_ids)
        if len(t_ids) < 2:
            result = self._searcher.find(counter, self._seq_names[t_ids[0]],
                                         templates[0][2], 
                                         self._max_mismatches, 
                                         products_finder)
        else:
            result = self._searcher.batch_find(counter, self._seq_names, 
                                               templates,  
                                               self._max_mismatches, 
                                               products_finder)
        return result._getvalue()
    #end def
    
    
    def _find_products_in_db(self, counter, PCR_Sim, P_Finder, seqs_info):
        print '\nPSR Simulation: searching for annealing sites in provided sequences...'
        #sort templates into short, suitable for plain search and long -- for mp
        short_templates = []
        long_templates  = []
        for t_id, t_len in seqs_info.iteritems():
            if mp_better(t_len):
                long_templates.append(t_id)
            else: short_templates.append(t_id)
        #setup work counters
        if short_templates and long_templates: 
            counter.set_subwork(2, (sum(seqs_info[t_id] for t_id in short_templates),
                                    sum(seqs_info[t_id] for t_id in long_templates)))
            short_counter  = counter[0]
            long_counter   = counter[1]
        else: long_counter = short_counter = counter
        if long_templates: long_counter.set_subwork(len(long_templates), [seqs_info[t_id] for t_id in long_templates])
        #If there's enough short templates, run a series of batch searches.
        #This is needed to lower memory usage pikes during search and to 
        #better utilize cpus.
        if short_templates:
            chunks = even_chunks(short_templates, 
                                 max(len(short_templates)/cpu_count, 1))
            start  = 0; short_counter.set_subwork(len(chunks))
            for i, chunk in enumerate(chunks):
                end = start+chunk
                results = self._find_products_in_templates(short_counter[i], 
                                                           short_templates[start:end], 
                                                           P_Finder)
                if results is None: 
                    if self._abort_event.is_set(): return False
                    else: continue
                for t_name, m_path in results.iteritems():
                    PCR_Sim.add_mixture(t_name, m_path)
        #if there're long templates, search sequentially
        for i, t_id in enumerate(long_templates):
            result = self._find_products_in_templates(long_counter[i], (t_id,), P_Finder)
            if result is None:
                if self._abort_event.is_set(): return False
                else: continue
            PCR_Sim.add_mixture(self._seq_names[t_id], result)
        return PCR_Sim.not_empty()
    #end def
    
    
    def _find_products(self, counter, PCR_Sim, P_Finder, seq_files, seq_ids=None):
        #try to connect to a database
        if not seq_files or not self._try_connect_db(seq_files): return False
        #get names and legths of sequences
        self._seq_names = self._get_names(seq_ids)
        if not self._seq_names:
            self._seq_db.close() 
            return False
        seq_lengths = self._seq_db.get_lengths(seq_ids)
        #find primer annealing sites
        result = self._find_products_in_db(counter, PCR_Sim, P_Finder, seq_lengths)
        self._seq_db.close()
        return result
    #end def

    
    def _format_header(self):
        header = iPCR_Interface._format_header(self)
        if self._max_mismatches != None:
            header += 'Number of mismatches allowed: %d\n\n' % self._max_mismatches
        return header
#end class