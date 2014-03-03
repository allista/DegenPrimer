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
Created on Mar 1, 2014

@author: Allis Tauri <allista@gmail.com>
'''

from PCR_Mixture import MixtureFactory
import SearchEngine


class PCR_ProductsFinder(MixtureFactory):
    
    @staticmethod
    def _combine_annealings(*annealings_list):
        fwd_annealings = []
        rev_annealings = []
        for fwd, rev in annealings_list:
            fwd_annealings.extend(fwd)
            rev_annealings.extend(rev)
        return fwd_annealings,rev_annealings
    #end def
    
    
    def find(self, t_name, template, mismatches):
        search_results  = []
        for primer in self._primers:
            if self._abort_event.is_set(): return None
            result = SearchEngine.find(template, primer, mismatches, self._abort_event)
            if result is None: return None
            search_results.append(result)
        return self.create_PCR_mixture(t_name, *self._combine_annealings(*search_results))
    #end def
    
    
    def batch_find(self, t_names, templates, mismatches):
        results        = dict()
        search_results = dict()
        for primer in self._primers:
            if self._abort_event.is_set(): return None
            result = SearchEngine.batch_find(templates, primer, mismatches, self._abort_event)
            if result is None: return None
            for t_id, annealings in result.iteritems():
                if t_id in search_results:
                    search_results[t_id] = self._combine_annealings(search_results[t_id], annealings)
                else: search_results[t_id] = annealings
        for t_id, (fwd_annealings, rev_annealings) in search_results.items():
            mixture = self.create_PCR_mixture(t_names[t_id], fwd_annealings, rev_annealings)
            if mixture is not None: results[t_names[t_id]] = mixture 
        if not results: return None
        return results
    #end def
#end class


#ProductsFnder wrapped in a Manager
def _find(t_name, template, mismatches, products_finder):
    mixture = products_finder.find(t_name, template, mismatches)
    if mixture is None: return None
    return mixture.save()
#end def

def _batch_find(t_names, templates, mismatches, products_finder):
    mixtures = products_finder.batch_find(t_names, templates, mismatches)
    if mixtures is None: return None
    return dict((k, m.save()) for k, m in mixtures.iteritems())
#end def
        
#mp computations
from MultiprocessingBase import FuncManager
PPFManager = FuncManager((_find, _batch_find))
