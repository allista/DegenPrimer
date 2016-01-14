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
from SearchEngine import SearchEngine
from BioUtils.Tools.Multiprocessing import MultiprocessingBase
from WorkCounter import WorkCounter

class PCR_ProductsFinder(MixtureFactory, MultiprocessingBase):
    
    def __init__(self, abort_event, *args, **kwargs):
        MixtureFactory.__init__(self, abort_event, *args, **kwargs)
        MultiprocessingBase.__init__(self, abort_event)
        self._searcher  = SearchEngine(self._abort_event)
        self._p_weights = [p.num_components for p in self._primers]
        self._pw_sum    = sum(self._p_weights)
        self._num_p     = len(self._primers)
    #end def 
    
    @property
    def searcher(self): return self._searcher
    
    @staticmethod
    def _combine_annealings(*annealings_list):
        fwd_annealings = []
        rev_annealings = []
        for fwd, rev in annealings_list:
            fwd_annealings.extend(fwd)
            rev_annealings.extend(rev)
        return fwd_annealings,rev_annealings
    #end def
    
    
    def matches_to_mixture(self, counter, t_name, matches_list):
        all_annealings = []
        m_weights = [len(m) for m in matches_list]
        counter.set_subwork(len(matches_list)+1,
                            m_weights+[0.01*sum(m_weights)])
        for i, matches in enumerate(matches_list):
            if self.aborted(): return None
            annealings = self._searcher.compile_duplexes_mp(counter[i], *matches)
            if annealings is None: return None
            all_annealings.append(annealings)
        mixture = self.create_PCR_mixture(counter[-1], t_name, 
                                          *self._combine_annealings(*all_annealings))
        counter.done()
        return mixture
    #end def
    
    
    def find_matches(self, counter, t_name, template, mismatches):
        all_matches = []
        counter.set_subwork(self._num_p, self._p_weights)
        for i, primer in enumerate(self._primers):
            if self.aborted(): return None
            matches = self._searcher.find_matches(counter[i], template, primer, mismatches)
            if matches is None: return None
            all_matches.append(matches)
        return all_matches
    #end def
    
    
    def find(self, counter, t_name, template, mismatches):
        all_annealings = []
        counter.set_subwork(self._num_p+1,
                            self._p_weights+[0.01*self._pw_sum])
        for i, primer in enumerate(self._primers):
            if self.aborted(): return None
            annealings = self._searcher.find(counter[i], template, primer, mismatches)
            if annealings is None: return None
            all_annealings.append(annealings)
        mixture = self.create_PCR_mixture(counter[-1], t_name, 
                                          *self._combine_annealings(*all_annealings))
        counter.done()
        return mixture
    #end def
    
    
    def batch_find(self, counter, t_names, templates, mismatches):
        results        = dict()
        search_results = dict()
        counter.set_subwork(2)
        counter[0].set_subwork(self._num_p)
        for i, primer in enumerate(self._primers):
            if self.aborted(): return None
            result = self._searcher.batch_find(counter[0][i], templates, primer, mismatches)
            if result is None: return None
            for t_id, annealings in result.iteritems():
                if t_id in search_results:
                    search_results[t_id] = self._combine_annealings(search_results[t_id], annealings)
                else: search_results[t_id] = annealings
        print '\nPSR Simulation: creating PCR mixtures...\n'
        @MultiprocessingBase.data_mapper
        def _worker(i, factory, search_results, t_names):
            t_id, (fwd_annealings, rev_annealings) = search_results[i]
            t_name = t_names[t_id]
            mixture = factory.create_PCR_mixture(WorkCounter(), t_name, fwd_annealings, rev_annealings)
            if not mixture: return t_name, None
            return t_name, mixture.save()
        @MultiprocessingBase.results_assembler
        def _assembler(index, result, results):
            if result[1]: results[result[0]] = result[1]
        work = self.Work(counter=counter[1])
        work.prepare_jobs(_worker, range(len(search_results)), None,
                          self, search_results.items(), t_names)
        work.set_assembler(_assembler, results)
        self.start_work(work)
        if not self.wait(work): return None
        if not results: return None
        return results
    #end def
#end class


#ProductsFnder wrapped in a Manager
def _find_matches(counter, t_name, template, mismatches, products_finder):
    return products_finder.find_matches(counter, t_name, template, mismatches)

def _matches_to_mixture(counter, t_name, matches_list, products_finder):
    mixture = products_finder.matches_to_mixture(counter, t_name, matches_list)
    if mixture is None: return None
    return mixture.save()
#end def

def _find(counter, t_name, template, mismatches, products_finder):
    mixture = products_finder.find(counter, t_name, template, mismatches)
    if mixture is None: return None
    return mixture.save()
#end def

def _batch_find(counter, t_names, templates, mismatches, products_finder):
    return products_finder.batch_find(counter, t_names, templates, mismatches)
#end def
        
#mp computations
from BioUtils.Tools.UMP import FuncManager
PPFManager = FuncManager('PPFManager', 
                         (_matches_to_mixture, 
                          _find_matches, 
                          _find, 
                          _batch_find))