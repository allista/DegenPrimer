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

from BioUtils.Tools.Multiprocessing import MultiprocessingBase
from BioUtils.Tools.tmpStorage import from_shelf, cleanup_files
from BioUtils.SeqUtils import pretty_rec_name 

from .PCR_Mixture import MixtureFactory
from .SearchEngine import SearchEngine
from .WorkCounter import WorkCounter

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
    def _extract_annealings(*annealings_list):
        fwd_annealings = []
        rev_annealings = []
        for filename in annealings_list:
            ann = from_shelf(filename)
            if not ann: continue
            if ann[0]: fwd_annealings.extend(ann[0]) 
            if ann[1]: rev_annealings.extend(ann[1])
        return fwd_annealings,rev_annealings
    #end def
    
    @staticmethod
    def _combine_annealings(*annealings_list):
        fwd_annealings = []
        rev_annealings = []
        for ann in annealings_list:
            if ann[0]: fwd_annealings.extend(ann[0])
            if ann[1]: rev_annealings.extend(ann[1])
        return fwd_annealings,rev_annealings
    #end def
    
    def matches_to_mixture(self, counter, tname, matches_list):
        all_annealings = []
        m_weights = [len(m) for m in matches_list]
        counter.set_subwork(len(matches_list)+1,
                            m_weights+[0.01*sum(m_weights)])
        for i, matches in enumerate(matches_list):
            if self.aborted(): return None
            annealings = self._searcher.compile_duplexes_mp(counter[i], *matches)
            if annealings is None: return None
            all_annealings.append(annealings)
        mixture = self.create_PCR_mixture(counter[-1], tname, 
                                          *self._combine_annealings(*all_annealings))
        counter.done()
        return mixture.save()
    #end def
    
    def find_matches(self, counter, template, mismatches):
        all_matches = []
        counter.set_subwork(self._num_p, self._p_weights)
        for i, primer in enumerate(self._primers):
            if self.aborted(): return None
            matches = self._searcher.find_matches(counter[i], template, primer, mismatches)
            if matches is None: return None
            all_matches.append(matches)
        return all_matches
    #end def
    
    def find(self, counter, template, mismatches):
        all_annealings = []
        counter.set_subwork(self._num_p+1,
                            self._p_weights+[0.01*self._pw_sum])
        for i, primer in enumerate(self._primers):
            if self.aborted(): return None
            annealings = self._searcher.find(counter[i], template, primer, mismatches)
            if annealings is None: continue
            all_annealings.append(annealings)
        if not all_annealings: return None
        tname = pretty_rec_name(template)
        mixture = self.create_PCR_mixture(counter[-1], tname, 
                                          *self._extract_annealings(*all_annealings))
        cleanup_files(all_annealings)
        counter.done()
        if not mixture: return None
        return {tname: mixture.save()}
    #end def
    
    def batch_find(self, counter, templates, mismatches, **kwargs):
        @MultiprocessingBase.data_mapper
        def worker(template):
            search_results = []
            for primer in self._primers:
                matches = self._searcher.find_matches(WorkCounter(), template, primer, mismatches)
                if not matches: continue
                duplexes = self._searcher.compile_duplexes(WorkCounter(), *matches)
                if duplexes: search_results.append(duplexes)
            if not search_results: return None
            all_annealings = self._combine_annealings(*search_results)
            tname = pretty_rec_name(template)
            mixture = self.create_PCR_mixture(WorkCounter(), tname, all_annealings[0], all_annealings[1])
            if not mixture: return None
            return tname, mixture.save()
        @MultiprocessingBase.results_assembler
        def assembler(index, result, results):
            if result: results[result[0]] = result[1]
        results = dict()
        counter.set_subwork(1)
        work = self.Work(counter=counter[0])
        work.start_work(worker, templates, None, **kwargs)
        work.assemble(assembler, results)
        if not self.wait(work): return None
        if not results: return None
        return results
    #end def
#end class


#ProductsFnder wrapped in a Manager
def _find_matches(counter, template, mismatches, products_finder):
    return products_finder.find_matches(counter, template, mismatches)

def _matches_to_mixture(counter, t_name, matches_list, products_finder):
    mixture = products_finder.matches_to_mixture(counter, t_name, matches_list)
    if mixture is None: return None
    return mixture
#end def

def _find(counter, template, mismatches, products_finder):
    mixture = products_finder.find(counter, template, mismatches)
    if mixture is None: return None
    return mixture
#end def

def _batch_find(counter, templates, mismatches, products_finder, **kwargs):
    return products_finder.batch_find(counter, templates, mismatches, **kwargs)
#end def
        
#mp computations
from BioUtils.Tools.UMP import FuncManager
PPFManager = FuncManager('PPFManager', 
                         (_matches_to_mixture, 
                          _find_matches, 
                          _find, 
                          _batch_find))