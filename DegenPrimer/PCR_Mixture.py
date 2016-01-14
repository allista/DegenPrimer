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

from AbortableBase import AbortableBase
from Product import Product, Region
from tmpStorage import tmpDict, roDict, tupleView


class PCR_Mixture(object):
    __slots__ = ['reaction_id', 'annealings', 'templates', 'products']
    def __init__(self, reaction_id):
        self.reaction_id = reaction_id
        self.annealings  = []       #all unique sites of primer annealing
        self.templates   = [[],[]]  #all template regions participating in the reaction
        self.products    = dict()   #all possible products of the reaction
    #end def
    

    def add_annealing(self, duplexes, template):
        self.annealings.append((duplexes, template))
        strand = template.forward
        self.templates[strand].append(template)
        self.templates[strand].sort(key=lambda(x): x.start)
        self._compact_templates(strand)
    #end def
    
    
    def add_annealings(self, annealings):
        self.annealings.extend(annealings)
        fwd_len = len(self.templates[1])
        rev_len = len(self.templates[0])
        for _dups, template in annealings:
            self.templates[template.forward].append(template)
        if len(self.templates[1]) > fwd_len:
            self.templates[1].sort(key=lambda(x): x.start)
            self._compact_templates(1)
        if len(self.templates[0]) > rev_len:
            self.templates[0].sort(key=lambda(x): x.start)
            self._compact_templates(0)
    #end def
    
    
    def _compact_templates(self, strand):
        compacted = [self.templates[strand][0]]
        for T in self.templates[strand]:
            if T.overlaps(compacted[-1]):
                compacted[-1] += T
            else: compacted.append(T)
        self.templates[strand] = compacted
    #end def
    
    
    def add_product(self, 
                      start,        #start position of the product on the target sequence
                      end,          #end position of the product on the target sequence
                      fwd_duplexes, #duplexes of forward primers with corresponding template sequence
                      rev_duplexes, #duplexes of reverse primers with corresponding template sequence
                     ):
        '''
        Add a new product.
        start, end -- positions in the sequence (starting form 1) corresponding 
        to the beginning and the end of PCR prduct
        *_duplexes -- lists of tuples of the form (duplex_object, id), where 
        duplex_object represents annealed primer and id is it's name 
        '''
        #initialize new product and reset nonzero flag
        new_product = Product(self.reaction_id, start, end, fwd_duplexes, rev_duplexes)
        #if no such product from this hit, add it to the list
        new_product_hash = hash(new_product)
        if new_product_hash not in self.products:
            self.products[new_product_hash]  = new_product
        else: #append primers
            self.products[new_product_hash] += new_product
    #end def
    
    
    def save(self):
        d = tmpDict(persistent=True)
        d['reaction_id']   = self.reaction_id
        d['fwd_templates'] = self.templates[1]
        d['rev_templates'] = self.templates[0]
        for p_id, product in self.products.iteritems():
            d['p%d'%p_id] = product
        for i, annealing in enumerate(self.annealings):
            d['a%d'%i] = annealing
        d.close()
        return d.filename
    #end def
#end class


class ShelvedMixture(object):
    '''A readonly interface to PCR_Mixture saved in a DB'''
    __slots__ = ['_d', 'reaction_id', 'annealings', 'templates', 'products']
    
    def __init__(self, fpath):
        self._d = roDict(fpath)
        self.reaction_id  = self._d['reaction_id']
        self.templates    = [[],[]]
        self.templates[1] = self._d['fwd_templates']
        self.templates[0] = self._d['rev_templates']
        self.products     = dict()
        annealing_keys    = []
        for k in self._d.iterkeys():
            if k.startswith('p'):
                self.products[int(k.strip('p'))] = self._d[k]
            elif k.startswith('a'):
                annealing_keys.append(k)
        self.annealings   = tupleView(annealing_keys, self._d)
    #end def
#end class


class MixtureFactory(AbortableBase):
    def __init__(self, 
                  abort_event,
                  primers,                #all primers (generally degenerate) that may be present in the system
                  min_amplicon,           #minimum amplicon length 
                  max_amplicon,           #maximum amplicon length
                  with_exonuclease=False, #if polymerase does have have 3' exonuclease activity, products of primers with 3'-mismatches will also be icluded
                  include_side_annealings=False, #if True, include annealings which do not give any products into equilibrium system as side reactions
                  ):
        AbortableBase.__init__(self, abort_event)
        self._primers                 = primers
        self._min_amplicon            = min_amplicon
        self._max_amplicon            = max_amplicon
        self._with_exonuclease        = with_exonuclease
        self._include_side_annealings = include_side_annealings
    #end def
    
    @staticmethod
    def _add_side_annealings(conditions, sorted_annealings):
        conditions.add_annealings(tuple((d, t) for _p, d, t in sorted_annealings))
    #end def
            
    
    def _sort_annealings(self, hit, strand, annealings, good, bad):
        for _pos, _duplexes in annealings:
            if strand:
                _template = Region(hit, _pos-_duplexes[0][0].fwd_len+1, _pos, strand)
            else:
                _template = Region(hit, _pos, _pos+_duplexes[0][0].fwd_len-1, strand)
            good_annealings = (_pos, [], _template)
            bad_annealings  = (_pos, [], _template)
            for _duplex, _id in _duplexes:
                #check if there are such primers in the system at all
                duplex_with_primer = False
                for primer in self._primers:
                    if _duplex.fwd_seq in primer:
                        duplex_with_primer = True
                        break
                if not duplex_with_primer: continue
                #check 3' mismatches
                if self._with_exonuclease: #if poly has 3'-5'-exonuclease, add to the good annealings
                    good_annealings[1].append((_duplex, _id))
                elif _duplex.have_3_matches: #if not, check 3' annealing
                    good_annealings[1].append((_duplex, _id))
                else:
                    bad_annealings[1].append((_duplex, _id))
            #add annealings
            if good_annealings[1]: good[strand].append(good_annealings)
            if bad_annealings[1]:  bad.append(bad_annealings)
    #end def
    

    def create_PCR_mixture(self, counter, hit, fwd_annealings, rev_annealings):
        '''Within annealing sites of primers find PCR products; add other 
        annealings as side reactions.
        hit - name of a target sequence
        fwd_annealings - a list of primer annealing sites on a direct strand 
        of the target sequence: [(position, [(duplex, primer_id), ...]), ...]
        rev_annealings - a list of primer annealing sites on a reverse strand 
        of the target sequence (the structure is the same).
        Return PCR_Mixture object if some products produced by the annealings 
        were found. None otherwise.'''
        if not fwd_annealings or not rev_annealings: 
            counter.done(); return None
#        print 'PCR Simulation: searching for possible PCR products in %s...' % hit
        if self._include_side_annealings: counter.set_work(3)
        else: counter.set_work(2)
        mixture = PCR_Mixture(hit)        
        hit = str(hit) #if hit is unicode
        #sort annealings into good (ones that are suitable for product generation) 
        #and bad (ones that are not)
        good_annealings = ([], []) #0 is reverse strand, 1 is forward
        bad_annealings  = [] 
        self._sort_annealings(hit, 1, fwd_annealings, 
                              good_annealings, bad_annealings); del fwd_annealings
        self._sort_annealings(hit, 0, rev_annealings, 
                              good_annealings, bad_annealings); del rev_annealings
        #sort good annealings by position
        good_annealings[1].sort(key=lambda x: x[0])
        good_annealings[0].sort(key=lambda x: x[0])
        counter.count()
        #find possible products in range [min_amplicon, max_amplicon]
        products_added  = False
        added_positons  = (set(), set())
        rev_start       = 0
        for fwd_pi, (fwd_pos, fwd_dups, fwd_templ) in enumerate(good_annealings[1]):
            start = fwd_pos+1
            for rev_pi, (rev_pos, rev_dups, rev_templ) in enumerate(good_annealings[0][rev_start:]):
                end = rev_pos-1
                if start >= end:
                    rev_start = rev_pi
                    continue
                if not (self._min_amplicon <= end-start+1 <= self._max_amplicon): 
                    break
                mixture.add_product(start, end, fwd_dups, rev_dups)
                if fwd_pi not in added_positons[1]:
                    mixture.add_annealing(fwd_dups, fwd_templ)
                    added_positons[1].add(fwd_pi)
                if rev_pi not in added_positons[0]:
                    mixture.add_annealing(rev_dups, rev_templ)
                    added_positons[0].add(rev_pi)
                if not products_added: products_added = True
        counter.count()
        if products_added: #if some products were found, set nonzero flag and save PCR mixture to the DB
            if self._include_side_annealings: #add side annealings if requested
                self._add_side_annealings(mixture, [ann for i, ann in enumerate(good_annealings[1])
                                                       if i not in added_positons[1]])
                self._add_side_annealings(mixture, [ann for i, ann in enumerate(good_annealings[0])
                                                       if i not in added_positons[0]])
                self._add_side_annealings(mixture, bad_annealings)
                counter.count()
#            print 'PCR Simulation: found some possible products in %s.' % hit
            return mixture
#        print 'PCR Simulation: no products were found in %s.' % hit
        return None
    #end def
#end class
