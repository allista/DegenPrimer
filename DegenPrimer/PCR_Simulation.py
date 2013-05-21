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
Created on Nov 16, 2012

@author: Allis Tauri <allista@gmail.com>
'''

from copy import deepcopy
from datetime import timedelta
from math import log, exp
from Equilibrium import Equilibrium
import StringTools
import TD_Functions
from TD_Functions import format_PCR_conditions
from StringTools import wrap_text, line_by_line, hr
from Product import Product, Region
from MultiprocessingBase import MultiprocessingBase


class PCR_Simulation(MultiprocessingBase):
    '''In silica PCR results filtration, computation and visualization.'''
    
    #histogram
    _all_products_title  = 'PCR-product'
    _hist_column_title   = '     final concentration     '
    #bar width should not be less than 2*len('100.00%')+1, or value will not fit
    _hist_width          = max(2*len('999.9 mM')+2, len(_hist_column_title)) 
    #electrophoresis
    _window_percent      = 0.05 #percentage of a length of the longest PCR product; it defines band width on electrophoresis
    #indirectly defines to what extent should lengths 
    #of two products differ in order of them to occupy 
    #separate bands on the electrophoresis.
    _precision           = 1e6
    #minimum equilibrium constant: reactions with EC less than this will not be taken into account
    _min_K               = 100
    #products with quantity less than maximum quantity multiplied by this factor will not be included in the report
    _min_quantity_factor = 0.001
    #maximum 3' matches to remove in alternative primer annealing
    _max_3_matches       = 3


    def __init__(self, 
                 primers,                #all primers (generally degenerate) that may be present in the system
                 min_amplicon,           #minimum amplicon length 
                 max_amplicon,           #maximum amplicon length
                 polymerase,             #polymerase concentration in Units per liter
                 with_exonuclease=False, #if polymerase does have have 3' exonuclease activity, products of primers with 3'-mismatches will also be icluded
                 num_cycles=20,          #number of PCR cycles
                 include_side_annealings=False, #if True, include annealings which do not give any products into equilibrium system as side reactions
                 ):
        MultiprocessingBase.__init__(self)
        self._primers               = primers
        self._min_amplicon          = min_amplicon
        self._max_amplicon          = max_amplicon
        self._polymerase            = polymerase
        self._with_exonuclease      = with_exonuclease
        self._include_side_annealings = include_side_annealings
        self._elongation_time       = max_amplicon/1000.0 #minutes
        self._num_cycles            = num_cycles

        self._annealings            = dict() #dictionary of all primer annealings
        self._templates             = dict() #dictionary of all templates
        self._used_primer_ids       = set()
        self._products              = dict() #dictionary of all possible products
        
        self._side_reactions        = dict()
        self._side_concentrations   = dict()
        self._primer_concentrations = dict()
        for primer in self._primers:
            self._primer_concentrations.update(dict.fromkeys(primer.str_sequences, 
                                                             primer.concentration))
        
        self._solutions             = dict() #solutions of equilibrium systems
        self._max_objective_value   = -1 #objective value of the worst solution
        self._reaction_ends         = dict()
        
        self._nonzero               = False  #true if there're some products and their quantities were calculated
    #end def
    
    
    def __nonzero__(self):
        return self._nonzero


    def hits(self):
        return list(self._products.keys())

    
    def product_concentration(self, hit_id, bounds):
        if self._nonzero and hit_id in self._products:
            lookup_region = Region(hit_id, *bounds, forward=True)
            total_concentration   = 0
            product_concentration = -1
            for _unused, product in self._products[hit_id].items():
                total_concentration += product.quantity
                if product == lookup_region:
                    product_concentration = product.quantity
            if total_concentration == 0 or product_concentration < 0: return 0, 0
            return (product_concentration, 
                    product_concentration/total_concentration)
        return 0, 0
    #end def
    
    
    def add_side_concentrations(self, concentrations):
        if concentrations: self._side_concentrations.update(concentrations)
        
        
    def add_side_reactions(self, reactions):
        if not reactions: return 
        self._side_reactions.update(dict([R for R in reactions.items() 
                                          if R[1]['constant'] >= self._min_K
                                          or R[1]['type'] == 'A']))
    #end def
    
    
    def _register_added_duplexes(self, hit, _pos, _duplexes, _template):
        self._annealings[hit].append((_pos, _duplexes, _template))
        self._add_template(hit, _template)
        alt_annealings  = (_pos, [], _template)
        for _duplex, _id in _duplexes:
            self._used_primer_ids.add(_id)
            #if no exonuclease, try to add alternative annealing without 3' matches to side_annealings
            if not self._with_exonuclease and not _duplex.fwd_3_mismatch:
                #try to add alternative annealing without 3' matches to side annealings
                new_duplex = deepcopy(_duplex)
                if new_duplex.strip_3_matches(self._max_3_matches) \
                and new_duplex.K >= self._min_K:
                    alt_annealings[1].append((new_duplex, _id))
        if alt_annealings[1]: self._annealings[hit].append(alt_annealings)
    #end def
    
    
    def _register_side_annealings(self, hit, sorted_annealings):
        for _pos, _duplexes, _template in sorted_annealings:
            registered_annealings = (_pos, [], _template)
            for _duplex, _id in _duplexes:
                if _id in self._used_primer_ids:
                    registered_annealings[1].append((_duplex, _id))
            if registered_annealings[1]: 
                self._annealings[hit].append(registered_annealings)
                self._add_template(hit, _template)
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
                #check equilibrium constant
                if _duplex.K < self._min_K: continue
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
                elif not _duplex.fwd_3_mismatch: #if not, check 3' annealing
                    good_annealings[1].append((_duplex, _id))
                else:
                    bad_annealings[1].append((_duplex, _id))
            #add annealings
            if good_annealings[1]: good[strand].append(good_annealings)
            if bad_annealings[1]:  bad.append(bad_annealings)
    #end def
    
    
    def add_annealings(self, hit, fwd_annealings, rev_annealings):
        '''Add primer annealing sites to the simulation. Within them find 
        PCR products, add other annealings as side reactions.
        hit - name of a target sequence
        fwd_annealings - a list of primer annealing sites on a direct strand 
        of the target sequence: [(position, [(duplex, primer_id), ...]), ...]
        rev_annealings - a list of primer annealing sites on a reverse strand 
        of the target sequence (the structure is the same).
        Return True if some products produced by the annealings were added. 
        False otherwise.'''
        if not fwd_annealings or not rev_annealings: return False
        self._annealings[hit] = []
        good_annealings       = ([], []) #0 is reverse strand, 1 is forward
        bad_annealings        = []
        #sort annealings into good (ones that are suitable for product generation) 
        #and bad (ones that are not) 
        self._sort_annealings(hit, 1, fwd_annealings, good_annealings, bad_annealings)
        self._sort_annealings(hit, 0, rev_annealings, good_annealings, bad_annealings)
        #sort good annealings by position
        good_annealings[1].sort(key=lambda x: x[0])
        good_annealings[0].sort(key=lambda x: x[0])
        #find possible products in range [min_amplicon, max_amplicon]
        products_added  = False
        added_positons  = (set(), set())
        rev_start       = 0
        for fwd_pos, fwd_dups, fwd_templ in good_annealings[1]:
            start = fwd_pos+1
            for rev_pi, (rev_pos, rev_dups, rev_templ) in enumerate(good_annealings[0][rev_start:]):
                end = rev_pos-1
                if start >= end:
                    rev_start = rev_pi
                    continue
                if not (self._min_amplicon <= end-start+1 <= self._max_amplicon): 
                    break
                self._add_product(hit, start, end, fwd_dups, rev_dups)
                if fwd_pos not in added_positons[1]:
                    self._register_added_duplexes(hit, fwd_pos, fwd_dups, fwd_templ)
                if rev_pos not in added_positons[0]:
                    self._register_added_duplexes(hit, rev_pos, rev_dups, rev_templ)
                if not products_added: products_added = True
        if not products_added: #delete all annealings if no products have been added
            del self._annealings[hit]
        else: 
            self._nonzero = False
            if self._include_side_annealings:
                self._register_side_annealings(hit, [ann for ann in good_annealings[1]
                                                     if ann[0] not in added_positons[1]])
                self._register_side_annealings(hit, [ann for ann in good_annealings[0]
                                                     if ann[0] not in added_positons[0]])
                self._register_side_annealings(hit, bad_annealings)
        return products_added
    #end def
    
    
    def _add_product(self, 
                       hit,          #name of the target sequence
                       start,        #start position of the product on the target sequence
                       end,          #end position of the product on the target sequence
                       fwd_duplexes, #duplexes of forward primers with corresponding template sequence
                       rev_duplexes, #duplexes of reverse primers with corresponding template sequence
                      ):
        '''
        Add a new product to the simulation.
        hit -- name of template sequence (termed as in BLAST)
        start, end -- positions in the sequence (starting form 1) corresponding 
        to the beginning and the end of PCR prduct
        *_duplexes -- lists of tuples of the form (duplex_object, id), where 
        duplex_object represents annealed primer and id is it's name 
        '''
        #initialize new product and reset nonzero flag
        new_product = Product(hit, start, end, fwd_duplexes, rev_duplexes)
        #if there was no such hit before, add it
        if hit not in self._products:
            self._products[hit] = dict()
        #if no such product from this hit, add it to the list
        new_product_hash = hash(new_product)
        if new_product_hash not in self._products[hit]:
            self._products[hit][new_product_hash]  = new_product
        else: #append primers
            self._products[hit][new_product_hash] += new_product
    #end def
    
    
    def _add_template(self, hit, new_template):
        #add new template
        if hit not in self._templates:
            self._templates[hit] = [[],[]]
        strand = new_template.forward
        self._templates[hit][strand].append(new_template)
        self._templates[hit][strand].sort(key=lambda(x): x.start)
        compacted = [self._templates[hit][strand][0]]
        for T in self._templates[hit][strand]:
            if T.overlaps(compacted[-1]):
                compacted[-1] += T
            else: compacted.append(T)
        self._templates[hit][strand] = compacted
    #end def
    
    
    def _find_template(self, hit, query_template):
        for T in self._templates[hit][query_template.forward]:
            if T.overlaps(query_template):
                return T
        return None
    
    
    def _construct_annealing_reactions(self, hit_id, annealings):
        reactions = dict()
        for _pos, _duplexes, _template in annealings:
            template = self._find_template(hit_id, _template)
            template_hash = hash(template)
            for _duplex, _id in _duplexes:
                r_hash = hash((_duplex.fwd_seq, template))
                reactions[r_hash] = Equilibrium.compose_reaction(_duplex.K, 
                                                                 _duplex.fwd_seq, 
                                                                 template_hash, 'AB')
        return reactions
    #end def
    
    
    def _calculate_equilibrium_for_hit(self, hit_id):
        #assemble reactions and concentrations for this hit
        reactions = self._construct_annealing_reactions(hit_id, self._annealings[hit_id])
        if not reactions: return hit_id, None
        reactions.update(self._side_reactions)
        concentrations = dict(self._primer_concentrations)
        concentrations.update([(hash(T), TD_Functions.PCR_P.DNA) 
                               for T in self._templates[hit_id][1]])
        concentrations.update([(hash(T), TD_Functions.PCR_P.DNA) 
                               for T in self._templates[hit_id][0]])
        concentrations.update(self._side_concentrations)
        #calculate equilibrium
        equilibrium = Equilibrium(reactions, concentrations)
        equilibrium.calculate()
        return hit_id, equilibrium
    #end def
        
            
    def _calculate_equilibrium(self, abort_e):
        results = self._parallelize_work(abort_e, 5e-3, 
                                         self._calculate_equilibrium_for_hit, 
                                         self._products.keys())
        for hit_id, equilibrium in results:
            self._solutions[hit_id]   = equilibrium
            self._max_objective_value = max(self._max_objective_value, 
                                            equilibrium.solution_objective_value)
    #end def
    
    
    def _calculate_quantities_for_hit(self, hit_id):
        equilibrium = self._solutions[hit_id]
        if not equilibrium: return hit_id, None, None
        solution  = equilibrium.solution
        cur_state = [TD_Functions.PCR_P.dNTP*4.0, #a matrix is considered to have equal quantities of each letter
                     dict(self._primer_concentrations), []]
        reaction_end = {'poly':[], 'cycles':3}
        cur_primers,cur_variants = cur_state[1],cur_state[2]
        products = self._products[hit_id]
        for product_id in products:
            product = products[product_id]
            #forward primer annealing, first cycle
            fwd_strands_1 = dict()
            fwd_template  = self._find_template(hit_id, product.fwd_template)
            for fwd_primer, _id in product.fwd_primers:
                r_hash = hash((fwd_primer.fwd_seq, fwd_template))
                if r_hash in solution:
                    fwd_strands_1[fwd_primer.fwd_seq] = [solution[r_hash],
                                                         equilibrium.get_product_concentration(r_hash)]
                    cur_primers[fwd_primer.fwd_seq] -= fwd_strands_1[fwd_primer.fwd_seq][1] 
                    cur_state[0] -= fwd_strands_1[fwd_primer.fwd_seq][1]*len(product)
            #reverse primer annealing, first cycle
            rev_strands_1 = dict()
            rev_template  = self._find_template(hit_id, product.rev_template)
            for rev_primer, _id in product.rev_primers:
                r_hash = hash((rev_primer.fwd_seq, rev_template))
                if r_hash in solution:
                    rev_strands_1[rev_primer.fwd_seq] = [solution[r_hash],
                                                         equilibrium.get_product_concentration(r_hash)]
                    cur_primers[rev_primer.fwd_seq] -= rev_strands_1[rev_primer.fwd_seq][1]
                    cur_state[0] -= rev_strands_1[rev_primer.fwd_seq][1]*len(product)
            #second and third cycles
            for f_id in fwd_strands_1:
                for r_id in rev_strands_1:
                    #second
                    fwd_strand_2  = (fwd_strands_1[f_id][0])*(rev_strands_1[r_id][1])
                    rev_strand_2  = (rev_strands_1[r_id][0])*(fwd_strands_1[f_id][1])
                    cur_primers[f_id] -= fwd_strand_2
                    cur_primers[r_id] -= rev_strand_2
                    cur_state[0] -= (fwd_strand_2+rev_strand_2)*len(product)
                    #third
                    var = [f_id, r_id, fwd_strand_2+rev_strand_2, product_id]
                    cur_variants.append(var)
                    cur_primers[f_id] -= rev_strand_2
                    cur_primers[r_id] -= fwd_strand_2
                    cur_state[0] -= (fwd_strand_2+rev_strand_2)*len(product)
        if not cur_variants: return hit_id, None, None
        #check if primers or dNTP were used in the first three cycles
        if cur_state[0] <= 0 or min(cur_state[1]) <= 0:
            print '\nPCR simulation warning:'
            print '   Template DNA: %s' % hit_id
            if cur_state[0] <= 0:
                print '   dNTP have been depleted in the first three cycles.'
            if min(cur_state[1]) <= 0:
                print '   Some primers have been depleted in the first three cycles.'
            print 'Try to change reaction conditions.'
            return hit_id, None, None 
        #all consequent PCR cycles
        cur_variants.sort(key=lambda(x): x[2])
        for cycle in range(4, self._num_cycles+1):
            idle = True
            #save current state
            prev_state = deepcopy(cur_state)
            #calculate a normal cycle
            for var in cur_variants:
                #if any primer is depleted, continue
                if cur_primers[var[0]] is None \
                or cur_primers[var[1]] is None: 
                    continue
                #count consumption
                cur_primers[var[0]] -= var[2]
                cur_primers[var[1]] -= var[2]
                cur_state[0] -= var[2]*2*len(products[var[3]])
                #count synthesis
                var[2] += var[2]
                #something was synthesized in this cycle
                idle = False
            if idle: break #all primers have been depleted
            self._correct_cycle(products, reaction_end, cycle, prev_state, cur_state)
            #count the cycle of each synthesized products
            for var0, var1 in zip(prev_state[2], cur_state[2]):
                if var0[2] != var1[2]: 
                    products[var1[3]].cycles = cycle
            #count the cycle of the reaction
            reaction_end['cycles'] = cycle
            #if no dNTP left, end the reaction
            if cur_state[0] <= 0: break
        #end cycles
        reaction_end['dNTP'] = cur_state[0]/4.0 if cur_state[0] >= 0 else 0
        #sum up all variants of each product
        for var in cur_variants:
            products[var[3]].quantity += var[2]
        #filter out products with zero quantity
        max_product_quantity   = max(prod[1].quantity 
                                     for prod in products.items())
        products = dict([prod for prod in products.items() 
                        if prod[1].quantity > max(TD_Functions.PCR_P.DNA, 
                                                  max_product_quantity*self._min_quantity_factor)])
        return hit_id, products, reaction_end
    #end def
    
    
    def _calculate_product_quantities(self, abort_e):
        results = self._parallelize_work(abort_e, 1e-3, 
                                         self._calculate_quantities_for_hit, 
                                         self._products.keys())
        for hit_id, products, reaction_end in results:
            self._products[hit_id]      = products
            self._reaction_ends[hit_id] = reaction_end
    #end def
    
    
    def _correct_cycle(self, products, reaction_end, cycle, prev_state, cur_state):
        prev_dNTP,prev_primers,prev_variants = prev_state
        cur_dNTP,cur_primers,cur_variants = cur_state
        #check if: some primers were depleted
        depleted_primers = dict()
        for p_id in cur_primers:
            if  cur_primers[p_id] != None \
            and cur_primers[p_id] < 0:
                consume_ratio = (prev_primers[p_id]/
                                 (prev_primers[p_id] 
                                  - cur_primers[p_id]))
                depleted_primers[p_id] = consume_ratio
        #if so, recalculate consumption for products which use them
        if depleted_primers:
            for var0, var1 in zip(prev_variants, cur_variants):
                if  not var1[0] in  depleted_primers \
                and not var1[1] in  depleted_primers:
                    continue
                if cur_primers[var1[0]] is None \
                or cur_primers[var1[1]] is None: 
                    continue
                consume_ratio0 = depleted_primers[var1[0]] if var1[0] in depleted_primers else 1
                consume_ratio1 = depleted_primers[var1[1]] if var1[1] in depleted_primers else 1
                real_addition = (var1[2]-var0[2])*min(consume_ratio0,
                                                      consume_ratio1)
                #restore subtracted
                cur_primers[var1[0]] += var0[2]
                cur_primers[var1[1]] += var0[2]
                cur_dNTP += var0[2]*2*len(products[var1[3]])
                #subtract real consumption
                cur_primers[var1[0]] -= real_addition
                cur_primers[var1[1]] -= real_addition
                cur_dNTP -= real_addition*2*len(products[var1[3]])
                #count synthesis
                var1[2] = var0[2] + real_addition
        #set depleted primers to exact zero (due to floating point errors they may still be small real numbers)
        for p_id in depleted_primers: cur_primers[p_id] = 0
        #check if: dNTP was used,
        #     and  polymerase activity could handle this cycle.                    
        dNTP_consumption = prev_dNTP - cur_dNTP
        max_consumptioin = self._polymerase*10e-9*self._elongation_time/30
        #save polymerase shortage information
        if dNTP_consumption > max_consumptioin:
            if reaction_end['poly']:
                #if previous cycle passed without a shortage, start new record
                if cycle - reaction_end['poly'][-1][1] > 1:
                    reaction_end['poly'].append([cycle, cycle])
                else: #append to the current record 
                    reaction_end['poly'][-1][1] = cycle
            else: #if there was no shortage, add current new record
                reaction_end['poly'].append([cycle, cycle])
        #correct consumption if needed
        if dNTP_consumption > max_consumptioin \
        or cur_dNTP < 0:
            consume_ratio = min(max_consumptioin/dNTP_consumption,
                                prev_dNTP/dNTP_consumption)
            cur_dNTP = prev_dNTP
            cur_primers.update(prev_primers)
            for var0, var1 in zip(prev_variants, cur_variants):
                if cur_primers[var1[0]] is None \
                or cur_primers[var1[1]] is None: 
                    continue
                #due to primer correction var1[2]-var0[2] is not always equal to var0[2]
                real_addition = (var1[2]-var0[2])*consume_ratio
                cur_primers[var1[0]] -= real_addition  
                cur_primers[var1[1]] -= real_addition
                cur_dNTP -= real_addition*2*len(products[var1[3]])
                var1[2] = var0[2] + real_addition
        cur_state[0] = cur_dNTP
        #check again if some primers were depleted
        for p_id in cur_primers:
            if cur_primers[p_id] <= 0:
                cur_primers[p_id] = None
    #end def


    def run(self):
        if not self._products: return
        abort_event = self._new_event()
        #calculate equilibrium in the system
        self._calculate_equilibrium(abort_event)
        if abort_event.is_set(): return
        self._calculate_product_quantities(abort_event)
        if abort_event.is_set(): return
        self._clean_abort_even(abort_event)

        #filter out hits without products
        self._products = dict([hit for hit in self._products.items() 
                               if hit[1] is not None])
        self._nonzero  = len(self._products) > 0
    #end def    
    
    
    @classmethod
    def _construct_histogram(cls, main_title, products, reaction_ends=None, with_titles=True, sort=True):
        #construct histogram
        histogram  = []
        max_start  = max(len(str(p.start)) for p in products)
        max_end    = max(len(str(p.end))   for p in products)
        max_len    = max(len(str(len(p)))  for p in products)
        for product in products:
            product_spec = ''
            #if the flag is set, print full hit title
            if with_titles:
                product_spec += product.name + '\n'
            #compose product specification
            product_spec  += '%d%s bp [%d%s-%s%d]' \
                           % (len(product),
                              ' '*(max_len-len(str(len(product)))),
                              product.start,
                              ' '*(max_start-len(str(product.start))),
                              ' '*(max_end-len(str(product.end))),
                              product.end)
            #print last cycle of the reaction where this product was still generating
            if reaction_ends:
                product_spec += ' %d cycles' % product.cycles
            product_spec += '\n'
            #print primers which give this product
            if not with_titles:
                primer_ids  = ', '.join(product.fwd_ids)
                if primer_ids: primer_ids += ' <> '
                primer_ids += ', '.join(product.rev_ids)
                if primer_ids: product_spec += primer_ids + '\n'
            #append histogram row
            histogram.append([product_spec, product.quantity, len(product)])
        #sort histogram by product length, from the longest to the shortest
        if sort: histogram.sort(key=lambda(x): x[2], reverse=True)
        #format histogram
        return cls._format_histogram(main_title, histogram)
    #end def
    
    
    @classmethod
    def _format_histogram(cls, title, histogram):
        #maximum column name width and value
        max_value = max(c[1] for c in histogram)
        widths    = [StringTools.text_width-cls._hist_width-1, cls._hist_width+1]
        histogram_string  = '-'*StringTools.text_width + '\n'
        #wrap column and hist titles
        histogram_string += line_by_line([title, cls._hist_column_title], 
                                         widths, j_center=True)
        histogram_string += '-'*StringTools.text_width + '\n'
        #histogram lines
        for col in histogram:
            #line value
            hist_value = int(round((cls._hist_width*col[1])/max_value))
            col_spacer = cls._hist_width - hist_value
            #value figure
            value_str  = TD_Functions.format_concentration(col[1])
            hist_bar  = ''
            if len(value_str) < col_spacer:
                _spacer = col_spacer-len(value_str)
                hist_bar += '#'*hist_value + ' '*_spacer + value_str 
            else:
                _bar = hist_value-len(value_str)
                hist_bar += '#'*(_bar/2) + value_str 
                hist_bar += '#'*(_bar-_bar/2) + ' '*col_spacer
            histogram_string += line_by_line([col[0],hist_bar], widths, divider=':')
            histogram_string += '-'*StringTools.text_width + '\n'
        histogram_string += '\n'
        return histogram_string
    #end def
    
    
    def _construct_electrophoresis(self, products):
        max_len      = max(len(p) for p in products)
        window       = int(max_len*self._window_percent)
        nearest_srip = max_len-window
        max_len_log  = int(log(max_len)*self._precision)
        min_len_log  = int(log(min(len(p) for p in products))*self._precision)
        window_log   = int((log(max_len)-log(nearest_srip))*self._precision)
        #construct phoresis
        phoresis    = [[l,0,int(exp(l/self._precision))] for l in range(min_len_log, max_len_log+window_log, window_log)]
        for product in products:
            l = int(log(len(product))*self._precision)
            p = (l - min_len_log)/window_log
            #fluorescence intensity is proportional to the length of a duplex
            phoresis[p][1] += product.quantity*len(product)
        phoresis = phoresis[::-1]
        #format electrophorogram
        return self._format_electrophoresis(window, phoresis)
    #end def
    
    
    def _format_electrophoresis(self, window, phoresis):
        text_width  = StringTools.text_width
        max_line    = max(min(p.concentration for p in self._primers)*self._max_amplicon, 
                          max(l[1] for l in phoresis))
        max_mark    = max(len(str(l[2]+window)) for l in phoresis)*2
        line_width  = text_width - max_mark - 7 #mark b :###   :
        #format phoresis
        phoresis_text = ''
        phoresis_text += ' '*(max_mark+5)+':'+'-'*line_width+':\n'
        phoresis_text += ' '*(max_mark+5)+':'+' '*line_width+':\n'
        for l in range(len(phoresis)):
            line = phoresis[l]
            prev_mark   = phoresis[l-1][2] if l > 0 else line[2]+window
            mark_spacer = max_mark - len(str(line[2])) - len(str(prev_mark))
            line_value  = int(round((line_width*line[1])/max_line))
            line_spacer = line_width - line_value 
            phoresis_text += '%d-%d%s bp :%s%s:\n' % (prev_mark,
                                                      line[2], 
                                                      ' '*mark_spacer,
                                                      '#'*line_value,
                                                      ' '*line_spacer) 
        phoresis_text += ' '*(max_mark+5)+':'+' '*line_width+':\n'
        phoresis_text += ' '*(max_mark+5)+':'+'-'*line_width+':\n'
        return phoresis_text
    #end def
    

    def per_hit_header(self, hit):
        if not self._nonzero: 
            return '\nNo PCR products have been found.\n'
        header_string  = ''
        header_string += 'Reaction ended in: %d cycles.\n' \
            % self._reaction_ends[hit]['cycles']
        header_string += 'C(dNTP) after reaction: ~%s\n' \
            % TD_Functions.format_concentration(self._reaction_ends[hit]['dNTP'])
        if self._reaction_ends[hit]['cycles'] < self._num_cycles:
            if self._reaction_ends[hit]['dNTP'] == 0:
                header_string += 'dNTP have been depleted.\n'
            else:
                header_string += 'Primers have been depleted.\n'
        if self._reaction_ends[hit]['poly']:
            header_string += wrap_text('Polymerase activity was '
                                       'insufficient in these cycles:\n')
            shortage_string = ''
            for shortage_period in self._reaction_ends[hit]['poly']:
                if shortage_string: shortage_string += ', '
                if shortage_period[0] == shortage_period[1]:
                    shortage_string += str(shortage_period[0])
                else:
                    shortage_string += '%d-%d' % tuple(shortage_period)
            header_string += wrap_text(shortage_string + '\n')
        header_string += '\n'
        return header_string

        
    def all_products_histogram(self):
        if not self._nonzero: 
            return '\nNo PCR products have been found.\n'
        all_products = []
        prods_list = self._products.values()
        prods_list.sort(key=lambda(d): max(p.quantity for p in d.values()), reverse=True)
        for p_dict in prods_list:
            prods = p_dict.values()
            prods.sort(key=len, reverse=True)
            all_products += prods
        return self._construct_histogram(self._all_products_title, all_products, 
                                         None, with_titles=True, sort=False)
    #end def
    
    
    def per_hit_histogram(self, hit):
        if not self._nonzero: 
            return '\nNo PCR products have been found.\n'
        products = self._products[hit].values()
        return self._construct_histogram('%s' % hit, products, 
                                         self._reaction_ends, with_titles=False)
    #end def
    
    
    def per_hit_electrophoresis(self, hit):
        if not self._nonzero: 
            return '\nNo PCR products have been found.\n'
        products = self._products[hit].values()
        return self._construct_electrophoresis(products)
    #end def
    
    
    def all_graphs_grouped_by_hit(self):
        if not self._nonzero: 
            return '\nNo PCR products have been found.\n'
        all_graphs = ''
        prods = self._products.items()
        prods.sort(key=lambda(d): max(p.quantity for p in d[1].values()), reverse=True)
        hits = [p[0] for p in prods]
        for hit in hits:
            all_graphs += self.per_hit_histogram(hit)
            all_graphs += self.per_hit_header(hit)
            all_graphs += hr(' electrophorogram of PCR products ')
            all_graphs += self.per_hit_electrophoresis(hit)
            all_graphs += hr('')
            all_graphs += '\n\n'
        return all_graphs
    
    
    def format_report_header(self):
        header_string  = ''
        header_string += hr(' PCR conditions ')
        if self._with_exonuclease:
            header_string += "DNA polymerase HAS 3'-5'-exonuclease activity.\n"
        else: header_string += "DNA polymerase doesn't have 3'-5'-exonuclease activity.\n"
        header_string += 'Speed of DNA polymerase is considered to be 1 kbs/min\n\n'
        header_string += 'Minimum amplicon size: %d\n' % self._min_amplicon
        header_string += 'Maximum amplicon size: %d\n' % self._max_amplicon
        header_string += 'Elongation time:       %s\n' % timedelta(minutes=self._elongation_time)
        header_string += 'Maximum cycles:        %d\n' % self._num_cycles
        header_string += '\n'
        header_string += format_PCR_conditions(self._primers, self._polymerase)+'\n'
        header_string += hr(' primers and their melting temperatures ')
        for primer in self._primers:
            header_string += repr(primer) + '\n'
        header_string += '\n'
        return header_string
    #end def
    
    
    def format_quantity_explanation(self):
        expl_string  = ''
        expl_string += hr(' estimation of PCR products concentrations ')
        expl_string += 'Value of an objective function at the solution ' + \
                       '(the lower the better):\n   %e\n' % \
                        self._max_objective_value
        expl_string += wrap_text('This value shows "distance" to the solution of '
                           'the system of equilibrium equations which were used '
                           'to calculate concentrations of PCR products.\n\n')
        expl_string += wrap_text(('Products with concentration less than %.2f%% '
                                  'of the concentration of the most abundant '
                                  'product or less than initial DNA concentration '
                                  'are not shown.'
                                 '\n\n') % (self._min_quantity_factor*100))
        expl_string += wrap_text('Boundaries of a product and it\'s length do '
                                 'not include primers.\n\n')
        expl_string += '\n'
        return expl_string
    #end def
    
    
    def format_products_report(self):
        prod_string  = ''
        prod_string += wrap_text('For each target sequence a list of possible '
                                 'PCR products which were predicted by the '
                                 'simulation is given. ' 
                                 'Information about each product includes:\n'
                                 '-start and end positions on a target sequence.\n'
                                 '   *the first nucleotide of a sequence is at the position 1\n'
                                 '   *primers are not included, so their 3\'-ends are located at '
                                 'start-1 and end+1 positions.\n'
                                 '-length in base pairs.\n'
                                 '   *the length too do not include primers.\n'
                                 '-concentration which was calculated during PCR simulation.\n'
                                 '-number of PCR cycles during which the product has been generated.\n'
                                 '   *if this number is lower than the number of simulated reaction cycles '
                                 'it means that primers producing the product were depleted.\n'
                                 '-lists of forward and reverse primers which produced the product'
                                 '\n\n\n')
        for hit in self._products:
            prod_string += hr(' %s ' % hit, '*')
            products = self._products[hit].values()
            products.sort(key=lambda x: x.start)
            for pi, product in enumerate(products):
                prod_string += hr(' product %d ' % (pi+1), '=')
                prod_string += product.pretty_print(with_name=False)
                prod_string += hr('', '=')
            prod_string += hr('', '*')
        prod_string += '\n'
        return prod_string
#end class



#test
if __name__ == '__main__':
    import StringTools
    from SecStructures import Duplex
    from Primer import Primer
    import os
    os.chdir('../')
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC
    TD_Functions.PCR_P.dNTP  = 0.02e-3
    TD_Functions.PCR_P.DNA   = 1e-10
    TD_Functions.PCR_P.PCR_T = 60
    primers = []
    primers.append(Primer.from_sequences((SeqRecord(Seq('ATATTCTACGACGGCTATCC', IUPAC.unambiguous_dna), id='F-primer'),
                                          SeqRecord(Seq('ATATGCTACGACGGCTATCC', IUPAC.unambiguous_dna)),
                                          SeqRecord(Seq('ATATTCTACGACGGCTATCC', IUPAC.unambiguous_dna))), 
                                         0.1e-6))
    primers.append(Primer.from_sequences((SeqRecord(Seq('CAAGGGCTAGAAGCGGAAG'[::-1], IUPAC.unambiguous_dna), id='R-primer'),
                                          SeqRecord(Seq('CAAGGGCTAGACGCGGAAG'[::-1], IUPAC.unambiguous_dna)),
                                          SeqRecord(Seq('CAAGGCCTAGAGGCGGAAG'[::-1], IUPAC.unambiguous_dna))),
                                         0.1e-6))
    PCR_res = PCR_Simulation(primers, 500, 3000, 0.01*1e6, False, 30)
    PCR_res.add_product('Pyrococcus yayanosii CH1', 982243, 983644,
                        Duplex(Seq('ATATTCTACGACGGCTATCC'), 
                               Seq('ATATTCTACGACGGCTATCC').reverse_complement()),
                        Duplex(Seq('CAAGGGCTAGAAGCGGAAG'[::-1]), 
                               Seq('CAAGGGCTAGAAGCGGAAG').complement())) 
    PCR_res.add_product('Pyrococcus yayanosii CH1', 962243, 963744,
                        Duplex(Seq('ATATTCTACGACGGCTATCC'), 
                               Seq('ATATTCTACGACGGCTATCC').reverse_complement()),
                        Duplex(Seq('CAAGGGCTAGAAGCGGAAG'[::-1]), 
                               Seq('CAAGGGCTAGAAGCGGAAG').complement()))
    PCR_res.add_product('Pyrococcus yayanosii CH2', 962243, 963644, 
                        Duplex(Seq('ATATGCTACGACGGCTATCC'), 
                               Seq('ATATGCTACGACGGCTATCC').reverse_complement()),
                        Duplex(Seq('CAAGGGCTAGACGCGGAAG'[::-1]), 
                               Seq('CAAGGGCTAGACGCGGAAG').complement()))
    PCR_res.add_product('Pyrococcus yayanosii CH3', 922243, 923644, 
                        Duplex(Seq('ATATTCTACGACGGCTATCC'), 
                               Seq('ATATTCTACGACGGCTATCC').reverse_complement()),
                        Duplex(Seq('CAAGGCCTAGAGGCGGAAG'[::-1]), 
                               Seq('CAAGGCCTAGAGGCGGAAG').complement()))

    PCR_res.run()
    print PCR_res.format_report_header()
    print PCR_res.format_quantity_explanation()
    print PCR_res.all_products_histogram()
    print ''
    print PCR_res.all_graphs_grouped_by_hit()
