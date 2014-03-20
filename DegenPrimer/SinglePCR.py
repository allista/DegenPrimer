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
Created on Feb 24, 2014

@author: Allis Tauri <allista@gmail.com>
'''


from abc import ABCMeta
from copy import deepcopy
from Equilibrium import Equilibrium, Reaction
from AbortableBase import AbortableBase
from MultiprocessingBase import MultiprocessingBase
from WorkCounter import WorkCounter
import TD_Functions as tdf


class PCR_Base(AbortableBase):
    '''Base class for PCR simulations'''
    __metaclass__ = ABCMeta
    
    #products with quantity less than maximum quantity multiplied by this factor will not be included in the report
    _min_quantity_factor = 0.001
    
    def __init__(self, 
                  abort_event,
                  polymerase,
                  with_exonuclease,
                  num_cycles,
                  ):
        AbortableBase.__init__(self, abort_event)
        self._polymerase       = polymerase
        self._with_exonuclease = with_exonuclease
        self._num_cycles       = num_cycles
        self._PCR_P            = deepcopy(tdf.PCR_P)
    #end def
#end class


class SinglePCR(PCR_Base, MultiprocessingBase):
    '''Simulates a single PCR'''
    def __init__(self, 
                  abort_event,
                  primers_concentrations,
                  elongation_time,
                  polymerase,
                  with_exonuclease,
                  num_cycles,
                  side_reactions,
                  side_concentrations,
                  ):
        MultiprocessingBase.__init__(self, abort_event)
        PCR_Base.__init__(self, abort_event, polymerase, with_exonuclease, num_cycles)
        self._elongation_time        = elongation_time
        self._primers_concentrations = primers_concentrations
        self._side_reactions         = side_reactions
        self._side_concentrations    = side_concentrations
    #end def
    
    
    @staticmethod
    def _find_template(query_template, templates):
        for T in templates[query_template.forward]:
            if T.overlaps(query_template): return T
        return None
    
    
    @staticmethod
    @MultiprocessingBase.data_mapper
    def _construct_annealing_reaction(annealing, templates):
        reactions = dict()
        _duplexes, _template = annealing
        template = SinglePCR._find_template(_template, templates)
        template_hash = hash(template)
        for _duplex, _id in _duplexes:
            for _dimer in _duplex.dimers:
                r_hash = hash((_duplex.fwd_seq, _dimer, template))
                reactions[r_hash] = Reaction(_dimer.K, 
                                             _duplex.fwd_seq, 
                                             template_hash, 'AB')
        return reactions
    #end def
    
    @staticmethod
    @MultiprocessingBase.results_assembler
    def _reactions_assembler(index, result, reactions):
        reactions.update(result)
    
    def _construct_annealing_reactions(self, pcr_mixture, counter):
        reactions = dict()
        work = self.Work(timeout=0.1, counter=counter)
        work.prepare_jobs(self._construct_annealing_reaction, 
                          pcr_mixture.annealings, None, 
                          pcr_mixture.templates)
        work.set_assembler(self._reactions_assembler, reactions)
        work.start()
        if not work.wait(): return None
        return reactions
    #end def
    
    
    def _calculate_equilibrium(self, pcr_mixture, counter):
        counter.set_subwork(2, (1, 10))
        #assemble reactions and concentrations for this hit
        reactions = self._construct_annealing_reactions(pcr_mixture, counter[0])
        if not reactions: return None
        reactions.update(self._side_reactions)
        concentrations = dict(self._primers_concentrations)
        concentrations.update((hash(T), self._PCR_P.DNA) 
                              for T in pcr_mixture.templates[1])
        concentrations.update((hash(T), self._PCR_P.DNA) 
                              for T in pcr_mixture.templates[0])
        concentrations.update(self._side_concentrations)
        #calculate equilibrium
        equilibrium = Equilibrium(self._abort_event, reactions, concentrations)
        equilibrium.calculate(counter[1])
        return equilibrium
    #end def
    

    def _calculate_first_cycle(self, primers, template, product_len, 
                                  equilibrium, cur_state):
        strand_products = dict()
        for primer, _id in primers:
            strand_products[primer.fwd_seq] = [0, 0]
            for dimer in primer.dimers:
                if dimer.fwd_mismatch and not self._with_exonuclease: continue
                r_hash = hash((primer.fwd_seq, dimer, template))
                if r_hash in equilibrium.solution:
                    strand_products[primer.fwd_seq][0] += equilibrium.solution[r_hash]
                    strand_products[primer.fwd_seq][1] += equilibrium.get_product_concentration(r_hash)
            if strand_products[primer.fwd_seq][1] != 0:
                cur_state[0] -= strand_products[primer.fwd_seq][1]*product_len
                cur_state[1][primer.fwd_seq] -= strand_products[primer.fwd_seq][1] 
            else: del strand_products[primer.fwd_seq]
        return strand_products
    #end def
    
    
    def _calculate_quantities(self, pcr_mixture, equilibrium, counter):
        if not equilibrium: return None, None
        cur_state = [self._PCR_P.dNTP*4.0, #a matrix is considered to have equal quantities of each letter
                     dict(self._primers_concentrations), []]
        reaction_end = {'poly':[], 'cycles':3}
        cur_primers,cur_variants = cur_state[1],cur_state[2]
        products = pcr_mixture.products
        counter.set_subwork(2, (2, self._num_cycles-2))
        #first two cycles
        counter[0].set_work(len(products))
        for product_id, product in products.iteritems():
            product_len = len(product)
            #forward primer annealing, first cycle
            fwd_template  = self._find_template(product.fwd_template, pcr_mixture.templates)
            fwd_strands_1 = self._calculate_first_cycle(product.fwd_primers, 
                                                        fwd_template, product_len, 
                                                        equilibrium, cur_state)
            #reverse primer annealing, first cycle
            rev_template  = self._find_template(product.rev_template, pcr_mixture.templates)
            rev_strands_1 = self._calculate_first_cycle(product.rev_primers, 
                                                        rev_template, product_len, 
                                                        equilibrium, cur_state)
            #second
            for f_id in fwd_strands_1:
                for r_id in rev_strands_1:
                    #second
                    fwd_strand_2  = (fwd_strands_1[f_id][0])*(rev_strands_1[r_id][1])
                    rev_strand_2  = (rev_strands_1[r_id][0])*(fwd_strands_1[f_id][1])
                    cur_primers[f_id] -= fwd_strand_2
                    cur_primers[r_id] -= rev_strand_2
                    cur_state[0] -= (fwd_strand_2+rev_strand_2)*product_len
                    #at this point concentrations of both strands of a product are equal
                    var = [f_id, r_id, fwd_strand_2, product_id] #a variant of *dimeric* product
                    cur_variants.append(var)
            counter[0].count()
        if not cur_variants: return None, None
        #check if primers or dNTP were used in the first three cycles
        if cur_state[0] <= 0 or min(cur_state[1]) <= 0:
            print '\nPCR simulation warning:'
            print '   Template DNA: %s' % pcr_mixture.reaction_id
            if cur_state[0] <= 0:
                print '   dNTP have been depleted in the first two cycles.'
            if min(cur_state[1]) <= 0:
                print '   Some primers have been depleted in the first two cycles.'
            print 'Try to change reaction pcr_mixture.'
            return None, None 
        #all consequent PCR cycles
        counter[1].set_work(self._num_cycles)
        cur_variants.sort(key=lambda(x): x[2])
        for cycle in xrange(3, self._num_cycles+1):
            if self.aborted(): break
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
            counter[1].count()
        #end cycles
        reaction_end['dNTP'] = cur_state[0]/4.0 if cur_state[0] >= 0 else 0
        #sum up all variants of each product
        for var in cur_variants:
            products[var[3]].quantity += var[2]
        #filter out products with zero quantity
        max_product_quantity = max(prod.quantity for prod in products.values())
        products = dict([prod for prod in products.items() 
                        if prod[1].quantity > max(self._PCR_P.DNA, 
                                                  max_product_quantity*self._min_quantity_factor)])
        counter[1].done()
        return products, reaction_end
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

    
    def __call__(self, pcr_mixture, counter=None):
        if not counter: counter = WorkCounter() 
        counter.set_subwork(2, (len(pcr_mixture.annealings), 
                                len(pcr_mixture.products)*self._num_cycles/2.0))
        equilibrium = self._calculate_equilibrium(pcr_mixture, counter[0])
        if equilibrium is None: return None
        products, reaction_end = self._calculate_quantities(pcr_mixture, equilibrium, counter[1])
        if products is None: return None
        return pcr_mixture.reaction_id, equilibrium.objective_value, products, reaction_end
    #end def
#end class