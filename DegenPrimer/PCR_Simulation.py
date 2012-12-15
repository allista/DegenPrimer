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
Created on Nov 16, 2012

@author: Allis Tauri <allista@gmail.com>
'''

from copy import deepcopy
from datetime import timedelta
from math import log, exp
from itertools import chain
from Equilibrium import Equilibrium
import StringTools
import TD_Functions
from TD_Functions import format_PCR_conditions
from StringTools import wrap_text, hr



class Region(object):
    '''Region of a sequence.'''
    def __init__(self, name, start, end):
        self._name  = name
        self._start = min(start, end)
        self._end   = max(start, end)
    #end def
    
    def __repr__(self):
        rep  = '\n'
        rep += 'name:   %s\n' % self._name
        rep += 'start:  %d\n' % self._start
        rep += 'end:    %d\n' % self._end
        rep += 'length: %d\n' % len(self)
        return rep
    #end def
    
    @property
    def name(self): return self._name
    
    @property
    def start(self): return self._start
    
    @property
    def end(self): return self._end

    def __len__(self): return self._end-self._start
        
    def __hash__(self):
        return hash((self._name, self._start, self._end))
        
    def __iadd__(self, T):
        if self._name != T.name: return self
        self._start = min(self._start, T.start)
        self._end   = max(self._end, T.end)
        return self
    #end def
    
    def overlaps(self, T):
        return (self._name == T.name 
                and
                (self._start <= T.start <= self._end)
                or
                (self._start <= T.end <= self._end))
    #end def
#end class


class Product(Region):
    '''Representation of a PCR product. A product is defined as a sequence 
    region bounded by start and end positions which correspond to 3' ends of 
    primers (forward and reverse) which produce this product.'''
    def __init__(self, template_name, start, end, 
                 fwd_primer_duplex=None, rev_primer_duplex=None):
        Region.__init__(self, template_name, start, end)
        #quantity
        self._quantity     = 0
        #primers and templates
        self._fwd_primers  = set()
        self._rev_primers  = set()
        self._fwd_template = Region(template_name, start, start)
        self._rev_template = Region(template_name, end, end)
        self.add_fwd_primer(fwd_primer_duplex)
        self.add_rev_primer(rev_primer_duplex)
    #end def
    
    def __repr__(self):
        rep  = Region.__repr__(self)
        rep += 'quantity: %s\n' % TD_Functions.format_concentration(self._quantity)
        rep += '------\n'
        rep += 'forward template: %s\n' % repr(self._fwd_template)
        rep += 'reverse template: %s\n' % repr(self._rev_template)
        rep += '------\n'
        rep += 'forward primers:\n%s\n' % repr(self._fwd_primers)
        rep += 'reverse primers:\n%s\n' % repr(self._rev_primers)
        return rep
    #end def
    
    @property
    def quantity(self): return self._quantity
    
    @quantity.setter
    def quantity(self, q): self._quantity = q
        
    @property
    def fwd_primers(self): return self._fwd_primers
    
    @property
    def rev_primers(self): return self._rev_primers
    
    @property
    def fwd_template(self): return self._fwd_template
    
    @property
    def rev_template(self): return self._rev_template


    def add_fwd_primer(self, primer_duplex):
        if not primer_duplex: return
        self._fwd_primers.add(primer_duplex)
        self._fwd_template += Region(self._name, 
                                     max(self._start-len(primer_duplex.fwd_seq), 0), 
                                     self._start)
    #end def
        
    def add_rev_primer(self, primer_duplex):
        if not primer_duplex: return
        self._rev_primers.add(primer_duplex)
        self._rev_template += Region(self._name, 
                                     self._end, 
                                     self._end+len(primer_duplex.fwd_seq))
    #end def
#end class
    


class PCR_Simulation(object):
    '''In silica PCR results filtration, computation and visualization.'''
    
    #histogram
    _all_products_title = 'PCR-product'
    _hist_column_title  = '-----final concentration-----'
    #bar width should not be less than 2*len('100.00%')+1, or value will not fit
    _hist_width         = max(2*len('999.9 mM')+2, len(_hist_column_title)) 
    #electrophoresis
    _window_percent     = 0.05 #percentage of a length of the longest PCR product; it defines band width on electrophoresis
    #indirectly defines to what extent should lengths 
    #of two products differ in order of them to occupy 
    #separate bands on the electrophoresis.
    _precision          = 1e6   


    def __init__(self, 
                 primers,                #all primers (generally degenerate) that may be present in the system
                 min_amplicon,           #minimum amplicon length 
                 max_amplicon,           #maximum amplicon length
                 polymerase,             #polymerase concentration in Units per liter
                 with_exonuclease=False, #if polymerase does have have 3' exonuclease activity, products of primers with 3'-mismatches will also be icluded
                 num_cycles=20,          #number of PCR cycles
                 ):
        self._primers             = primers
        self._min_amplicon        = min_amplicon
        self._max_amplicon        = max_amplicon
        self._polymerase          = polymerase
        self._with_exonuclease    = with_exonuclease
        self._elongation_time     = max_amplicon/1000.0 #minutes
        self._num_cycles          = num_cycles
        self._reaction_ends       = dict()

        self._nonzero   = False  #true if there're some products and their quantities were calculated
        
        self._templates = []     #list of all templates
        self._products  = dict() #dictionary of all possible products
        self._side_reactions = dict() #dictionary of all reactions
        self._side_concentrations   = dict()
        self._primer_concentrations = dict()
        for primer in self._primers:
            self._primer_concentrations.update(dict.fromkeys(primer.str_sequences, primer.concentration))
        
        self._solutions = dict() #solutions of equilibrium systems
        self._max_objective_value = -1 #objective value of the worst solution
    #end def
    
    
    def __nonzero__(self):
        return self._nonzero

    def hits(self):
        return list(self._products.keys())
    
    def add_side_concentrations(self, concentrations):
        if concentrations: self._side_concentrations.update(concentrations)
        
    def add_side_reactions(self, reactions):
        if reactions: self._side_reactions.update(reactions)
    

    def add_product(self, 
                    hit,        #name of the target sequence
                    start,      #start position of the product on the target sequence
                    end,        #end position of the product on the target sequence
                    fwd_duplex, #duplex of forward primer with corresponding template sequence
                    rev_duplex, #duplex of reverse primer with corresponding template sequence
                    ):
        '''add a new product to the list'''
        #check if there are such primers in the system at all
        fwd_valid = False
        rev_valid = False
        for primer in self._primers:
            fwd_valid |= fwd_duplex.fwd_seq in primer
            rev_valid |= rev_duplex.fwd_seq in primer
        if not (fwd_valid and rev_valid): return False
        #check 3' mismatches
        if not self._with_exonuclease \
        and (fwd_duplex.fwd_3_mismatch or rev_duplex.fwd_3_mismatch): 
            return False
        #construct product
        new_product = Product(hit, start, end, fwd_duplex, rev_duplex)
        #check amplicon length
        if not (self._min_amplicon <= len(new_product) <= self._max_amplicon): 
            return False
        #if new product is to be added, reset nonzero flag
        self._nonzero = False
        #if there was no such hit before, add it
        if hit not in self._products.keys():
            self._products[hit] = dict()
        #if no such product from this hit, add it to the list
        new_product_hash = hash(new_product)
        if new_product_hash not in self._products[hit].keys():
            self._products[hit][new_product_hash]  = new_product
        else: #append primers
            self._products[hit][new_product_hash] += new_product
        #add templates
        self._add_template(self._products[hit][new_product_hash].fwd_template)
        self._add_template(self._products[hit][new_product_hash].rev_template)
        return True
    #end def
    
    
    def _add_template(self, new_template):
        #add new template
        self._templates.append(new_template)
        self._templates.sort(key=lambda(x): x.start)
        #compact overlapping templates
        compacted = [self._templates[0]]
        for T in self._templates:
            if T.overlaps(compacted[-1]):
                compacted[-1] += T
            else: compacted.append(T)
        self._templates = compacted
    #end def
    
    
    def _find_template(self, query_template):
        for T in self._templates:
            if T.overlaps(query_template):
                return T
        return None
    
    
    def _construct_reactions(self, hit_ids):
        reactions = dict()
        for hit_id in hit_ids:
            for product in self._products[hit_id].values():
                fwd_template = self._find_template(product.fwd_template)
                for fwd_primer in product.fwd_primers:
                    r_hash = hash((fwd_primer.fwd_seq, fwd_template))
                    reactions[r_hash] = Equilibrium.compose_reaction(fwd_primer.K, 
                                                                     fwd_primer.fwd_seq, 
                                                                     hash(fwd_template), 'AB')
                rev_template = self._find_template(product.fwd_template)
                for rev_primer in product.rev_primers:
                    r_hash = hash((rev_primer.fwd_seq, rev_template))
                    reactions[r_hash] = Equilibrium.compose_reaction(rev_primer.K, 
                                                                     rev_primer.fwd_seq, 
                                                                     hash(rev_template), 'AB')
        return reactions
    #end def
    
            
    def _calculate_equilibrium(self):
        #compose a list of reactant concentrations
        concentrations = dict()
        concentrations.update(self._primer_concentrations)
        concentrations.update(dict.fromkeys([hash(T) for T in self._templates], TD_Functions.C_DNA))
        concentrations.update(self._side_concentrations)
        #compose list(s) of reactions and
        #use Equilibrium to compute dimer quantities
        self._max_objective_value = 0
        for hit_id in self._products.keys():
            reactions = self._construct_reactions((hit_id,))
            reactions.update(self._side_reactions)
            equilibrium = Equilibrium(reactions, concentrations)
            equilibrium.calculate()
            self._solutions[hit_id] = equilibrium
            self._max_objective_value = max(self._max_objective_value, 
                                            equilibrium.solution_objective_value)
    #end def
    
    
    def run(self):
        if not self._products: return
        #calculate equilibrium in the system
        self._calculate_equilibrium()
        #compute quantities of products, filter out those with low quantity
        new_hits = dict()
        for hit_id in self._products:
            equilibrium = self._solutions[hit_id]
            solution    = equilibrium.solution
            primer_concentrations = dict(self._primer_concentrations)
            dNTP_concentration    = TD_Functions.C_dNTP*4.0 #a matrix is considered to have equal quantities of each letter
            dNTP_consumption      = 0
            all_variants = []
            new_products = dict()
            for product_id in self._products[hit_id]:
                product = self._products[hit_id][product_id]
                #forward primer annealing, first cycle
                rev_strands_1 = dict()
                fwd_template = self._find_template(product.fwd_template)
                for fwd_primer in product.fwd_primers:
                    r_hash = hash((fwd_primer.fwd_seq, fwd_template))
                    if r_hash in solution:
                        rev_strands_1[fwd_primer.fwd_seq] = [solution[r_hash],
                                                             equilibrium.get_product_concentration(r_hash)]
                        primer_concentrations[fwd_primer.fwd_seq] -= equilibrium.reactants_consumption(fwd_primer.fwd_seq)[0]
                        dNTP_concentration -= equilibrium.reactants_consumption(fwd_primer.fwd_seq)[0]*len(product)
                #reverse primer annealing, first cycle
                fwd_strands_1 = dict()
                rev_template = self._find_template(product.fwd_template)
                for rev_primer in product.rev_primers:
                    r_hash = hash((rev_primer.fwd_seq, rev_template))
                    if r_hash in solution:
                        fwd_strands_1[rev_primer.fwd_seq] = [solution[r_hash],
                                                             equilibrium.get_product_concentration(r_hash)]
                        primer_concentrations[rev_primer.fwd_seq] -= equilibrium.reactants_consumption(rev_primer.fwd_seq)[0]
                        dNTP_concentration -= equilibrium.reactants_consumption(rev_primer.fwd_seq)[0]*len(product)
                #second and third cycles
                for f_id in fwd_strands_1:
                    for r_id in rev_strands_1:
                        #second
                        fwd_strand_2 = (fwd_strands_1[f_id][0])*(rev_strands_1[r_id][1])
                        rev_strand_2 = (rev_strands_1[r_id][0])*(fwd_strands_1[f_id][1])
                        #third
                        var = [f_id, r_id, fwd_strand_2+rev_strand_2, product_id]
                        all_variants.append(var)
                        primer_concentrations[f_id] -= fwd_strand_2+rev_strand_2
                        primer_concentrations[r_id] -= fwd_strand_2+rev_strand_2
                        dNTP_concentration -= (fwd_strand_2+rev_strand_2)*len(product)
                #add variants of this product to the common list
                #and product itself to the new list of products of current hit
                if fwd_strands_1 and rev_strands_1:
                    new_products[product_id] = product
            if not new_products: continue
            #all consequent PCR cycles
            all_variants.sort(key=lambda(x): x[2])
            self._reaction_ends[hit_id] = {'poly':[],
                                           'cycles':3,
                                           'products':dict().fromkeys(new_products.keys(),0)}
            for cycle in range(4, self._num_cycles+1):
                #save current state
                idle = True
                dNTP_previous = dNTP_concentration
                prev_variants = deepcopy(all_variants)
                prev_primers  = deepcopy(primer_concentrations)
                #calculate a normal cycle
                for var in all_variants:
                    #if any primer is depleted, continue
                    if primer_concentrations[var[0]] is None \
                    or primer_concentrations[var[1]] is None: 
                        continue
                    #count consumption
                    primer_concentrations[var[0]] -= var[2]
                    primer_concentrations[var[1]] -= var[2]
                    dNTP_concentration -= var[2]*len(new_products[var[3]])
                    #count synthesis
                    var[2] += var[2]
                    #something was synthesized in this cycle
                    idle = False
                if idle: break #all primers have been depleted
                #check if: some primers were depleted
                depleted_primers = dict()
                for p_id in primer_concentrations:
                    if  primer_concentrations[p_id] != None \
                    and primer_concentrations[p_id] < 0:
                        consume_ratio = (prev_primers[p_id]/
                                         (prev_primers[p_id] 
                                          - primer_concentrations[p_id]))
                        depleted_primers[p_id] = consume_ratio
                #if so, recalculate consumption for products which use them
                if depleted_primers:
                    for var0, var1 in zip(prev_variants, all_variants):
                        if  not var1[0] in  depleted_primers \
                        and not var1[1] in  depleted_primers:
                            continue
                        if primer_concentrations[var1[0]] is None \
                        or primer_concentrations[var1[1]] is None: 
                            continue
                        consume_ratio0 = depleted_primers[var1[0]] if var1[0] in depleted_primers else 1
                        consume_ratio1 = depleted_primers[var1[1]] if var1[1] in depleted_primers else 1
                        real_addition = (var1[2]-var0[2])*min(consume_ratio0,
                                                              consume_ratio1)
                        #restore subtracted
                        primer_concentrations[var1[0]] += var0[2]
                        primer_concentrations[var1[1]] += var0[2]
                        dNTP_concentration += var0[2]*len(new_products[var1[3]])
                        #subtract real consumption
                        primer_concentrations[var1[0]] -= real_addition
                        primer_concentrations[var1[1]] -= real_addition
                        dNTP_concentration -= real_addition*len(new_products[var1[3]])
                        #count synthesis
                        var1[2] = var0[2] + real_addition
                #at this point concentrations of all primers are >= 0
                #check if: dNTP was used,
                #     and  polymerase activity could handle this cycle.                    
                dNTP_consumption = dNTP_previous - dNTP_concentration
                max_consumptioin = self._polymerase*10e-9*self._elongation_time/30
                #save polymerase shortage information
                if dNTP_consumption > max_consumptioin:
                    if self._reaction_ends[hit_id]['poly']:
                        #if previous cycle passed without a shortage, start new record
                        if cycle - self._reaction_ends[hit_id]['poly'][1] > 1:
                            self._reaction_ends[hit_id]['poly'].append([cycle, cycle])
                        else: #append to the current record 
                            self._reaction_ends[hit_id]['poly'][1] = cycle
                    else: #if there was no shortage, add current new record
                        self._reaction_ends[hit_id]['poly'].append([cycle, cycle])
                #correct consumption if needed
                if dNTP_consumption > max_consumptioin \
                or dNTP_concentration < 0:
                    #ratio is always < 1
                    consume_ratio = min(max_consumptioin/dNTP_consumption,
                                        dNTP_previous/dNTP_consumption)
                    dNTP_concentration    = dNTP_previous
                    primer_concentrations = prev_primers
                    for var0, var1 in zip(prev_variants, all_variants):
                        if primer_concentrations[var1[0]] is None \
                        or primer_concentrations[var1[1]] is None: 
                            continue
                        #due to primer correction var1[2]-var0[2] is not always equal to var0[2]
                        real_addition = (var1[2]-var0[2])*consume_ratio 
                        primer_concentrations[var1[0]] -= real_addition  
                        primer_concentrations[var1[1]] -= real_addition
                        dNTP_concentration -= real_addition*len(new_products[var1[3]])
                        var1[2] = var0[2] + real_addition
                #check again if some primers were depleted
                for p_id in primer_concentrations:
                    if primer_concentrations[p_id] <= 0:
                        primer_concentrations[p_id] = None
                #count the cycle of the reaction and for each product separately
                self._reaction_ends[hit_id]['cycles'] = cycle
                for var0, var1 in zip(prev_variants, all_variants):
                    if var0[2] != var1[2]:
                        self._reaction_ends[hit_id]['products'][var1[3]] = cycle
                #if no dNTP left, end the reaction
                if dNTP_concentration <= 0: break
            #end cycles
            self._reaction_ends[hit_id]['dNTP'] = dNTP_concentration/4.0 if dNTP_concentration >= 0 else 0
            #sum up all variants of each product
            for var in all_variants:
                new_products[var[3]].quantity += var[2]
            #add new products list to new hits list
            new_hits[hit_id] = new_products
        #store new hits if there are any
        if new_hits:
            self._products = new_hits
            self._nonzero  = True
        else: self._products = dict()
    #end def

    
    @classmethod
    def _construct_histogram(cls, main_title, products, reaction_ends=None, with_titles=True):
        text_width = StringTools.text_width
        #construct histogram
        histogram  = dict()
        max_start  = max(len(str(p.start)) for p in products)
        max_end    = max(len(str(p.end))   for p in products)
        max_len    = max(len(str(len(p)))  for p in products)
        max_title  = 0
        if with_titles:
            max_title  = max(len(p.name) for p in products)
        spec_limit = text_width-cls._hist_width-2
        for product in products:
            product_spec   = '%d%s bp [%d%s-%s%d]' \
                           % (len(product),
                              ' '*(max_len-len(str(len(product)))),
                              product.start,
                              ' '*(max_start-len(str(product.start))),
                              ' '*(max_end-len(str(product.end))),
                              product.end)
            if reaction_ends:
                product_spec += ' %d cycles' \
                              % reaction_ends[product.name]['products'][hash(product)]
            if with_titles: 
                title_limit  = spec_limit - len(product_spec) -1
                title_spacer = max_title - len(product.name)
                if title_limit < 0: colname = product_spec
                elif title_limit < len(product.name):
                    colname = product.name
                    colname = (colname)[:title_limit/2-1] + '-' \
                            + (colname)[len(colname)-title_limit/2:] + ' '\
                            + product_spec
                else:
                    colname = (product.name + ' '*title_spacer)[:title_limit] \
                            + ' ' + product_spec
            else:
                colname = product_spec+' '*(spec_limit - len(product_spec))
            histogram[colname] = [product.quantity, len(product)]
        #sort histogram by product length, from the longest to the shortest
        histogram = sorted(zip(histogram, histogram.values()), key=lambda x: x[1][1], reverse=True)
        histogram = tuple((line[0], line[1][0]) for line in histogram)
        #format histogram
        return cls._format_histogram(main_title, histogram)
    #end def
    
    
    @classmethod
    def _format_histogram(cls, title, histogram):
        histogram_string = ''
        #maximum column name width and value
        max_name  = StringTools.text_width-cls._hist_width-2
        max_value = max(c[1] for c in histogram)
        #cut name column title if necessary 
        if len(title) > max_name:
            title = title[:max_name]
        #spacers
        names_spacer    = max_name - len(title)
        coltitle_spacer = cls._hist_width-len(cls._hist_column_title)
        #column_titles 
        histogram_string += '-'*(names_spacer/2)+title        + \
                            '-'*(names_spacer-names_spacer/2) + '|'+ \
                            ' '*(coltitle_spacer/2) + cls._hist_column_title + \
                            ' '*(coltitle_spacer-coltitle_spacer/2)+'|\n'
        #histogram lines
        for col in histogram:
            #line title
            if len(col[0]) > max_name:
                name = col[0][:max_name]
            else: name = col[0]
            name_spacer = max_name - len(name)
            histogram_string += name + ' '*name_spacer
            #line value
            hist_value  = int(round((cls._hist_width*col[1])/max_value))
            col_spacer  = cls._hist_width - hist_value
            #value figure
            value_str   = TD_Functions.format_concentration(col[1])
            #value bar
            if len(value_str) < col_spacer:
                _spacer = col_spacer-len(value_str)
                histogram_string += ':' + '#'*hist_value + ' '*_spacer 
                histogram_string += value_str + ':' + '\n'
            else:
                _bar = hist_value-len(value_str)
                histogram_string += ':' + '#'*(_bar/2) + value_str 
                histogram_string += '#'*(_bar-_bar/2) + ' '*col_spacer + ':' + '\n'
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
        max_mark    = max(len(str(l[2])) for l in phoresis)*2
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
                    header_string += str(shortage_period[0])
                else:
                    header_string += '%d-%d' % shortage_period
            header_string += wrap_text(shortage_string + '\n')
        header_string += '\n\n'
        return header_string

        
    def all_products_histogram(self):
        if not self._nonzero: 
            return '\nNo PCR products found, or quantities have not been calculated.\n'
        all_products = list(chain(*(p.values() for p in self._products.values())))
        return self._construct_histogram(self._all_products_title, all_products)
    #end def
    
    
    def per_hit_histogram(self, hit):
        if not self._nonzero: 
            return '\nNo PCR products found, or quantities have not been calculated.\n'
        products = self._products[hit].values()
        return self._construct_histogram('%s' % hit, products, self._reaction_ends, with_titles=False)
    #end def
    
    
    def all_products_electrophoresis(self):
        if not self._nonzero: 
            return '\nNo PCR products found, or quantities have not been calculated.\n'
        all_products = list(chain(*(p.values() for p in self._products.values())))
        return self._construct_electrophoresis(all_products)
    #end def
    
    
    def per_hit_electrophoresis(self, hit):
        if not self._nonzero: 
            return '\nNo PCR products found, or quantities have not been calculated.\n'
        products = self._products[hit].values()
        return self._construct_electrophoresis(products)
    #end def
    
    
    def all_graphs_grouped_by_hit(self):
        if not self._nonzero: 
            return '\nNo PCR products found, or quantities have not been calculated.\n'
        all_graphs = ''
        for hit in self._products.keys():
            all_graphs += self.per_hit_histogram(hit)
            all_graphs += self.per_hit_header(hit)
            all_graphs += self.per_hit_electrophoresis(hit)
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
        expl_string += hr(' estimation of PCR product quantities ')
        expl_string += 'Value of an objective function at the solution ' + \
                       '(the lower the better):\n   %e\n' % \
                        self._max_objective_value
        expl_string += wrap_text('This value shows "distance" to the solution of '
                           'the system of equilibrium equations which were used '
                           'to calculate quantities of PCR products.\n\n')
        expl_string += '\n'
        return expl_string
    #end def    
#end class



#test
if __name__ == '__main__':
    from SecStructures import Duplex
    from Primer import Primer
    import os
    os.chdir('../')
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC
    TD_Functions.C_DNA  = 1e-10
    TD_Functions.PCR_T  = 60
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
    #print 'products:'
    #print PCR_res._products
    #print bool(PCR_res)
    print PCR_res.format_report_header()
    print PCR_res.format_quantity_explanation()
    print PCR_res.all_products_histogram()
    print ''
    print PCR_res.all_graphs_grouped_by_hit()