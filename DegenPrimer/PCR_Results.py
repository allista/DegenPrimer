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

from math import log, exp
from itertools import chain
from scipy.optimize import fmin_slsqp
from OligoFunctions import dimer
from TD_Functions import dimer_dG, equilibrium_constant
import TD_Functions
import StringTools

class PCR_Results(object):
    '''In silica PCR results standardization 
    and methods for their visualization.'''

    #histogram
    _all_products_title = 'PCR-product'
    _hist_column_title  = '---relative concentration---'
    _hist_width         = len(_hist_column_title)
    #electrophoresis
    _window_percent     = 0.05
    _precision          = 1e6
    #quantity solution accuracy
    _acc                = 1e-10
    #equilibrium constant ratio threshold
    _constant_ratio_threshold = 1e-6

    def __init__(self, 
                 min_amplicon,      #minimum amplicon length 
                 max_amplicon,      #maximum amplicon length
                 no_exonuclease,    #if polymerase does not have 3' exonuclease activity, products of primers with 3'-mismatches will be discarded   
                 quantity_threshold #products with relative quantity (ranged [0,1]) lower than the threshold will be discarded
                 ):
        self._min_amplicon = min_amplicon
        self._max_amplicon = max_amplicon
        self._no_exonuclease = no_exonuclease
        self._quantity_threshold = quantity_threshold
        self._products  = dict()
        self._constants = dict() #dictionary of all equilibrium constants 
        self._primers   = [] #list of all primers
        self._templates = [] #list of all template sequences
        self._nonzero   = False
    #end def
    
    
    def __nonzero__(self):
        return self._nonzero
    #end def
    
    
    def add_product(self, 
                    hit,        #name of the target sequence
                    start,      #start position of the product on the target sequence
                    end,        #end position of the product on the target sequence
                    length,     #length of the product
                    fwd_seq,    #part of the target sequence (5'->3') complementary to the forward primer
                    rev_seq,    #part of the target sequence (5'->3') complementary to the reverse primer
                    fwd_primer, #sequence of the forward primer (5'->3')
                    rev_primer  #sequence of the reverse primer (5'->3')
                    ):
        '''add a new product to the list'''
        #construct product record
        new_product = {'title':hit,
                       'start':start,
                       'end':end,
                       'length':length,
                       'fwd_seq':fwd_seq,
                       'rev_seq':rev_seq,
                       'fwd_primer':fwd_primer, 
                       'rev_primer':rev_primer
                       }
        #calculate quantity, check thresholds
        if not self._prepare_product(new_product): return
        #if new product is to be added, reset nonzero flag
        self._nonzero = False
        #if there was no such hit befor, add it
        if hit not in self._products.keys():
            self._products[hit] = dict()
        #calculate product hash
        new_product_hash = hash((new_product['title'],
                                 new_product['start'],
                                 new_product['end'],
                                 new_product['length']))
        #if no such product from this hit, add it to the list
        if new_product_hash not in self._products[hit].keys():
            #delete sequences
            del new_product['fwd_primer'], \
                new_product['rev_primer'], \
                new_product['fwd_seq'], \
                new_product['rev_seq']
            #add new product
            self._products[hit][new_product_hash] = new_product
        else: #only add sequences
            self._products[hit][new_product_hash]['fwd_primers']   |= new_product['fwd_primers'] 
            self._products[hit][new_product_hash]['rev_primers']   |= new_product['rev_primers']
            self._products[hit][new_product_hash]['fwd_templates'] |= new_product['fwd_templates']
            self._products[hit][new_product_hash]['rev_templates'] |= new_product['rev_templates']
    #end def
    
    
    def _prepare_product(self, product):
        '''
        * filter product using thresholds
        * add sequences to global lists
        * calculate equilibrium constants
        * prepare sets of primers and templates which give this particular product
        '''
        #check thresholds:
        #amplicon length
        if  product['length']   < self._min_amplicon \
        or  product['length']   > self._max_amplicon:
            return False
        #3' mismatches
        fwd_dimer = dimer(product['fwd_primer'], product['fwd_seq'])
        rev_dimer = dimer(product['rev_primer'], product['rev_seq'])
        if self._no_exonuclease:
            if fwd_dimer.fwd_matches()[-1] < len(product['fwd_primer'])-1 \
            or rev_dimer.fwd_matches()[-1] < len(product['rev_primer'])-1:
                return False
        #add primers and templates if nesessary
        if not product['fwd_primer'] in self._primers:
            self._primers.append(product['fwd_primer'])
        if not product['rev_primer'] in self._primers:
            self._primers.append(product['rev_primer'])
        if not product['fwd_seq'] in self._templates:
            self._templates.append(product['fwd_seq'])
        if not product['rev_seq'] in self._templates:
            self._templates.append(product['rev_seq'])
        #calculate equilibrium constants for each pair [primer*template]
        fwd_dG = dimer_dG(fwd_dimer, product['fwd_primer'], product['fwd_seq'])
        rev_dG = dimer_dG(rev_dimer, product['rev_primer'], product['rev_seq'])
        
        fwd_hash = hash(str((product['fwd_primer']), str(product['fwd_seq'])))
        rev_hash = hash(str((product['rev_primer']), str(product['rev_seq'])))
        
        self._constants[fwd_hash] = equilibrium_constant(fwd_dG, TD_Functions.PCR_T)
        self._constants[rev_hash] = equilibrium_constant(rev_dG, TD_Functions.PCR_T)
        #prepare primer-template sets
        product['fwd_primers']   = set(str(product['fwd_primer']))
        product['rev_primers']   = set(str(product['rev_primer']))
        product['fwd_templates'] = set(str(product['fwd_seq']))
        product['rev_templates'] = set(str(product['rev_seq']))
        return True
    #end def
    
    
    def _construct_histogram(self, main_title, products, with_titles=True):
        text_width = StringTools.text_width
        #construct histogram
        histogram  = dict()
        max_start  = max(len(str(p['start']))  for p in products)
        max_end    = max(len(str(p['end']))    for p in products)
        max_len    = max(len(str(p['length'])) for p in products)
        max_title  = 0
        if with_titles:
            max_title  = max(len(p['title']) for p in products)
        spec_limit = text_width-self._hist_width-2
        for product in products:
            product_spec   = '%d%s bp [%d%s-%s%d]' % (product['length'],
                                                      ' '*(max_len-len(str(product['length']))),
                                                      product['start'],
                                                      ' '*(max_start-len(str(product['start']))),
                                                      ' '*(max_end-len(str(product['end']))),
                                                      product['end'])
            if with_titles: 
                title_limit  = spec_limit - len(product_spec) - 1
                title_spacer = max_title - len(product['title'])
                if title_limit < 0: title_limit = 0
                colname = (product['title']+' '*title_spacer)[:title_limit] + ' ' + product_spec
            else:
                colname = product_spec+' '*(spec_limit - len(product_spec))
            histogram[colname] = [product['quantity'],product['length']]
        #sort histogram by product length, from the longest to the shortest
        histogram = sorted(zip(histogram, histogram.values()), key=lambda x: x[1][1], reverse=True)
        histogram = tuple((line[0], line[1][0]) for line in histogram)
        #format histogram
        return self._format_histogram(main_title, histogram)
    #end def
    
    
    def _format_histogram(self, title, histogram):
        text_width = StringTools.text_width
        histogram_string = ''
        #maximum column name width and value
        max_name  = text_width-self._hist_width-2
        max_value = max(c[1] for c in histogram)
        #cut name column title if necessary 
        if len(title) > max_name:
            title = title[:max_name]
        #spacers
        names_spacer = max_name - len(title)
        #column_titles 
        histogram_string += '-'*(names_spacer/2)+title        + \
                            '-'*(names_spacer-names_spacer/2) + \
                            '|'+self._hist_column_title+'|\n'
        #histogram lines
        for col in histogram:
            if len(col[0]) > max_name:
                name = col[0][:max_name]
            else: name = col[0]
            name_spacer = max_name - len(name)
            hist_value  = int(round((self._hist_width*col[1])/max_value))
            col_spacer  = self._hist_width - hist_value
            histogram_string += name + ' '*name_spacer
            histogram_string += ':' + '#'*hist_value + ' '*col_spacer + ':' + '\n'
        histogram_string += '\n'
        return histogram_string
    #end def
    
    
    def _construct_electrophoresis(self, products):
        max_len      = max(p['length'] for p in products)
        window       = int(max_len*self._window_percent)
        nearest_srip = max_len-window
        max_len_log  = int(log(max_len)*self._precision)
        min_len_log  = int(log(min(p['length'] for p in products))*self._precision)
        window_log   = int((log(max_len)-log(nearest_srip))*self._precision)
        #construct phoresis
        phoresis    = [[l,0,int(exp(l/self._precision))] for l in range(min_len_log, max_len_log+window_log, window_log)]
        for product in products:
            l = int(log(product['length'])*self._precision)
            p = (l - min_len_log)/window_log
            phoresis[p][1] += product['quantity']
        phoresis = phoresis[::-1]
        #format electrophorogram
        return self._format_electrophoresis(window, phoresis)
    #end def
    
    
    def _format_electrophoresis(self, window, phoresis):
        text_width  = StringTools.text_width
        max_line    = max(l[1] for l in phoresis)
        max_mark    = max(len(str(l[2])) for l in phoresis)*2
        line_width  = text_width - max_mark - 7 #mark b :###   :
        #format phoresis
        phoresis_text = ''
        phoresis_text += ' '*(max_mark+5)+':'+'-'*line_width+':\n'
        phoresis_text += ' '*(max_mark+5)+':'+' '*line_width+':\n'
        for l in range(len(phoresis)):
            line = phoresis[l]
            prev_mark   = str(phoresis[l-1][2]) if l > 0 else str(line[2]+window)
            mark_spacer = max_mark - len(str(line[2])) - len(prev_mark)
            line_value  = int(round((line_width*line[1])/max_line))
            line_spacer = line_width - line_value 
            phoresis_text += '%s-%d%s bp :%s%s:\n' % (prev_mark,
                                                      line[2], 
                                                      ' '*mark_spacer,
                                                      '#'*line_value,
                                                      ' '*line_spacer) 
        phoresis_text += ' '*(max_mark+5)+':'+' '*line_width+':\n'
        phoresis_text += ' '*(max_mark+5)+':'+'-'*line_width+':\n'
        return phoresis_text
    #end def
    
    
    def _all_products(self):
        return list(chain(*(p.values() for p in self._products.values())))
    
    def hits(self):
        return list(self._products.keys())
    
    
    def compute_quantities(self):
        if not self._products: return
        #initial concentrations
        P    = TD_Functions.C_Prim*1e-6 #M
        D    = TD_Functions.C_DNA *1e-9 #M
        DUP  = min(P,D) #maximum concentration of a duplex
        
        #remove products with small constants
        max_K    = max(chain(*(p['constants'] for p in self._all_products())))
        filtered_hits = dict()
        for hit in self._products:
            filtered_products = dict()
            for product_hash in self._products[hit]:
                product = self._products[hit][product_hash]
                filtered_constants = []
                for c in product['constants']:
                    if c/max_K >= self._constant_ratio_threshold:
                        filtered_constants.append(c)
                if filtered_constants:
                    product['constants'] = filtered_constants
                    filtered_products[product_hash] = product
            if filtered_products:
                filtered_hits[hit] = filtered_products
        self._products = filtered_hits 
        #all constants list            
        K        = tuple(chain(*(p['constants'] for p in self._all_products())))
        print 'K: %d items' % len(K)  #TODO: remove
        print K     #TODO: remove
        
        #factory that generates left side function of the i-th's equation of the system
        def function_factory(i):
            def func(r):
                DUPi = DUP*r[i]
                return (DUPi-((P-sum(DUP*r[j] for j in range(len(r))))*(D-DUPi))*K[i])*1e8/(len(K)**4)
            return func        
        #all left-side functions
        l_funcs  = list(function_factory(i) for i in range(len(K)))
        #objective function
        objective_function =  lambda(r): sum(l_funcs[i](r)**2 for i in range(len(K)))
        #initial estimation of a solution
        
        r0    = tuple(k/max_K/3 for k in K)
        print 'r0:'  #TODO: remove
        print r0     #TODO: remove
        print 'obj_func(r0):'   #TODO: remove
        print objective_function(r0)   #TODO: remove
        #solution bounds
        bounds = ((0,1),)*len(K)
        #try to find the solution
        sol    = fmin_slsqp(objective_function, r0, acc=self._acc*len(K), bounds=bounds, full_output=True)#, disp=0)
        print 'solution' #TODO: remove
        print sol #TODO: remove
        if sol[3] != 0 or sol[1] > self._acc*10*len(K): #TODO: need to handle this situation properly
            raise UserWarning('''No solution for the set of converion degrees 
            of all annealing products has been found''')
        #else, compute quantities from the solution
        print 'max %f; min %f' % (max(sol[0]), min(sol[0])) #TODO: remove 
        sol = sol[0]
        for hit in self._products:
            for product_hash in self._products[hit]:
                product = self._products[hit][product_hash]
                #calculate PCR product quantity as a sum 
                #of conversion degrees of singular annealing products
                for i in range(len(product['constants'])):
                    if i == 0: product['quantity'] = sol.pop()
                    else: product['quantity'] += sol.pop()
                #if a quantity is smaller than a threshold, remove such product 
                if product['quantity'] < self._quantity_threshold:
                    del self._products[hit][product_hash]
                    if not self._products[hit]:
                        del self._products[hit]
        #if there're stil some products left, set nonzero flag to True
        if self._products: self._nonzero = True
    #end def
    
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
        return self._construct_histogram('%s' % hit, products, with_titles=False)
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
            all_graphs += self.per_hit_electrophoresis(hit)
        return all_graphs
#end class



#test
if __name__ == '__main__':
#    from scipy.optimize import fmin_slsqp, fsolve
#
#    k = [1,2,3,4,5]
#    bounds = ((0,1),)*5
#    a = 5
#    b = -3
#
#    def f_factory(i):
#        def func(r):
#            global a,b,k
#            return a*r[i]**2+(b+sum(r))*r[i]-k[i]
#        return func
#    
#    functions = tuple(f_factory(i) for i in range(len(k)))
#    
#    sys_func = lambda(r): \
#                    tuple(functions[i](r) for i in range(len(k)))
#    
#    obj_func = lambda(r): \
#        sum(functions[i](r)**2 for i in range(len(k)))
#    
#    for i in range(5,10):
#        r0 = (i/10.0,)*5
#        acc=1e-10
#        sol1 = fmin_slsqp(obj_func, r0, acc=acc, bounds=bounds, full_output=True)
#        sol2 = fsolve(sys_func, r0).tolist()
#        if sol1[3] == 0 and sol1[1] <= acc:
#            print 'sol1:\n', sol1[0]
#            print sys_func(sol1[0])
#            print obj_func(sol1[0])
#        print 'sol2:\n', sol2
#        print sys_func(sol2)
#        print obj_func(sol2)
    
    
    
    from Bio.Seq import Seq
    PCR_res = PCR_Results(50, 3000, False, 0)
    PCR_res.add_product('Pyrococcus yayanosii CH1', 982243, 983644, 1420, 
                        Seq('ATATTCTACGACGGCTATCC').reverse_complement(), 
                        Seq('CAAGGGTTAGAAGCGGAAG').complement(), 
                        Seq('ATATTCTACAACGGCTATCC'), 
                        Seq('CAAGGGCTAGAAACGGAAG'[::-1]))
    PCR_res.add_product('Pyrococcus yayanosii CH1', 962243, 963644, 1420, 
                        Seq('ATATGCTACGACGGCTATCC').reverse_complement(), 
                        Seq('GAAGGGTTAGACGCGGAAG').complement(), 
                        Seq('ATATTCTACAACGGTTATCC'), 
                        Seq('CAAGAGCTAGAAACGGAAG'[::-1]))
    PCR_res.add_product('Pyrococcus yayanosii CH1', 922243, 923644, 1420, 
                        Seq('ATATTCTACGACAGCTATCC').reverse_complement(), 
                        Seq('CAAGGCTTAGAAGCGGAAG').complement(), 
                        Seq('ATATTCTACAACGTCTATCC'), 
                        Seq('CAAGGGCTAGAGACGGAAG'[::-1]))
    print 'products:'
    print PCR_res._products
    PCR_res.compute_quantities()
    print bool(PCR_res)
    print PCR_res.all_products_histogram()
    print ''
    print PCR_res.all_products_electrophoresis()
    print ''
    print PCR_res.all_graphs_grouped_by_hit()