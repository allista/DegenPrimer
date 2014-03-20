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

from datetime import timedelta
from math import log, exp
from collections import Iterable
from SecStructures import min_K
from SinglePCR import PCR_Base, SinglePCR
from tmpStorage import roDict, cleanup_file, register_tmp_file
from StringTools import wrap_text, line_by_line, hr
import TD_Functions as tdf
import StringTools


class PCR_Simulation_Interface(PCR_Base):
    def __init__(self,
                  abort_event,
                  primers,       #all primers (generally degenerate) that may be present in the system 
                  min_amplicon,  #minimum amplicon length 
                  max_amplicon,  #maximum amplicon length 
                  polymerase,    #polymerase concentration in Units per liter
                  with_exonuclease=False, #if polymerase does have have 3' exonuclease activity, products of primers with 3'-mismatches will also be included
                  num_cycles=20, #number of PCR cycles
                  side_reactions=None,
                  side_concentrations=None,
                  include_side_annealings=False #if True, include annealings which do not give any products into equilibrium system as side reactions
                  ):
        assert isinstance(primers, Iterable), \
        ('PCR_Simulation_Interface: primers should be an iterable. '
         'Given %s instead') % str(primers)
        assert len(primers), 'PCR_Simulation_Interface: no primers given.'

        PCR_Base.__init__(self, abort_event, polymerase, with_exonuclease, num_cycles)
        self._primers                 = primers
        self._min_amplicon            = min_amplicon
        self._max_amplicon            = max_amplicon
        self._include_side_annealings = include_side_annealings
   
        self._elongation_time         = max_amplicon/1000.0 #minutes
        self._primers_concentrations  = dict()
        for primer in self._primers:
            self._primers_concentrations.update(dict.fromkeys(primer.str_sequences, 
                                                              primer.concentration))

        self._side_reactions          = side_reactions if side_reactions else dict()
        self._side_concentrations     = side_concentrations if side_concentrations else dict()
#end class


#manager to isolate forking point in a low-memory process
from PCR_Mixture import ShelvedMixture
from MultiprocessingBase import Parallelizer
from UMP import FuncManager, at_manager

class run_pcr_from_file(object):
    def __init__(self, pcr, counter=None):
        self._pcr     = pcr
        self._counter = counter
    #end def
        
    def __call__(self, mixture_path):
        return self._pcr(ShelvedMixture(mixture_path), self._counter)
#end class

#manager to isolate forking point in a low-memory process
def _ParallelPCR(counter, abort_event, pcr, mixture_paths):
    parallelize = at_manager(Parallelizer, 'parallelize_work_to_shelf')
    return parallelize(abort_event, False, 0.1, run_pcr_from_file(pcr), 
                       mixture_paths, counter=counter)
#end def

PCR_Manager = FuncManager('PCR_Manager', (_ParallelPCR,))


class PCR_Simulation(PCR_Simulation_Interface):
    '''Parallel simulation of multiple PCRs. Visualization of the results.'''
    
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

    def __init__(self, *args, **kwargs): 
        PCR_Simulation_Interface.__init__(self, *args, **kwargs)
        #conditions
        self._reactions_ids           = set()
        self._PCR_mixtures            = dict()
        #results
        self._max_objective_value     = -1 #objective value of the worst solution
        self._reaction_ends           = dict()
        self._products                = dict()
        self._nonzero                 = False  #true if there're some products and their quantities were calculated
        #parallel PCR manager process
        self._pcrM = PCR_Manager()
        self._pcrM.start()
    #end def
    
    def __del__(self):
        try: self._mpm.shutdown()
        except: pass
    #end def
    
    def __nonzero__(self): return self._nonzero
    
    
    def not_empty(self): return bool(self._PCR_mixtures)


    def hits(self): return list(self._reactions_ids)

    
    def add_side_concentrations(self, concentrations):
        if concentrations: self._side_concentrations.update(concentrations)
        
        
    def add_side_reactions(self, reactions):
        if not reactions: return 
        self._side_reactions.update(dict([R for R in reactions.items() 
                                          if R[1].constant >= min_K
                                          or R[1].type == 'A']))
    #end def
    
    
    def add_mixture(self, reaction_id, mixture_path):
        self._reactions_ids.add(reaction_id)
        self._PCR_mixtures[reaction_id] = mixture_path
        register_tmp_file(mixture_path)
    #end def
    
    
    def _new_PCR(self):
        return SinglePCR(self._abort_event, 
                   self._primers_concentrations, 
                   self._elongation_time, 
                   self._polymerase, 
                   self._with_exonuclease, 
                   self._num_cycles, 
                   self._side_reactions, 
                   self._side_concentrations)
    #end def
    

    def run(self, counter):
        if not self._PCR_mixtures: return
        print('PCR Simulation: calculating equilibrium in the system '
              'and running the simulation...')
        pcr = self._new_PCR()
        if len(self._PCR_mixtures) > 1:
            counter.set_subwork(2, (len(self._PCR_mixtures), 0.01))
            results_path = self._pcrM.ParallelPCR(counter[0],
                                                  self._abort_event, pcr,
                                                  self._PCR_mixtures.values())
            results = roDict(results_path)['result']
            cleanup_file(results_path)
            counter[1].done()
        else:
            pcr = run_pcr_from_file(pcr, counter)
            results = (pcr(self._PCR_mixtures.itervalues().next()),)
        if results is None or self.aborted(): return
        for hit_id, obj_val, products, reaction_end in results:
            self._max_objective_value   = max(self._max_objective_value, obj_val)
            self._products[hit_id]      = products
            self._reaction_ends[hit_id] = reaction_end
        #filter out hits without products
        self._products = dict(hit for hit in self._products.items() 
                              if hit[1])
        self._nonzero  = len(self._products) > 0
        print 'PCR Simulation: done.'
    #end def    
    
    
    @classmethod
    def _construct_histogram(cls, main_title, products, reaction_ends=None, with_titles=True, sort=True):
        #construct histogram
        histogram  = []
        max_start  = max(len(str(p.start)) for p in products)
        max_end    = max(len(str(p.end))   for p in products)
        max_len    = max(len(str(len(p)))  for p in products)
        for product in products:
            product_len  = len(product)
            product_spec = ''
            #if the flag is set, print full hit title
            if with_titles:
                product_spec += product.name + '\n'
            #compose product specification
            product_spec  += '%d%s bp [%d%s-%s%d]' \
                           % (product_len,
                              ' '*(max_len-len(str(product_len))),
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
            histogram.append([product_spec, product.quantity, product_len])
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
            value_str  = tdf.format_concentration(col[1])
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
        phoresis    = [[l,0,int(exp(l/self._precision))] for l in xrange(min_len_log, max_len_log+window_log, window_log)]
        for product in products:
            product_len = len(product)
            l = int(log(product_len)*self._precision)
            p = (l - min_len_log)/window_log
            #fluorescence intensity is proportional to the length of a duplex
            phoresis[p][1] += product.quantity*product_len
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
        for l in xrange(len(phoresis)):
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
            % tdf.format_concentration(self._reaction_ends[hit]['dNTP'])
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
        header_string += tdf.format_PCR_conditions(self._primers, self._polymerase)+'\n'
        header_string += hr(' primers and their melting temperatures ')
        for primer in self._primers:
            header_string += repr(primer) + '\n'
        header_string += '\n'
        return header_string
    #end def
    
    
    def format_quantity_explanation(self):
        expl_string  = ''
        expl_string += hr(' estimation of PCR products concentrations ')
        expl_string += 'The value of the objective function of at the solution ' + \
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
            prod_string += hr(' %d products have been found ' % len(self._products[hit]), '*')
            products = self._products[hit].values()
            products.sort(key=lambda x: x.start)
            for pi, product in enumerate(products):
                prod_string += hr(' product %d ' % (pi+1), '=')
                prod_string += product.pretty_print(with_name=False, include_fwd_3_mismatch=self._with_exonuclease)
                prod_string += hr('', '=')
            prod_string += hr('', '*')
        prod_string += '\n'
        return prod_string
#end class