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
Created on Feb 15, 2014

@author: Allis Tauri <allista@gmail.com>
'''

from BioUtils.Tools.Multiprocessing import MultiprocessingBase, raise_tb_on_error

from .Equilibrium import Equilibrium
from .ReporterInterface import ReporterInterface
from .SecStructures import SecStructures
from BioUtils.Tools.Text import hr, time_hr, wrap_text
from . import TD_Functions

class AllSecStructures(ReporterInterface, MultiprocessingBase):
    '''
    Given two lists of non-degenerate primers calculate all possible secondary 
    structures (hairpins, dimers, cross-dimers) and equilibrium state between them.
    '''
    def __init__(self, abort_event, job_id, primers):
        #initial check
        try:
            if len(primers) == 0: 
                raise ValueError('AllSecStructures: no primers given.')
        except TypeError:
            raise TypeError(('AllSecStructures: primers should be an iterable. '
                            'Given %s instead') % str(primers))
        ReporterInterface.__init__(self)
        MultiprocessingBase.__init__(self, abort_event)
        #all primers list and concentrations dict
        self._primers        = primers
        self._all_primers    = []
        self._concentrations = dict()
        self._reactions      = dict()
        self._self           = []
        self._cross          = []
        self._all_structures = dict()
        self._equilibrium_concentrations = None
        #reports
        self._job_id = job_id
        self._short_structs_filename = job_id+'-structures-short.txt'
        self._full_structs_filename  = job_id+'-structures-full.txt'
    #end def
    
    
    #property functions
    def reactions(self): return self._reactions
    def concentrations(self): return self._concentrations
    def equilibrium_concentrations(self): return self._equilibrium_concentrations
    
    
    @staticmethod
    @MultiprocessingBase.data_mapper
    def _self_structures(seq):
        structures = SecStructures(seq)
        reactions  = structures.compose_reactions() 
        return structures,reactions
    #end def
    
    @staticmethod
    @MultiprocessingBase.data_mapper
    def _cross_structures(index, sequences): 
        structures = SecStructures(sequences[index[0]], sequences[index[1]])
        reactions  = structures.compose_reactions() 
        return structures,reactions
    #end def
    
    @MultiprocessingBase.results_assembler_method
    def _structure_assembler(self, index, result, structures):
        structures[index] = result[0]
        self._reactions.update(result[1])
        self._all_structures.update(dict().fromkeys(result[1].keys(), result[0]))
    #end def
    
    
    def find_structures(self):
        #self structures jobs
        self_works  = []
        for primer in self._primers:
            #list of all primers
            self._all_primers.extend(primer.seq_records)
            #self structure work
            work = self.Work()
            work.start_work(self._self_structures, primer.seq_records, None)
            self_works.append(work)
            self._self.append([None]*len(primer.seq_records))
            work.assemble(self._structure_assembler, self._self[-1])
        #cross structures jobs
        num_primers = len(self._all_primers)
        pair_index  = []
        for i in xrange(num_primers):
            for j in xrange(i+1, num_primers):
                pair_index.append((i,j))
        cross_work = self.Work()
        cross_work.start_work(self._cross_structures, 
                              pair_index, None, self._all_primers)
        self._cross = [None]*len(pair_index)
        cross_work.assemble(self._structure_assembler, self._cross)
        #wait for all works to be finisheds
        if not self.wait(cross_work, *self_works): return False
        #prepare primer concentrations dictionary
        for primer in self._primers:
            self._concentrations.update(dict().fromkeys((p for p in primer.str_sequences), 
                                                        primer.concentration))
        self._have_results = True
        return True
    #end def
    
    @MultiprocessingBase.raise_tb_on_error
    def calculate_equilibrium(self, counter):
        if not self._reactions: return False
        counter.set_subwork(3)
        #calculate equilibrium
        equilibrium = Equilibrium(self._abort_event, self._reactions, self._concentrations)
        equilibrium.calculate(counter[0])
        if equilibrium.solution is None: return False
        #calculate equilibrium primer concentrations
        self._equilibrium_concentrations = dict()
        counter[1].set_work(len(self._concentrations))
        for primer in self._concentrations:
            _C = equilibrium.reactants_consumption(primer)[0]
            self._equilibrium_concentrations[primer] = self._concentrations[primer] - _C
            counter[1].count()
        #set conversion degrees
        counter[2].set_work(len(equilibrium.solution))
        for r_hash in equilibrium.solution:
            self._all_structures[r_hash].set_conversion_degree(r_hash, equilibrium.solution[r_hash])
            counter[2].count()
    #end def
    
    
    def _print_structures(self, structures, full):
        def worker(struct):
            if not struct: return ''
            if full: return struct.formatFull()
            else: return struct.formatShort()
        strings = self.parallelize_work(1, worker, structures)
        return ''.join(strings)
    #end def
    
    
    def print_structures(self, full=True):
        strings = []
        for i, primer in enumerate(self._primers):
            structs = self._print_structures(self._self[i], full)
            if structs:
                strings.append(hr(' secondary structures of %s primers ' % primer.id, '#'))
                strings.append(structs) 
        if self._cross:
            cross_structs = self._print_structures(self._cross, full)
            if cross_structs:
                strings.append(hr(' cross-dimers ', '#'))
                strings.append(cross_structs)
        if strings:
            strings.insert(0, hr(' stable secondary structures ', symbol='#'))
        else: strings.append(hr(' no stable secondary structures found ', symbol='#'))
        return ''.join(strings)
    #end def
    
    def print_structures_short(self):
        return self.print_structures(full=False)
    
    
    def _format_primers_report_header(self):
        header_string  = ''
        header_string += time_hr()
        header_string += wrap_text('For each degenerate primer provided, a set '
                                   'of unambiguous primers is generated. '
                                   'For each such set the minimum, maximum and '
                                   'mean melting temperatures are calculated. '
                                   'For each primer in each set stable self-'
                                   'dimers and hairpins are predicted. '
                                   'For every possible combination of two '
                                   'unambiguous primers cross-dimers are also '
                                   'predicted. If an unambiguous primer is '
                                   'provided, it is treated as a set with a '
                                   'single element.\n\n')
        header_string += hr(' PCR conditions ')
        header_string += TD_Functions.format_PCR_conditions(self._primers)+'\n'
        header_string += hr(' primers and their melting temperatures ')
        for primer in self._primers:
            header_string += repr(primer) + '\n'
        #warning
        if len(self._primers) > 1:
            if abs(self._primers[0].Tm_min - self._primers[1].Tm_min) >= 5:
                header_string += '\nWarning: lowest melting temperatures of sense and antisense primes \n'
                header_string += '         differ more then by 5C\n'
        header_string += '\n'
        return header_string
    #end def
    
    
    def write_report(self):
        if not self._have_results: return
        full_file  = self._open_report('Full list of secondary structures', self._full_structs_filename)
        short_file = self._open_report('Secondary structures', self._short_structs_filename)
        #write header
        structs_header = self._format_primers_report_header() 
        full_file.write(structs_header)
        short_file.write(structs_header)
        #write secondary structures information
        full_file.write(self.print_structures())
        short_file.write(self.print_structures_short())
        #close files
        full_file.close()
        short_file.close()
        print '\nFull report with all secondary structures was written to:\n   %s' % self._full_structs_filename
        print '\nShort report with a summary of secondary structures was written to:\n   %s' % self._short_structs_filename
        self._add_report('Secondary structures', self._short_structs_filename)
    #end def
    
    @raise_tb_on_error
    def write_reports(self): self.write_report()
#end class