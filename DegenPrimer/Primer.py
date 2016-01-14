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
Created on Jul 11, 2012

@author: Allis Tauri <allista@gmail.com>
'''

import os
import re
from math import log
try:
    from Bio.Alphabet import IUPAC
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
except ImportError:
    print'The BioPython must be installed in your system.'
    raise
from MultiprocessingBase import parallelize_work
from SeqUtils import unambiguous_sequences
import TD_Functions as tdf


class Primer(object):
    '''Representation of a degenerate primer'''

    def __init__(self, seq_rec, concentration, generate_components=False):
        self._seq_record    = seq_rec
        self._concentration = concentration
        self._n_components  = None
        self._unambiguous   = None
        self._initialized   = False 
        #components and annealing temperatures
        self._Tm_max  = None
        self._Tm_mean = None
        self._Tm_min  = None
        if generate_components: 
            self.generate_components()
            self.calculate_Tms()
    #end def
    
    
    def _check_initialized(self):
        assert self._initialized, ('%s primer was not initialized' 
                                   % self._seq_record.id)
        
    def _check_components(self):
        assert self._n_components, ('%s primer was not initialized' 
                                    % self._seq_record.id)


    def __len__(self): return len(self._seq_record)

    def __contains__(self, seq):
        return str(seq) in self.str_sequences

    def __nonzero__(self):
        return len(self._seq_record.seq) != 0

    def __str__(self):
        return str(self._seq_record.seq)
    
    def __repr__(self):
        self._check_initialized()
        rep  = '%(id)s:    %(len)db    5\'-%(seq)s-3\'\n'
        if self._unambiguous:
            rep += '   Tm max  = %(Tm_max).1f\n'
            rep += '   Tm mean = %(Tm_mean).1f\n'
            rep += '   Tm min  = %(Tm_min).1f\n'
            rep += '   Unambiguous primers:   %(unambs)d \n'
            rep += '   Concentration of each: %(cons)s\n'
        else: 
            rep += '   Tm:       %(Tm_max).1f\n'
            rep += '   Concentration: %(cons)s\n'
        return rep % {'id':      self.id,
                      'seq':     str(self._seq_record.seq),
                      'len':     len(self._seq_record.seq),
                      'unambs':  self._n_components,
                      'Tm_max':  self._Tm_max,
                      'Tm_mean': self._Tm_mean,
                      'Tm_min':  self._Tm_min,
                      'cons':    tdf.format_concentration(self.concentration)}
    #end def
    
    def __hash__(self): return hash(str(self._seq_record.seq))

    
    @property
    def id(self): return self._seq_record.id
    
    @property
    def num_components(self): return self._n_components


    def generate_components(self):
        try:
            self._unambiguous  = self.generate_unambiguous(self._seq_record)
            if self._unambiguous: self._n_components = len(self._unambiguous)
            else: self._n_components = 1
            return True
        except ValueError, e:
            print 'Primer: unable to generate unambiguous components:'
            print e
        except Exception, e:
            print e
        return False
    #end def
        
    def calculate_Tms(self, abort_event=None):
        self._check_components()
        if self.self_complement: return False
        #iterative algorithm is faster or no abort event is given to control parallelization
        if self._n_components < 256 or abort_event is None:
            temperatures  = []
            concentration = self.concentration
            for seq in (sequence.seq for sequence in self.seq_records):
                temperatures.append(tdf.primer_template_Tm(seq, concentration))
        else: #parallel algorithm is faster
            temperatures  = parallelize_work(abort_event, 
                                             1, tdf.primer_template_Tm, 
                                             self.sequences, self.concentration)
        if not temperatures: return False
        self._Tm_max  = max(temperatures)
        self._Tm_min  = min(temperatures)
        self._Tm_mean = (self._Tm_max + self._Tm_min)/2.0
        self._initialized = True
        return True
    #end def
    

    @property
    def master_sequence(self): return self._seq_record

    @property
    def all_seq_records(self):
        self._check_components() 
        if self._unambiguous:
            return [self._seq_record,] + self._unambiguous
        else: return [self._seq_record,]
    #end def

    @property
    def seq_records(self):
        self._check_components()
        if self._unambiguous:
            return self._unambiguous
        else: return [self._seq_record,]
    #end def
    
    @property
    def sequences(self):
        return [sequence.seq for sequence in self.seq_records]
    
    @property
    def str_sequences(self):
        return [str(sequence) for sequence in self.sequences]
    
    
    @property
    def total_concentration(self): return self._concentration
    
    @total_concentration.setter
    def total_concentration(self, C): 
        self._concentration = C
        self.calculate_Tms()
    #end def
    
    @property
    def concentration(self):
        self._check_components()
        if self._unambiguous:
            return self._concentration/float(len(self._unambiguous))
        return self._concentration
    #end def


    @property
    def Tm_max(self): return self._Tm_max
    
    @property
    def Tm_mean(self): return self._Tm_mean
    
    @property
    def Tm_min(self): return self._Tm_min


    @property
    def self_complement(self):
        '''return True if sequence is self-complement'''
        return str(self._seq_record.seq) == str(self._seq_record.seq.reverse_complement())
    

    def has_subsequence(self, seq):
        '''Return True if any of unambiguous primers contains seq sequence''' 
        for sequence in self.str_sequences:
            if sequence.count(str(seq)) > 0:
                return True
        return False
    #end def
    
    
    def find_sequences(self, seq):
        '''Return all unambiguous primers which contain seq sequence'''
        sequences = []
        for sequence in self.sequences:
            if str(sequence).count(str(seq)) > 0:
                sequences.append(sequence)
        return sequences
    #end def
    
    
    def find_records(self, seq):
        '''Return all records of unambiguous primers which contain seq sequence'''
        records = []
        for record in self.seq_records:
            if str(record.seq).count(str(seq)) > 0:
                records.append(record)
        return records
    #end def
    
    
    def find_id(self, seq):
        '''Return id of an unubmiguous primer defined by seq'''
        for record in self.seq_records:
            if str(record.seq) == str(seq):
                return record.id
        return ''
    #end def
             

    @classmethod
    def generate_unambiguous(cls, seq_rec):
        '''generate a list of all possible combinations from an ambiguous sequence'''
        unambiguous = unambiguous_sequences(seq_rec.seq)
        unambiguous_seq_list = list()
        num_unambiguous = len(unambiguous)
        id_template = '%s_%0' + '%d' % int(log(num_unambiguous, 10)+1) + 'd'
        if num_unambiguous < 2: return unambiguous_seq_list
        for s in xrange(num_unambiguous):
            un_seq = Seq(unambiguous[s], IUPAC.unambiguous_dna)
            unambiguous_seq_list.append(SeqRecord(un_seq, description=seq_rec.description, 
                                                  id=(id_template % (seq_rec.id, s+1))))
        return unambiguous_seq_list
    #end def
    
    
    @classmethod
    def from_sequences(cls, seq_records, concentration):
        '''Generate new Primer object from a list of unambiguous sequences. 
        If no sequences are provided, or provided sequences deffer in length, 
        return None.'''
        if not seq_records: return None
        #check length equality and alphabeth
        seq_len = len(seq_records[0].seq)
        for seq_rec in seq_records:
            if len(seq_rec.seq) != seq_len:
                return None
            if seq_rec.seq.alphabet != IUPAC.unambiguous_dna:
                return None
        #compose degenerate sequence
        master_sequence = [set() for i in xrange(seq_len)]
        for seq_rec in seq_records:
            for i in xrange(seq_len):
                master_sequence[i].add(seq_rec.seq[i])
        for i in xrange(seq_len):
            if len(master_sequence[i]) == 1:
                master_sequence[i] = master_sequence[i].pop()
                continue
            for letter in cls.IUPAC_ambiguous:
                if master_sequence[i] == set(cls.IUPAC_ambiguous[letter]):
                    master_sequence[i] = letter
                    break
        master_record = SeqRecord(Seq(''.join(l for l in master_sequence), IUPAC.ambiguous_dna),
                                  description=seq_records[0].description,
                                  id=seq_records[0].id)
        return cls(master_record, concentration)
    #end def
#end class


#tests
if __name__ == '__main__':
#    import cProfile
#    primer = Primer(SeqRecord(Seq('ATARTCTYCGAMGGCTATKCAGNCTGGGANGGNTACGNGGGTAAANAAACG'),id='primer1'), 0.9e-6)
    primer = Primer(SeqRecord(Seq('ATARTCTYCGAMGGCNATKCAGGNCTGRGGA'),id='primer1'), 0.9e-6)
    primer.generate_components()
    primer.calculate_Tms()
    print primer.str_sequences
    print primer.concentration
    print repr(primer)
    
#    from timeit import timeit
#    from tests.violin_plot import violin_plot
#    from matplotlib.pyplot import figure, show
#    
#    gen = [timeit('primer.generate_components()', 'from __main__ import primer', number=1) for _i in xrange(1)]
#    gen_mp = [timeit('primer.generate_components_mp()', 'from __main__ import primer', number=1) for _i in xrange(1)]
#    data = [gen, gen_mp]
#    print data
#    
#    fig=figure()
#    ax = fig.add_subplot(111)
#    violin_plot(ax,data,range(len(data)),bp=1)
#    show()
    
    
#    cProfile.run('''
#for _n in xrange(100): primer.generate_components_mp()
#for _n in xrange(100): primer.generate_components()
#''',
#    'gen_components.profile')
    print 'Done.'
    
