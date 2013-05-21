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
import TD_Functions
from StringTools import print_exception


class Primer(object):
    '''Representation of a degenerate primer'''

    IUPAC_ambiguous = {'R': ['A', 'G'],
                       'Y': ['C', 'T'],
                       'S': ['G', 'C'],
                       'W': ['A', 'T'],
                       'K': ['G', 'T'],
                       'M': ['A', 'C'],
                       'B': ['C', 'G', 'T'],
                       'D': ['A', 'G', 'T'],
                       'H': ['A', 'C', 'T'],
                       'V': ['A', 'C', 'G'],
                       'N': ['A', 'T', 'C', 'G']}

    
    def __init__(self, seq_rec, concentration):
        self._seq_record    = seq_rec
        self._concentration = concentration
        self._unambiguous   = self.generate_unambiguous(seq_rec)
        #annealing temperatures
        self._Tm_max  = None
        self._Tm_mean = None
        self._Tm_min  = None
        self.calculate_Tms()
    #end def

    def __len__(self): return len(self._seq_record)

    def __contains__(self, seq):
        return str(seq) in self.str_sequences

    def __nonzero__(self):
        return len(self._seq_record.seq) != 0

    def __str__(self):
        return str(self._seq_record.seq)
    
    def __repr__(self):
        rep  = '%(id)s:    %(len)db    5\'-%(seq)s-3\'\n'
        if self._unambiguous:
            rep += '   Tm max  = %(Tm_max).1f\n'
            rep += '   Tm mean = %(Tm_mean).1f\n'
            rep += '   Tm min  = %(Tm_min).1f\n'
            rep += '   Unambiguous primers:   %(unambs)d \n'
            rep += '   Concentration of each: %(cons)s\n'
        else: 
            rep += '   Tm:       %(Tm_max).1f\n'
            rep += '   Concentration: %(cons)s'
        return rep % {'id':      self.id,
                      'seq':     str(self._seq_record.seq),
                      'len':     len(self._seq_record.seq),
                      'unambs':  len(self._unambiguous),
                      'Tm_max':  self._Tm_max,
                      'Tm_mean': self._Tm_mean,
                      'Tm_min':  self._Tm_min,
                      'cons':    TD_Functions.format_concentration(self.concentration)}
    #end def


    def calculate_Tms(self):
        if self.self_complement: return
        temperatures  = []
        concentration = self.concentration 
        for sequence in self.sequences:
            temperatures.append(TD_Functions.primer_template_Tm(sequence, concentration))
        self._Tm_max  = max(temperatures)
        self._Tm_min  = min(temperatures)
        self._Tm_mean = (self._Tm_max + self._Tm_min)/2.0
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
        master_sequence = [set() for i in range(seq_len)]
        for seq_rec in seq_records:
            for i in range(seq_len):
                master_sequence[i].add(seq_rec.seq[i])
        for i in range(seq_len):
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
    
    @property
    def id(self): return self._seq_record.id

    @property
    def master_sequence(self): return self._seq_record

    @property
    def total_concentration(self): return self._concentration
    
    @total_concentration.setter
    def total_concentration(self, C): 
        self._concentration = C
        self.calculate_Tms()
    #end def

    @property
    def all_seq_records(self):
        if self._unambiguous:
            return [self._seq_record,] + self._unambiguous
        else: return [self._seq_record,]
    #end def

    @property
    def seq_records(self):
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
    def concentration(self):
        if self._unambiguous:
            return self._concentration/float(len(self._unambiguous))
        else: return self._concentration
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
        unambiguous_strings = ['']
        #generate list of sequence strings
        for letter in seq_rec.seq:
            n_sequences = len(unambiguous_strings)
            if letter in ['A', 'T', 'G', 'C']:
                for s in range(n_sequences):
                    unambiguous_strings[s] += letter
            elif letter in cls.IUPAC_ambiguous:
                replacements = cls.IUPAC_ambiguous[letter]
                n_replacements = len(replacements) 
                unambiguous_strings *= n_replacements
                for r in range(n_replacements):
                    for s in range(n_sequences):
                        unambiguous_strings[s+n_sequences*r] += replacements[r]
            else: raise ValueError('Unknown letter %s in sequence:\n%s' % (letter, str(seq_rec)))
        #make a list of SeqRecords
        unambiguous_seq_list = list()
        num_unambiguous = len(unambiguous_strings)
        id_template = '%s_%0' + '%d' % int(log(num_unambiguous, 10)+1) + 'd'
        if num_unambiguous < 2: return unambiguous_seq_list
        for s in range(num_unambiguous):
            un_seq = Seq(unambiguous_strings[s], IUPAC.unambiguous_dna)
            unambiguous_seq_list.append(SeqRecord(un_seq, description=seq_rec.description, 
                                                  id=(id_template % (seq_rec.id, s+1))))
        return unambiguous_seq_list
    #end def
#end class


def load_sequence(seq_string, rec_id='', desc=''):
    '''generate a SeqRecord object with sequence from a raw string or a file'''
    #check seq_string to be a sequence
    _is_sequence  = True
    full_alphabet = ['A', 'T', 'C', 'G'] + Primer.IUPAC_ambiguous.keys() 
    for letter in seq_string:
        _is_sequence &= letter in full_alphabet
        if not _is_sequence: break
    #if it is, make a SeqRecord out of it
    if _is_sequence:
        try:
            record = SeqRecord(Seq(seq_string, IUPAC.ambiguous_dna).upper(), 
                               id=rec_id, name=desc, description=desc)
        except Exception, e:
            print_exception(e)
            raise ValueError('Primer.load_sequence: unable to interpret %s as a sequence data.' % seq_string)
        return record
    #else, check if seq_string is an existing file
    elif os.path.isfile(seq_string):
        filename = seq_string
        #check file format by extension
        genbank_pattern = re.compile(".*(\.gb|\.gbk)$")
        fasta_pattern   = re.compile(".*(\.fa|\.faa|\.fasta)$")
        if genbank_pattern.match(filename): filetype = "gb"
        elif fasta_pattern.match(filename): filetype = "fasta"
        else:
            raise ValueError(('Unable to guess format of %s'
                              '*it is expected that GenBank '
                              'files have .gb or .gbk extension '
                              'and FASTA files have .fa, .faa or '
                              '.fasta extension.') % filename)
        #parse file
        #only the first record is loaded
        try:
            record = SeqIO.parse(filename, filetype).next().upper()
        except Exception, e:
            print_exception(e)
            raise ValueError('Primer.load_sequence: unable to parse %s.' % filename)
        #set record id, name and description
        if rec_id: 
            record.id = rec_id
        if not record.name:
            record.name = desc
        if not record.description:
            record.description = desc
        return record
    else:
        raise ValueError('No such file o directory: %s' % seq_string)
#end def