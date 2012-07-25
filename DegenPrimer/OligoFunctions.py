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
Created on Jul 11, 2012

@author: Allis Tauri <allista@gmail.com>
'''

import os
import re
from StringTools import print_exception
try:
    from Bio.Alphabet import IUPAC
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
except Exception, e:
    print_exception(e)
    raise ImportError('The BioPython must be installed in your system.')


IUPAC_ambiguous={'R': ['A', 'G'],
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


def generate_unambiguous(seq_rec):
    '''generate a list of all possible combinations from an ambiguous sequence'''
    global IUPAC_ambiguous
    unambiguous_strings = ['']
    #generate list of sequence strings
    for letter in seq_rec.seq:
        n_sequences = len(unambiguous_strings)
        if letter in ['A', 'T', 'G', 'C']:
            for s in range(n_sequences):
                unambiguous_strings[s] += letter
        elif letter in IUPAC_ambiguous:
            replacements = IUPAC_ambiguous[letter]
            n_replacements = len(replacements) 
            unambiguous_strings *= n_replacements
            for r in range(n_replacements):
                for s in range(n_sequences):
                    unambiguous_strings[s+n_sequences*r] += replacements[r]
        else: raise ValueError('Unknown letter %s in sequence:\n%s' % (letter, str(seq_rec)))
    #make a list of SeqRecords
    unambiguous_seq_list = list()
    if len(unambiguous_strings) < 2: return unambiguous_seq_list
    for s in range(len(unambiguous_strings)):
        un_seq = Seq(unambiguous_strings[s], IUPAC.unambiguous_dna)
        unambiguous_seq_list.append(SeqRecord(un_seq, description=seq_rec.description, 
                                              id=seq_rec.id+"_"+str(s)))
    return unambiguous_seq_list
#end def


def load_sequence(seq_string, rec_id=None, desc=None):
    '''generate a SeqRecord object with sequence from a raw string or a file'''
    #check if seq_string is a filename
    if os.path.isfile(seq_string):
        filename = seq_string
        #check file format by extension
        genbank_pattern = re.compile(".*(\.gb|\.gbk)$")
        fasta_pattern   = re.compile(".*(\.fa|\.fasta)$")
        if genbank_pattern.match(filename): filetype = "gb"
        elif fasta_pattern.match(filename): filetype = "fasta"
        else:
            print 'Unable to guess format of', filename
            print '*it is expected that GenBank files have .gb or .gbk extension '
            'and FASTA files have .fa or .fasta extension.'
        #parse file
        print 'Parsing', filename
        #only the first record is loaded
        record = SeqIO.parse(filename, filetype).next().upper()
        if rec_id: 
            record.id = rec_id
        if not record.name:
            record.name = desc
        if not record.description:
            record.description = desc
        return record
    #if not, treat it as a raw sequence
    else: return SeqRecord(Seq(seq_string, IUPAC.ambiguous_dna).upper(), 
                           id=rec_id, name=desc, description=desc)
#end def