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
Created on Jan 14, 2016

@author: Allis Tauri <allista@gmail.com>
'''

import os, re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

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
    
IUPAC_unambiguous = ['A', 'T', 'G', 'C']

def unambiguous_sequences(s):
    '''generate a list of all possible unambiguous DNA sequences from an ambiguous sequence'''
    unambiguous = ['']
    #generate list of sequence strings
    for letter in s:
        n_sequences = len(unambiguous)
        if letter in IUPAC_ambiguous:
            replacements = IUPAC_ambiguous[letter]
            n_replacements = len(replacements) 
            unambiguous *= n_replacements
            for r in xrange(n_replacements):
                for s in xrange(n_sequences):
                    unambiguous[s+n_sequences*r] += replacements[r]
        else:
            for s in xrange(n_sequences):
                unambiguous[s] += letter
    return unambiguous

def load_sequence(seq_string, rec_id='', desc=''):
    '''generate a SeqRecord object with sequence from a raw string or a file'''
    #check if seq_string is an existing file
    if os.path.isfile(seq_string):
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
            print e
            raise ValueError('Primer.load_sequence: unable to parse %s.' % filename)
        #set record id, name and description
        if rec_id: 
            record.id = rec_id
        if not record.name:
            record.name = desc
        if not record.description:
            record.description = desc
        return record
    else: #check if seq_string is a valid sequence
        _is_sequence  = True
        sequence      = seq_string.replace(' ', '').upper()
        full_alphabet = IUPAC_unambiguous + IUPAC_ambiguous.keys() 
        for letter in sequence:
            _is_sequence &= letter in full_alphabet
            if not _is_sequence: break
        #if it is, make a SeqRecord out of it
        if not _is_sequence: 
            raise ValueError(('Primer.load_sequence: given string is not a '
                              'valid sequence nor a valid path '
                              'to a file: %s') % seq_string)
        try:
            record = SeqRecord(Seq(sequence, IUPAC.ambiguous_dna).upper(), 
                               id=rec_id, name=desc, description=desc)
        except Exception, e:
            print e
            raise ValueError('Primer.load_sequence: unable to interpret %s as a sequence data.' % sequence)
        return record
#end def