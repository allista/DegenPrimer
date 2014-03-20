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
Created on Mar 15, 2013

@author: Allis Tauri <allista@gmail.com>
'''

from copy import deepcopy
from StringTools import hr
import TD_Functions as tdf


class Region(object):
    '''Region of a sequence.'''
    __slots__ = ['name', 'start', 'end', 'forward']
    
    def __init__(self, name, start, end, forward):
        if start < 1 or end < 1:
            raise ValueError('Region: start and end of a sequence region should be grater than 1.')
        if start > end:
            raise ValueError('Region: start of a sequence region should be less than it\'s end')
        self.name    = name
        self.start   = start
        self.end     = end
        self.forward = forward
    #end def

    def pretty_print(self, with_name=True):
        rep  = ''
        if with_name:
            rep += 'target: %s\n' % self.name
        rep += 'strand: %s\n' % ('forward' if self.forward else 'reverse')
        rep += 'start:  %d\n' % self.start
        rep += 'end:    %d\n' % self.end
        rep += 'length: %d\n' % len(self)
        return rep
    #end def
    
    def __str__(self): return self.pretty_print()
    
    def __len__(self): return self.end - self.start +1
        
    def __hash__(self): return hash((self.name, self.start, self.end, self.forward))
    
    def __eq__(self, other):
        return self.__hash__() == other.__hash__()
    
    def __ne__(self, other):
        return self.__hash__() != other.__hash__()
        
    def __iadd__(self, T):
        if self.name != T.name \
        or self.forward != T.forward: return self
        self.start = min(self.start, T.start)
        self.end   = max(self.end, T.end)
        return self
    #end def
    
    def overlaps(self, T):
        return (self.name == T.name 
                and
                self.forward == T.forward
                and
                (self.start <= T.start <= self.end)
                or
                (self.start <= T.end <= self.end))
    #end def
#end class


class Product(Region):
    '''Representation of a PCR product. A product is defined as a sequence 
    region bounded by start and end positions which correspond to 3' ends of 
    primers (forward and reverse) which produce this product.'''
    __slots__ = ['quantity', 'cycles',
                 '_fwd_ids', '_rev_ids', 
                 'fwd_primers', 'rev_primers', 
                 '_fwd_margin', '_rev_margin',
                 'fwd_template', 'rev_template',
                 ]
    
    def __init__(self, template_name, start, end, 
                 fwd_primers=None, rev_primers=None):
        Region.__init__(self, template_name, start, end, forward=True)
        #quantity and number of cycles
        self.quantity     = 0
        self.cycles       = 0
        #primers and templates
        self._fwd_ids     = set()
        self._rev_ids     = set()
        self.fwd_primers  = set()
        self.rev_primers  = set()
        self._fwd_margin  = self.start-1
        self._rev_margin  = self.end+1
        if self._fwd_margin < 1: self._fwd_margin = 1
        self.fwd_template = Region(template_name, self._fwd_margin, self._fwd_margin, forward=True)
        self.rev_template = Region(template_name, self._rev_margin, self._rev_margin, forward=False)
        for primer in fwd_primers:
            self.add_fwd_primer(primer)
        for primer in rev_primers:
            self.add_rev_primer(primer)
    #end def
    
    
    def pretty_print(self, with_name=True, include_fwd_3_mismatch=True):
        rep  = Region.pretty_print(self, with_name=with_name)
        rep += '\n'
        rep += 'concentration:    %s\n' % tdf.format_concentration(self.quantity)
        rep += 'number of cycles: %d\n' % self.cycles
        rep += '\n'
        rep += hr(' forward annealing site ')
        rep += self.fwd_template.pretty_print(with_name=False)
        rep += '\n'
        rep += hr(' reverse annealing site ')
        rep += self.rev_template.pretty_print(with_name=False)
        rep += '\n'
        rep += hr(' forward primers ')
        fwd_primers = list(self.fwd_primers)
        fwd_primers.sort(key=lambda x: x[1])
        for primer, _id in fwd_primers:
            rep += _id + ':\n' + primer.print_most_stable(include_fwd_3_mismatch) + '\n'
        rep += '\n'
        rep += hr(' reverse primers ')
        rev_primers = list(self.rev_primers)
        rev_primers.sort(key=lambda x: x[1])
        for primer, _id in rev_primers:
            rep += _id + ':\n' + primer.print_most_stable(include_fwd_3_mismatch) + '\n'
        return rep
    #end def
    
    
    def __str__(self):
        rep  = Region.pretty_print(self)
        rep += 'concentration:    %s\n' % tdf.format_concentration(self.quantity)
        rep += 'number of cycles: %d\n' % self.cycles
        rep += '\nforward annealing site:\n'
        rep += str(self.fwd_template)
        rep += '\nreverse annealing site:\n'
        rep += str(self.rev_template)
        rep += '\nforward primers:\n'
        fwd_primers = list(self.fwd_primers)
        fwd_primers.sort(key=lambda x: x[1])
        for primer, _id in fwd_primers:
            rep += _id + ':\n' + repr(primer) + '\n'
        rep += '\nreverse primers:\n'
        rev_primers = list(self.rev_primers)
        rev_primers.sort(key=lambda x: x[1])
        for primer, _id in rev_primers:
            rep += _id + ':\n' + repr(primer) + '\n'
        return rep
    #end def
    
    
    def __iadd__(self, T):
        if self.name != T.name \
        or self.forward != T.forward: return self
        if self.start > T.start:
            self.start        = T.start
            self._fwd_margin  = T._fwd_margin
            self.fwd_template = deepcopy(T.fwd_template)
            self.fwd_primers  = deepcopy(T.fwd_primers)
        elif self.start == T.start:
            for _primer in T.fwd_primers:
                self.add_fwd_primer(_primer)
        if self.end < T.end:
            self.end          = T.end
            self._rev_margin  = T._rev_margin
            self.rev_template = deepcopy(T.rev_template)
            self.rev_primers  = deepcopy(T.rev_primers)
        elif self.end == T.end:
            for _primer in T.rev_primers:
                self.add_rev_primer(_primer)
        return self
    #end def
    
    
    @property
    def fwd_ids(self):
        _ids = [primer[1] for primer in self.fwd_primers]
        _ids.sort() 
        return _ids
    
    @property
    def rev_ids(self): 
        _ids = [primer[1] for primer in self.rev_primers]
        _ids.sort() 
        return _ids

    def add_fwd_primer(self, primer):
        if primer in self.fwd_primers: return
        self.fwd_primers.add(primer)
        self.fwd_template += Region(self.name, 
                                    max(self.start-primer[0].fwd_len, 1), 
                                    max(self.start-1, 1), forward=True)
    #end def
        
    def add_rev_primer(self, primer):
        if primer in self.rev_primers: return
        self.rev_primers.add(primer)
        self.rev_template += Region(self.name, 
                                    self.end+1, 
                                    self.end+primer[0].fwd_len, forward=False)
    #end def
#end class