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
Created on Jun 23, 2012

@author: Allis Tauri <allista@gmail.com>
'''

import numpy as np
from array import array

from BioUtils.Tools.JSON import JSONattrs
from BioUtils.Tools.Text import hr

from .Equilibrium import Reaction
from . import TD_Functions as tdf
#default filtration parameters
max_dimer_dG   = -3  #kcal/mol #corresponds to equilibrium constant ~100  and conversion degree ~1e-3
max_hairpin_dG =  3  #kcal/mol #corresponds to equilibrium constant ~0.01 and conversion degree ~1e-2
min_K          = 100 #minimum equilibrium constant: annealing reactions with EC less than this will not be taken into account

rc_map = {'A':'T',
          'T':'A',
          'G':'C',
          'C':'G',
          'R':'Y',
          'Y':'R',
          'S':'S',
          'W':'W',
          'K':'M',
          'M':'K',
          'B':'V',
          'V':'B',
          'D':'H',
          'H':'D',
          'N':'N'}

def reverse_complement(seq): #this is much faster than the Bio.Seq.reverse_complement() method
    '''Make a reverse complement of a DNA sequence'''
    rc = ''
    for l in seq.upper(): rc += rc_map[l]
    return rc[::-1]
#end def


class Dimer(JSONattrs):
    '''
    Structure of a dimer without bulges (i.e. asymmetrical loops):
    forward matches is a list of nucleotide indices of forward sequence 
    directed 5'->3' that are paired to reverse strand.
    The offset is a forward-strand index of the reverse strand 3'-end. 
    For example:
                          234  7 :fwd_matches 
                        ATGCACGA :offset=2
                          |||  |
                          CGTTATGT
                          012  5 :rev_matches[n] = fwd_matches[n]-offset
    '''
    __slots__ = ['fwd_matches', 'offset', 'num_matches',
                 'conversion_degree', 'dG', 'K',
                 'fwd_mismatch']
    
    def __init__(self, fwd_matches=None, offset=0): 
        #dG, K and CD are properties of a Duplex, not a Dimer, but since 
        #a Duplex may have several associated dimers, it's convenient to store them here 
        if fwd_matches:
            self.fwd_matches  = tuple(sorted(fwd_matches))
            self.offset       = offset
            self.num_matches  = len(self.fwd_matches)
        else:
            self.fwd_matches  = tuple()
            self.offset       = None
            self.num_matches  = 0    
        #thermodynamic parameters
        self.conversion_degree = None
        self.dG = None
        self.K  = None
        #topological parameters
        self.fwd_mismatch = None
    #end def
    
    def __hash__(self):
        return hash((self.fwd_matches, self.offset))
    
    def __eq__(self, other):
        return self.__hash__() == other.__hash__()
    
    def __ne__(self, other):
        return self.__hash__() != other.__hash__()
    
    def __nonzero__(self): return self.num_matches > 0
    
    @property
    def rev_matches(self):
        return np.array(self.fwd_matches)-self.offset
    
    @property
    def fwd_min(self):
        return self.fwd_matches[0]
    
    @property
    def rev_min(self):
        return self.fwd_matches[0]-self.offset
    
    @property
    def fwd_max(self):
        return self.fwd_matches[-1]
    
    @property
    def rev_max(self):
        return self.fwd_matches[-1]-self.offset
    
    
    @property
    def alternatives(self):
        """remove terminal internal loops from a dimer one by one in all 
        possible combinations to produce a bunch of new structures"""
        #scan for loops from both ends
        if self.num_matches < 2: return []
        margins = [set()]
        for i in xrange(self.num_matches-1):
            margins[-1].add(self.fwd_matches[i])
            if self.fwd_matches[i+1]-self.fwd_matches[i] > 1:
                margins.append(set())
        margins[-1].add(self.fwd_matches[-1])
        num_loops = len(margins)-1
        #of no loops have been found, return empty list
        if not num_loops: return []
        #remove boundary loops one by one
        dimers  = []
        l_dimer = self
        for l in xrange(num_loops):
            r_dimer = l_dimer
            for r in xrange(num_loops,l,-1):
                r_dimer = Dimer(set(r_dimer.fwd_matches)-margins[r], self.offset)
                dimers.append(r_dimer)
            l_dimer = Dimer(set(l_dimer.fwd_matches)-margins[l], self.offset)
            dimers.append(l_dimer)
        return dimers
    #end def
    
    def __enter__(self):
        self.fwd_matches = list(self.fwd_matches)
        return self
    
    def __exit__(self, *ext_info):
        self.fwd_matches = tuple(sorted(self.fwd_matches))
    
    def add(self, fwd, rev):
        if self.offset is None:
            self.offset = fwd-rev
        elif self.offset != fwd-rev:
            raise ValueError('Dimer.add: fwd and rev indexes should be at the '
                             'same distance as all other match pairs. '
                             'Asymmetrical loops are not allowed.')
        self.fwd_matches.append(fwd)
        self.num_matches += 1
    #end def
    
        
    @classmethod
    def from_sequences(cls, forward, revcomp, offset=0):
        '''Form a dimer structure from two sequences using a pair of indexes 
        as an anchor; if anchor is not provided, sequences are left-aligned.
        The forward sequence should be 5'->3' oriented.
        The revcom sequence should be 3'->5' oriented'''
        dimer = cls()
        dimer.offset = offset
        dimer.fwd_matches = []
        if offset > 0: _range = xrange(offset, min(len(forward), len(revcomp)+offset))
        else: _range = xrange(0, min(len(forward)-offset, len(revcomp)+offset))
        for i in _range:
            if forward[i] == revcomp[i-offset]:
                #add indexes directly is much faster than to call Dimer.add method
                dimer.fwd_matches.append(i)
                dimer.num_matches += 1
        dimer.fwd_matches = tuple(dimer.fwd_matches)
        return dimer
    #end def
#end class

class Hairpin(JSONattrs):
    '''
    Structure of a hairpin:
    sequence directed 5'->3' 
    forward matches is a list of near-5'-end nucleotide indices
    reverse matches is a *reversed* list of near-3'-end nucleotide indices
    nucleotide fwd_match[n] canonically matches nucleotide rev_matches[n]
    e.g. fwd_matches = 5'-[2, 4, 7,  9]...-3'
         rev_matches = 3'-[22,20,17,15]...-5'
    '''
    
    __slots__ = ['fwd_matches', 'rev_matches',
                 'conversion_degree', 'dG', 'K']
    
    def __init__(self, fwd_matches, rev_matches): 
        if len(fwd_matches) != len(rev_matches):
            raise ValueError('SecStructures.Dimer: number of forward matches should equal that of reverse matches.')
        self.fwd_matches = tuple(sorted(fwd_matches))
        self.rev_matches = tuple(sorted(rev_matches, reverse=True))
        #thermodynamic parameters
        self.conversion_degree = None
        self.dG = None
        self.K  = None
    
    @classmethod
    def _default(cls): return cls((), ())
#end class


class Duplex(JSONattrs):
    '''Representation of a particular DNA:DNA duplex'''
    __slots__ = ['_fwd_sequence', '_rev_sequence',
                 'fwd_len', 'rev_len', 'name',
                 '_dimers',
                 '_nonzero',
                 '_fwd_matches', '_fwd_mismatches', 'fwd_3_overhang']
    
    def __init__(self, fwd_seq, rev_seq, name='', dimer=None, revcomp=False):
        '''
        fwd_seq -- string, 5'->3' sequence of forward strand
        rev_seq -- string, 5'->3' or 3'->5' sequence of reverse strand
        dimer, optinal -- an initial Dimer structure
        revcom  -- indicates rev_seq orientation: False means 5'->3' 
        '''
        self.name            = name
        self.fwd_len         = len(fwd_seq)
        self.rev_len         = len(rev_seq)
        self._fwd_sequence   = fwd_seq
        self._rev_sequence   = reverse_complement(rev_seq) if revcomp else rev_seq 
        if not dimer: dimer  = Dimer.from_sequences(fwd_seq, 
                                                    rev_seq if revcomp 
                                                    else reverse_complement(rev_seq))
        self._dimers         = [dimer]+dimer.alternatives
        self._fwd_matches    = []
        self._fwd_mismatches = []
        self.fwd_3_overhang  = None
        self._nonzero        = False 
        self._recalculate()
    #end def
    
    @classmethod
    def _default(cls): return cls('', '')
        
    def __nonzero__(self): return self._nonzero
    
    def __hash__(self):
        return hash((self._fwd_sequence, self._rev_sequence, self._dimers[0]))
    
    def __eq__(self, other):
        return self.__hash__() == other.__hash__()
    
    def __ne__(self, other):
        return self.__hash__() != other.__hash__()
    
    def __str__(self):
        duplex_str = '%d variants of annealing sorted by stability:\n' % len(self._dimers)
        for i, dimer in enumerate(self._dimers):
            duplex_str += '\n%d:\n' % (i+1)
            duplex_str += self.print_dimer(self._fwd_sequence, self._rev_sequence, dimer)
            duplex_str += '\n'
        return duplex_str+'\n'
    #end def
    
    def __repr__(self): return str(self)
    
    def print_most_stable(self, include_fwd_3_mismatch=True):
        '''Print the most stable duplex variant.
        If include_fwd_3_mismatch is False, the most stable variant *without*
        3' mismatches on forward strand will be printed'''
        duplex_str = ''
        if not self._dimers: return duplex_str
        if not include_fwd_3_mismatch:
            if self._fwd_matches:
                return self.print_dimer(self._fwd_sequence, self._rev_sequence, 
                                        self._dimers[self._fwd_matches[0]])
            duplex_str += 'There\'s no annealing variant without 3\' mismatches\n'
        return duplex_str+self.print_dimer(self._fwd_sequence, self._rev_sequence, 
                                           self._dimers[0])
    #end def
    
    def _recalculate(self):
        filtered_dimers = []
        for i, dimer in enumerate(self._dimers):
            if not dimer: continue
            #thermodynamical parameters
            dimer.dG = tdf.dimer_dG_corrected(dimer, 
                                              self._fwd_sequence,
                                              self._rev_sequence)
            if dimer.dG > max_dimer_dG: continue
            dimer.K  = tdf.equilibrium_constant(dimer.dG, tdf.PCR_P.PCR_T)
            if dimer.K < min_K: continue
            #mismatch at 3' end of forward sequence
            dimer.fwd_mismatch = dimer.fwd_max < self.fwd_len-1
            filtered_dimers.append(dimer)
        #check if some dimers left
        self._nonzero = bool(filtered_dimers)
        if not self._nonzero:
            self._dimers = tuple() 
            return
        #sort dimers by stability
        self._dimers  = tuple(sorted(filtered_dimers, key=lambda x: x.dG))
        #overhang of 3' forward sequence (no matter what dimer, it's always the same)
        fwd_tail = self.fwd_len-1 - self._dimers[0].fwd_max
        rev_tail = self.rev_len-1 - self._dimers[0].rev_max
        self.fwd_3_overhang = fwd_tail - rev_tail if fwd_tail > rev_tail else 0 
        #sort dimers by 3' mismatch
        for i, dimer in enumerate(self._dimers):
            if dimer.fwd_max < self.fwd_len-1:
                self._fwd_mismatches.append(i)
            else: self._fwd_matches.append(i)
    #end def
    
    @property
    def fwd_seq(self): return self._fwd_sequence
    
    @property
    def rev_seq(self): return self._rev_sequence
    
    @property
    def dimer(self): return self._dimers[0]
    
    @property
    def dimers(self): return self._dimers
    
    @property
    def dG(self): return self._dimers[0].dG
    
    @property
    def K(self): return self._dimers[0].K
    
    @property
    def fwd_3_matches(self): 
        return [self._dimers[i] for i in self._fwd_matches]
    
    @property
    def fwd_3_mismatches(self): 
        return [self._dimers[i] for i in self._fwd_mismatches]
    
    @property
    def have_3_matches(self): return bool(self._fwd_matches)
    
    @classmethod
    def print_dimer(cls, fwd_sequence, rev_sequence, dimer):
        duplex_string = ''
        if not dimer:
            duplex_string += fwd_sequence + '\n'
            duplex_string += rev_sequence[::-1] + '\n'
            return duplex_string
        fwd_matches  = dimer.fwd_matches
        #construct matches string
        matches = array('c', ' '*len(fwd_sequence))
        for m in fwd_matches: matches[m] = '|'
        #calculate spacers
        if dimer.offset > 0:
            duplex_string += "5' "+fwd_sequence+" 3'\n   "
            duplex_string += matches.tostring()+'\n'
            duplex_string += ' '*dimer.offset+"3' "+rev_sequence[::-1]+" 5'\n"
        else:
            offset_str = ''
            if dimer.offset < 0:
                offset_str = ' '*abs(dimer.offset)
                duplex_string += offset_str 
            duplex_string += "5' "+fwd_sequence+" 3'\n   "
            duplex_string += offset_str+matches.tostring()+'\n'
            duplex_string += "3' "+rev_sequence[::-1]+" 5'\n"
        #dG
        if dimer.dG is not None:
            duplex_string += 'dG(%.1fC) = %.2f kcal/mol\n' % (tdf.PCR_P.PCR_T, dimer.dG)
        #conversion degree
        if dimer.conversion_degree is not None:
            duplex_string += 'conversion degree = %.4f%%\n' % (dimer.conversion_degree*100)
        return duplex_string
    #end def
#end class


class SecStructures(object):
    '''
    Given one or two sequence objects compatible with BioPython SeqRecord objects, 
    computes sets of self-dimers, hairpins and (in case of two sequences) cross-dimers.
    To get formatted output just convert to a str: 
    >>> print str(MySecStructs)
    '''
    
    def __init__(self, seq_rec1, seq_rec2=None):
        '''
        seq_rec1 -- Bio.SeqRecord object containing first sequence
        seq_rec2 -- Bio.SeqRecord object containing second sequence. May be None. 
        If so, hairpins and self-dimers for seq_req1 will be calculated.
        Otherwise, only corss-dimers between seq_req1 and seq_req2 will be calculated.
        '''
        #attributes
        #sequences
        self._seq_rec1 = seq_rec1
        self._seq1 = seq_rec1.seq
        #self or cross
        if seq_rec2 and str(seq_rec1.seq) != str(seq_rec2.seq):
            self._seq_rec2 = seq_rec2
            self._seq2 = seq_rec2.seq
        else: 
            self._seq_rec2 = None
            self._seq2 = None
        #structures
        self._seq1_dimers    = None
        self._seq1_hairpins  = None
        self._cross_dimers   = None
        self._3prim_dimers   = 0
        self._3prim_hairpins = 0
        #minimum energy structures
        self._min_3prim_dimer_dG   = None
        self._min_3prim_dimer      = None
        self._min_3prim_hairpin_dG = None
        self._min_3prim_hairpin    = None
        #find secondary structures
        self._recalculate()
    #end def
    
    
    def __nonzero__(self):
        if self._seq2 and self._cross_dimers: return True
        elif self._seq1_dimers or self._seq1_hairpins: return True
        else: return False
    #end def
        
        
    def compose_reactions(self):
        reactions = dict()
        #if selfdimers and hairpins
        if not self._seq_rec2:
            if self._seq1_dimers:
                for d_hash in self._seq1_dimers:
                    D = self._seq1_dimers[d_hash]
                    reactions[d_hash] = Reaction(tdf.equilibrium_constant(D.dG, tdf.PCR_P.PCR_T), 
                                                 str(self._seq1), 
                                                 None, 
                                                 '2A')
            if self._seq1_hairpins:
                for h_hash in self._seq1_hairpins:
                    H = self._seq1_hairpins[h_hash]
                    reactions[h_hash] = Reaction(tdf.equilibrium_constant(H.dG, tdf.PCR_P.PCR_T), 
                                                 str(self._seq1), 
                                                 None, 
                                                 'A')
        else: #crossdimers
            if self._cross_dimers:
                for d_hash in self._cross_dimers:
                    D = self._cross_dimers[d_hash]
                    reactions[d_hash] = Reaction(tdf.equilibrium_constant(D.dG, tdf.PCR_P.PCR_T), 
                                                 str(self._seq1), 
                                                 str(self._seq2), 
                                                 'AB')
        return reactions
    #end def
    
    
    def set_conversion_degree(self, s_hash, r):
        if not self._seq_rec2:
            try: self._seq1_dimers[s_hash].conversion_degree   = r
            except: pass
            else: return True
            try: self._seq1_hairpins[s_hash].conversion_degree = r
            except: pass
            else: return True
        else:
            if s_hash in self._cross_dimers:
                self._cross_dimers[s_hash].conversion_degree   = r
                return True
            else: return False
    #end def
    
    
    def max_conversion_degree(self):
        if self._seq_rec2 and self._cross_dimers:
            return max(self._cross_dimers.values(), key=lambda(s): s.conversion_degree).conversion_degree
        else:
            CG = 0
            if self._seq1_dimers:
                CG += max(self._seq1_dimers.values(), key=lambda(s): s.conversion_degree).conversion_degree
            if self._seq1_hairpins:
                CG += max(self._seq1_hairpins.values(), key=lambda(s): s.conversion_degree).conversion_degree
            return CG if CG else None
    #end def
    
    
    def dimerMin_dG(self):
        if self._seq_rec2 and self._cross_dimers:
            return min(self._cross_dimers.values(), key=lambda(s): s.dG).dG
        elif self._seq1_dimers:
            return min(self._seq1_dimers.values(), key=lambda(s): s.dG).dG
        else: return None
    #end def
    
    
    def dimerMin_3prim_dG(self):
        if self._min_3prim_dimer_dG:
            return self._min_3prim_dimer_dG
        else: return None
    #end def
    
    
    def hairpinMin_dG(self):
        if self._seq_rec2 or not self._seq1_hairpins: 
            return None
        else: return min(self._seq1_hairpins.values(), key=lambda(s): s.dG).dG
    #end def
    
    
    def hairpinMin_3prim_dG(self):
        if self._min_3prim_hairpin_dG:
            return self._min_3prim_hairpin_dG
        else: return None
    #end def
    
    
    def formatFull(self):
        return self._format_str(full=True)
    
    def formatShort(self):
        return self._format_str(full=False)
    
    def __str__(self):
        return self._format_str()
    
    
    def _format_str(self, full=True):
        """If 'full', for single sequence output self-dimers and hairpins.
        For two sequences output only cross-dimers
        If not 'full', for single sequence output COUNTS and min dG of self-dimers and hairpins.
        For two sequences output the same only for cross-dimers"""
        if full:
            format_dimers = self._format_dimers
            format_hairpins = self._format_hairpins
        else:
            format_dimers = self._format_dimers_short
            format_hairpins = self._format_hairpins_short
        string  = '\n'
        if self._seq2:
            string += hr(' %s vs %s cross-dimers ' % (self._seq_rec1.id, self._seq_rec2.id), symbol='=')
            dimers  = self._cross_dimers.values()
            dimers.sort(key=lambda(s): s.dG)
            string += format_dimers(dimers, self._seq1, self._seq2)
        else:
            string += hr(' %s: %s ' % (self._seq_rec1.id, str(self._seq1)), symbol='=')
            string += hr(' self-dimers ')
            dimers  = self._seq1_dimers.values()
            dimers.sort(key=lambda(s): s.dG)
            string += format_dimers(dimers, self._seq1, self._seq1)
            string += hr(' hairpins ')
            hairpins = self._seq1_hairpins.values()
            hairpins.sort(key=lambda(s): s.dG)
            string += format_hairpins(hairpins, self._seq1)
        string += hr('',symbol='=')
        return string
    #end def
    
    
    def _recalculate(self):
        """For single sequence calculate self-dimers and hairpins.
        For two sequences calculate only cross-dimers"""
        self._3prim_dimers         = 0
        self._min_3prim_dimer_dG   = None
        self._min_3prim_dimer      = None
        self._3prim_hairpins       = 0
        self._min_3prim_hairpin_dG = None
        self._min_3prim_hairpin    = None 
        if self._seq2:
            self._seq1_dimers    = None
            self._seq1_hairpins  = None
            self._cross_dimers   = self._find_dimers(self._seq1, self._seq2)['dimers']
        else:
            all_structures       = self._find_dimers(self._seq1, self._seq1)
            self._seq1_dimers    = all_structures['dimers']
            self._seq1_hairpins  = self._find_hairpins(all_structures['hairpins'], self._seq1)
            self._cross_dimers   = None
    #end def
    
    
    def __contains__(self, s_hash):
        if not self._seq_rec2:
            return (s_hash in self._seq1_dimers or s_hash in self._seq1_hairpins) 
        else:
            return (s_hash in self._cross_dimers)
    #end def
    
    
    def _structure_hash(self, struct):
        if not self._seq_rec2:
            return hash((str(self._seq1), 
                         tuple(struct.fwd_matches+struct.rev_matches)))
        else:
            return hash((str(self._seq1), str(self._seq2), 
                         tuple(struct.fwd_matches+struct.rev_matches)))
    #end def
    
    
    def _find_dimers(self, seq1, seq2):
        """Search for all dimers"""
        all_structures = {'dimers': dict(), 'hairpins': dict()}
        self_dimers    = str(seq1) == str(seq2)
        dimers   = all_structures['dimers']
        hairpins = all_structures['hairpins']
        fwd_str  = str(seq1)
        rev_str  = reverse_complement(str(seq2))
        fwd_len  = len(fwd_str)
        rev_len  = len(rev_str)
        full_len = fwd_len+rev_len-1
        for front in xrange(full_len):
            #margins
            fwd_left  = fwd_len - front - 1
            if fwd_left < 0: fwd_left = 0
            #
            fwd_right = fwd_len
            if front >= rev_len: 
                fwd_right -= front - rev_len - 1
            #
            rev_left  = 0
            if front >= fwd_len:
                rev_left  += front - fwd_len + 1
            #   
            rev_right = front + 1
            if front >= rev_len:
                rev_right = rev_len
            #scan for matches
            dimer = Dimer()
            with dimer:
                for fwd,rev in zip(xrange(fwd_left,fwd_right),xrange(rev_left,rev_right)):
                    if(fwd_str[fwd] == rev_str[rev]):
                        dimer.add(fwd, rev)
            #if found some
            if dimer:
                #add a copy to the hairpins dict to search for hairpin structures later
                dimer_hash = self._structure_hash(dimer)
                if self_dimers and dimer_hash not in hairpins: 
                    hairpins[dimer_hash] = dimer
                #find the most stable dimer
                min_dG = max_dimer_dG
                min_dG_dimer = None
                new_dimers   = [dimer,] + dimer.alternatives
                for dimer in new_dimers:
                    #calculate dG for every dimer
                    dimer.dG = tdf.dimer_dG_corrected(dimer, seq1, seq2)
                    #find the most stable
                    if min_dG > dimer.dG:
                        min_dG = dimer.dG
                        min_dG_dimer = dimer
                #if there is a stable dimer, append it
                if min_dG_dimer:
                    dimer_hash = self._structure_hash(min_dG_dimer)
                    if dimer_hash not in dimers:
                        dimers[dimer_hash] = min_dG_dimer
                        #check for 3' dimer
                        if min_dG_dimer.fwd_matches[-1] > fwd_len-4 \
                        or min_dG_dimer.rev_min < 3:
                            self._3prim_dimers += 1
                            if self._min_3prim_dimer_dG == None \
                            or self._min_3prim_dimer_dG > min_dG:
                                self._min_3prim_dimer_dG = min_dG
                                self._min_3prim_dimer    = min_dG_dimer 
        return all_structures
    #end def


    def _find_hairpins(self, dimers, seq):
        """Find stable hairpins in given set of SYMMETRIC SELF-dimers"""
        seq_len = len(seq)
        hairpins   = dict()
        if not dimers: return hairpins
        for dimer in dimers.values():
            fwd_matches = set(dimer.fwd_matches)
            rev_matches = set(dimer.fwd_matches)
            while (rev_matches and fwd_matches)       and \
                  (min(rev_matches) < max(fwd_matches) or \
                   min(rev_matches) - max(fwd_matches) < 4): #at least 3 bases separation
                fwd_matches.remove(max(fwd_matches))
                rev_matches.remove(min(rev_matches))
            if fwd_matches:
                #find the most stable
                min_dG = max_hairpin_dG
                min_dG_hairpin = None
                dimer = Dimer(fwd_matches, min(fwd_matches)-(seq_len-1-max(rev_matches)))
                new_dimers = [dimer,] + dimer.alternatives
                for new_dimer in new_dimers:
                    hairpin = Hairpin(new_dimer.fwd_matches, set(seq_len-1-m for m in new_dimer.rev_matches))
                    dG = tdf.hairpin_dG_corrected(hairpin, seq)
                    hairpin.dG = dG
                    #find the most stable
                    if min_dG > dG:
                        min_dG = dG
                        min_dG_hairpin = hairpin
                #if there is a stable structure, append it
                if min_dG_hairpin:
                    h_hash = self._structure_hash(min_dG_hairpin)
                    if h_hash not in hairpins:
                        hairpins[h_hash] = min_dG_hairpin
                        #check for 3' structure
                        if min_dG_hairpin.rev_matches[0] > seq_len-4:
                            self._3prim_hairpins += 1
                            if self._min_3prim_hairpin_dG == None \
                            or self._min_3prim_hairpin_dG > min_dG:
                                self._min_3prim_hairpin_dG = min_dG
                                self._min_3prim_hairpin    = min_dG_hairpin 
        return hairpins
    #end def

    
    #formatting for secondary structures
    @classmethod
    def _format_dimer(cls, dimer, seq1, seq2):
        dimer_string = str(Duplex.print_dimer(seq1, seq2, dimer))
        #check for 3' dimer
        if dimer.fwd_matches[-1] > len(seq1)-4 or \
           dimer.rev_min < 3:
            dimer_string += '(3\'-dimer)\n'
        dimer_string += '\n'
        return dimer_string
    #end def


    def _format_dimers_header(self, dimers):
        header_string  = ''
        if not dimers: return 'No dimers found.\n\n'
        #dimers count
        header_string += 'Dimers count: %d\n' % len(dimers)
        #3'-count
        if self._3prim_dimers: 
            header_string += '3\'-dimers count: %d\n' % self._3prim_dimers
        #minimum dG
        header_string += 'Most stable dimer: %.2f kcal/mol\n' % self.dimerMin_dG()
        if self._3prim_dimers: 
            header_string += 'Most stable 3\' dimer: %.2f kcal/mol\n' % self.dimerMin_3prim_dG()
        header_string += '\n'
        return header_string
    #end def


    def _format_dimers_short(self, dimers, seq1, seq2):
        """short plain text representation of dimers"""
        dimers_string  = ''
        if not dimers: return 'No dimers found.\n\n'
        #print header
        dimers_string += self._format_dimers_header(dimers)
        #print the most stable dimer
        if self.dimerMin_dG():
            if dimers[0] != self._min_3prim_dimer: #if it is the 3' the next block will print it
                dimers_string += self._format_dimer(dimers[0], seq1, seq2)
        #print the most stable 3'-dimer
        if self.dimerMin_3prim_dG():
            dimers_string += self._format_dimer(self._min_3prim_dimer, seq1, seq2)
        return dimers_string
    #end def


    def _format_dimers(self, dimers, seq1, seq2):
        """plain text representation of dimers"""
        dimers_string  = ''
        if not dimers: return 'No dimers found.\n\n'
        #print header
        dimers_string += self._format_dimers_header(dimers)
        #print dimers
        for dimer in dimers:
            dimers_string += self._format_dimer(dimer, seq1, seq2)
        return dimers_string
    #end def


    @classmethod
    def _format_hairpin(cls, hairpin, seq):
        hairpin_string = ''
        fwd_matches = hairpin.fwd_matches
        rev_matches = hairpin.rev_matches
        loop        = (rev_matches[-1] - fwd_matches[-1]) - 1
        odd_loop    = loop%2
        half_loop   = (loop-odd_loop)/2
        fwd_tail    = fwd_matches[0] 
        rev_tail    = len(seq) - rev_matches[0] - 1
        break_pos   = fwd_matches[-1] + half_loop + 1
        #construct matches string
        matches = ''
        for m in xrange(len(fwd_matches)):
            if matches: matches += ' '*(fwd_matches[m]-fwd_matches[m-1]-1)
            matches += '|'
        #calculate spacers
        if fwd_tail > rev_tail:
            spacer1 = ''
            spacerM = ' '*(fwd_tail+3)
            spacer2 = ' '*(fwd_tail-rev_tail)
        else:
            spacer1 = ' '*(rev_tail-fwd_tail)
            spacerM = ' '*(rev_tail+3)
            spacer2 = ''
        loopM = ' '*(half_loop)
        link  = seq[break_pos] if odd_loop else ']'
        #construct hairpins string
        hairpin_string += "%(spacer1)s5' %(fwd)s\n" \
            % {'spacer1':spacer1,
               'fwd'    :str(seq[:break_pos])}
        hairpin_string += "%(spacerM)s%(matches)s%(loopM)s%(link)s\n"  \
            % {'spacerM':spacerM,
               'loopM'  :loopM,
               'matches':matches,
               'link'   :link}
        hairpin_string += "%(spacer2)s3' %(rev)s\n" \
            % {'spacer2':spacer2,
               'rev'    :str((seq[break_pos+odd_loop:])[::-1])}
        #dG
        if hairpin.dG != None:
            hairpin_string += 'dG(%.1fC) = %.2f kcal/mol\n' % (tdf.PCR_P.PCR_T, hairpin.dG)
        #conversion degree
        if hairpin.conversion_degree != None:
            hairpin_string += 'conversion degree = %.4f%%\n' % (hairpin.conversion_degree*100)
        #check for 3' hairpin
        if rev_matches[0] > len(seq)-4:
            hairpin_string += '(3\'-hairpin)\n'
        hairpin_string += '\n'
        return hairpin_string
    #end def


    def _format_hairpins_header(self, hairpins):
        header_string  = ''
        #hairpins count
        header_string += 'Hairpins count: %d\n' % len(hairpins)
        if self._3prim_hairpins: 
            header_string += '3\'-hairpins count: %d\n' % self._3prim_hairpins
        #minimum dG
        header_string += 'Most stable hairpin: %.2f kcal/mol\n' % self.hairpinMin_dG()
        if self._3prim_hairpins: 
            header_string += 'Most stable 3\' hairpin: %.2f kcal/mol\n' % self.hairpinMin_3prim_dG()
        header_string += '\n'
        return header_string
    #end def


    def _format_hairpins_short(self, hairpins, seq):
        hairpins_string  = ''
        if not hairpins: return 'No hairpins found\n\n'
        #print header
        hairpins_string += self._format_hairpins_header(hairpins)
        if self.hairpinMin_dG():
            if hairpins[0] != self._min_3prim_hairpin: #if it is the 3' the next block will print it
                hairpins_string += self._format_hairpin(hairpins[0], seq)
        #print the most stable hairpin
        if self.hairpinMin_3prim_dG():
            hairpins_string += self._format_hairpin(self._min_3prim_hairpin, seq)
        return hairpins_string
    #end def


    def _format_hairpins(self, hairpins, seq):
        hairpins_string  = ''
        if not hairpins: return 'No hairpins found\n\n'
        #print header
        hairpins_string += self._format_hairpins_header(hairpins)
        #print hairpins
        for hairpin in hairpins:
            hairpins_string += self._format_hairpin(hairpin, seq)
        return hairpins_string
    #end def
#end class