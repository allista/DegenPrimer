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
from Bio.Seq import Seq
import TD_Functions
from TD_Functions import dimer_dG_corrected, hairpin_dG_corrected, equilibrium_constant
from StringTools import hr
from Equilibrium import Reaction


#default filtration parameters
max_dimer_dG   = -3 #kcal/mol #corresponds to equilibrium constant ~100  and conversion degree ~1e-3
max_hairpin_dG =  3 #kcal/mol #corresponds to equilibrium constant ~0.01 and conversion degree ~1e-2


class Dimer(object):
    '''
    Structure of a dimer without bulges (i.e. asymmetrical loops):
    forward matches is a list of nucleotide indices of forward sequence 
    directed 5'->3' that are paired to reverse strand.
    reverse matches is a list of nucleotide indices of reverse 
    directed 3'->5' that are paired to reverse strand.
    Nucleotide fwd_match[n] canonically matches nucleotide rev_matches[n]
    e.g. fwd_matches = 5'-[2,4,7, 9]-3'
         rev_matches = 3'-[4,6,9,11]-5'
    '''
    __slots__ = ['fwd_matches', '_offset', 'num_matches',
                 'conversion_degree', 'dG', 'K',
                 'fwd_mismatch']
    
    def __init__(self, fwd_matches=None, offset=0): 
        #dG, K and CD are properties of a Duplex, not a Dimer, but since 
        #a Duplex may have several associated dimers, it's convenient to store them here 
        if fwd_matches:
            self.fwd_matches  = tuple(sorted(list(fwd_matches)))
            self._offset      = offset
            self.num_matches  = len(self.fwd_matches)
        else:
            self.fwd_matches  = tuple()
            self._offset      = None
            self.num_matches  = 0    
        #thermodynamic parameters
        self.conversion_degree = None
        self.dG = None
        self.K  = None
        #topological parameters
        self.fwd_mismatch   = None
    #end def
    
    def __hash__(self):
        return hash((self.fwd_matches, self._offset))
    
    def __eq__(self, other):
        return self.__hash__() == other.__hash__()
    
    def __ne__(self, other):
        return self.__hash__() != other.__hash__()
    
    def __nonzero__(self):
        if self.num_matches > 0: return True
        else: return False
    
    
    @property
    def rev_matches(self):
        return tuple(sorted(list(i - self._offset for i in self.fwd_matches)))
    
    @property
    def fwd_min(self):
        return min(self.fwd_matches)
    
    @property
    def rev_min(self):
        return min(self.fwd_matches)-self._offset
    
    @property
    def fwd_max(self):
        return max(self.fwd_matches)
    
    @property
    def rev_max(self):
        return max(self.fwd_matches)-self._offset
    
    
    @property
    def alternatives(self):
        """remove terminal internal loops from a dimer one by one in all 
        possible combinations to produce a bunch of new structures"""
        #scan for loops from both ends
        left_loops   = []
        right_loops  = []
        left_margin  = set()
        right_margin = set()
        for i,j in zip(xrange(self.num_matches-1), xrange(self.num_matches-1, 0, -1)):
            #left
            left_margin.add(self.fwd_matches[i])
            if self.fwd_matches[i+1]-self.fwd_matches[i] > 1:
                left_loops.append(left_margin)
                left_margin  = set()
            #right
            #matches and gaps
            right_margin.add(self.fwd_matches[j])
            if self.fwd_matches[j]-self.fwd_matches[j-1] > 1:
                right_loops.append(right_margin)
                right_margin = set()
        #of no loops have been found, return empty list
        if not left_loops: return []
        #else, remove boundary loops one by one
        dimers  = []
        l_dimer = self
        for l in xrange(len(left_loops)):
            r_dimer = l_dimer
            for r in xrange(len(right_loops)-l):
                r_dimer = Dimer(set(r_dimer.fwd_matches)-right_loops[r], self._offset)
                dimers.append(r_dimer)
            l_dimer = Dimer(set(l_dimer.fwd_matches)-left_loops[l], self._offset)
            dimers.append(l_dimer)
        return dimers
    #end def
       
    
    def add(self, fwd, rev):
        if self._offset is None:
            self._offset = fwd-rev
        elif self._offset != fwd-rev:
            raise ValueError('Dimer.add: fwd and rev indexes should be at the '
                             'same distance as all other match pairs. '
                             'Asymmetrical loops are not allowed.')
        self.fwd_matches += (fwd,)
        self.num_matches += 1
    #end def
    
        
    @classmethod
    def from_sequences(cls, seq1, seq2, anchor=None):
        '''Form a dimer structure from two sequences using a pair of indexes 
        as an anchor; if anchor is not provided, sequences are left-aligned;
        both sequences should be 5'->3' oriented'''
        forward  = seq1
        revcomp  = str(Seq(seq2).complement())[::-1]
        dimer    = cls()
        offset   = (anchor[0]-anchor[1]) if anchor else 0
        dimer._offset = offset
        dimer.fwd_matches = []
        if offset > 0: _range = range(offset, min(len(seq1), len(seq2)+offset))
        else: _range = range(0, min(len(seq1)-offset, len(seq2)+offset))
        for i in _range:
            if forward[i] == revcomp[i-offset]:
                #add indexes directly is much faster than to call Dimer.add method
                dimer.fwd_matches.append(i)
                dimer.num_matches += 1
        dimer.fwd_matches = tuple(dimer.fwd_matches)
        return dimer
    #end def
#end class


class Hairpin(object):
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
        self.fwd_matches = tuple(sorted(list(fwd_matches)))
        self.rev_matches = tuple(sorted(list(rev_matches), reverse=True))
        #thermodynamic parameters
        self.conversion_degree = None
        self.dG = None
        self.K  = None
    #end def
#end class


class Duplex(object):
    '''Representation of a particular DNA:DNA duplex'''
    __slots__ = ['_fwd_sequence', '_rev_sequence',
                 'fwd_len', 'rev_len', 
                 '_dimers',
                 '_nonzero',
                 '_fwd_matches', '_fwd_mismatches', 'fwd_3_overhang']
    
    def __init__(self, fwd_seq, rev_seq, dimer=None):
        '''
        fwd_seq -- Seq object, 5'->3' sequence of forward strand
        rev_seq -- Seq object, 5'->3' sequence of reverse strand
        '''
        self._fwd_sequence   = fwd_seq
        self._rev_sequence   = rev_seq
        self.fwd_len         = len(fwd_seq)
        self.rev_len         = len(rev_seq)
        if not dimer: dimer  = Dimer.from_sequences(fwd_seq, rev_seq)
        self._dimers         = [dimer]+dimer.alternatives
        self._fwd_matches    = []
        self._fwd_mismatches = []
        self.fwd_3_overhang  = None
        self._nonzero        = False 
        self._recalculate()
    #end def
    
    
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
            dimer.dG = dimer_dG_corrected(dimer, 
                                          self._fwd_sequence,
                                          self._rev_sequence)
            if dimer.dG > max_dimer_dG: continue
            dimer.K  = equilibrium_constant(dimer.dG, 
                                            TD_Functions.PCR_P.PCR_T)
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
    

    def filter_dimers_by_K(self, min_K):
        if len(self._dimers) < 2: return
        new_dimers = [self._dimers[0]]
        for i, dimer in enumerate(self._dimers[1:]):
            if dimer.K >= min_K:
                new_dimers.append(dimer)
            else:
                try: self._fwd_matches.remove(i+1)
                except: self._fwd_mismatches.remove(i+1)
        self._dimers = tuple(new_dimers)
    #end def
    
    @classmethod
    def print_dimer(cls, fwd_sequence, rev_sequence, dimer):
        duplex_string = ''
        if not dimer:
            duplex_string += fwd_sequence + '\n'
            duplex_string += rev_sequence + '\n'
            return duplex_string
        fwd_matches = dimer.fwd_matches
        rev_matches = dimer.rev_matches
        #construct matches string
        matches = ''
        for m in xrange(dimer.num_matches):
            if matches: matches += ' '*(fwd_matches[m]-fwd_matches[m-1]-1)
            matches += '|'            
        #calculate spacers
        if fwd_matches[0] > rev_matches[0]:
            spacer1 = ''
            spacerM = ' '*(fwd_matches[0]+3)
            spacer2 = ' '*(fwd_matches[0]-rev_matches[0])
        else:
            spacer1 = ' '*(rev_matches[0]-fwd_matches[0])
            spacerM = ' '*(rev_matches[0]+3)
            spacer2 = ''
        #construct structures string
        duplex_string += "%(spacer1)s5' %(seq1)s 3'\n" \
            % {'spacer1':spacer1,
               'seq1'   :fwd_sequence}
        duplex_string += "%(spacerM)s%(matches)s\n"    \
            % {'spacerM':spacerM,
               'matches':matches}
        duplex_string += "%(spacer2)s3' %(seq2)s 5'\n" \
            % {'spacer2':spacer2,
               'seq2'   :rev_sequence[::-1]}
        #dG
        if dimer.dG != None:
            duplex_string += 'dG(%.1fC) = %.2f kcal/mol\n' % (TD_Functions.PCR_P.PCR_T, 
                                                              dimer.dG)
        #conversion degree
        if dimer.conversion_degree != None:
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
                    reactions[d_hash] = Reaction(equilibrium_constant(D.dG, TD_Functions.PCR_P.PCR_T), 
                                                 str(self._seq1), 
                                                 None, 
                                                 '2A')
            if self._seq1_hairpins:
                for h_hash in self._seq1_hairpins:
                    H = self._seq1_hairpins[h_hash]
                    reactions[h_hash] = Reaction(equilibrium_constant(H.dG, TD_Functions.PCR_P.PCR_T), 
                                                 str(self._seq1), 
                                                 None, 
                                                 'A')
        else: #crossdimers
            if self._cross_dimers:
                for d_hash in self._cross_dimers:
                    D = self._cross_dimers[d_hash]
                    reactions[d_hash] = Reaction(equilibrium_constant(D.dG, TD_Functions.PCR_P.PCR_T), 
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
        rev_str  = str(seq2.reverse_complement())
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
                    dimer.dG = dimer_dG_corrected(dimer, seq1, seq2)
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
                        if max(min_dG_dimer.fwd_matches) > fwd_len-4 \
                        or min(min_dG_dimer.rev_matches) < 3:
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
                    dG = hairpin_dG_corrected(hairpin, seq)
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
                        if max(min_dG_hairpin.rev_matches) > seq_len-4:
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
           dimer.rev_matches[0]  < 3:
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
        loop        = (min(rev_matches) - max(fwd_matches)) - 1
        odd_loop    = loop%2
        half_loop   = (loop-odd_loop)/2
        fwd_tail    = min(fwd_matches) 
        rev_tail    = len(seq) - max(rev_matches) - 1
        break_pos   = max(fwd_matches) + half_loop + 1
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
            hairpin_string += 'dG(%.1fC) = %.2f kcal/mol\n' % (TD_Functions.PCR_P.PCR_T, hairpin.dG)
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


#tests
if __name__ == '__main__':
    from tests.asizeof import heappy

    d = Dimer()
    h = Hairpin([], [])
    du = Duplex('A', 'T')
    
    print heappy.iso(d)
    print heappy.iso(h)
    print heappy.iso(du)

    print heappy.iso(set())
    print heappy.iso(list())
    print heappy.iso(tuple())