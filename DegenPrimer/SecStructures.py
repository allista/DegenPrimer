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
Created on Jun 23, 2012

@author: Allis Tauri <allista@gmail.com>
'''
import TD_Functions
from TD_Functions import dimer_dG_corrected, hairpin_dG_corrected, C_Prim_M, equilibrium_constant
from StringTools import hr
from Equilibrium import Equilibrium

    
def all_combinations(lst):
    if not lst: return []
    if len(lst) == 1: return [[lst[0]]]
    combinations = []
    combinations.append([lst[0]])
    for comb in all_combinations(lst[1:]):
        combinations.append(comb)
        combinations.append([lst[0]]+comb)
    return combinations
#end def


class Dimer(object):
    '''
    Structure of a dimer:
    forward matches is a list of nucleotide indices of forward sequence directed 5'->3'
    reverse matches is a list of nucleotide indices of reverse sequence directed 3'->5'
    nucleotide fwd_match[n] canonically matches nucleotide rev_matches[n]
    e.g. fwd_matches = 5'-[2,4,7, 9]-3'
         rev_matches = 3'-[4,6,9,11]-5'
    '''
    def __init__(self, fwd_matches=None, rev_matches=None, dG=0):
        if fwd_matches:
            self._fwd_matches = set(fwd_matches)
        else:
            self._fwd_matches = set()    
        if rev_matches:
            self._rev_matches = set(rev_matches)
        else:
            self._rev_matches = set()
        self.dG = None
        self.conversion_degree = None
    
    def __nonzero__(self):
        if  len(self._fwd_matches) > 0 \
        and len(self._rev_matches) == len(self._fwd_matches):
            return True
        else: return False
    
    def __getitem__(self, i):
        if   i == 0: return self._fwd_matches
        elif i == 1: return self._rev_matches
        elif i == 2: return self.dG
        else: raise IndexError('SecStructures.Dimer: index out of bounds.')
        
    def fwd_matches(self):
        l = list(self._fwd_matches)
        l.sort()
        return l
    
    def rev_matches(self):
        l = list(self._rev_matches)
        l.sort()
        return l
        
    def add(self, fwd, rev):
        self._fwd_matches.add(fwd)
        self._rev_matches.add(rev)
#end class


class Hairpin(Dimer):
    '''
    Structure of a hairpin:
    described sequence directed 5'->3' 
    forward matches is a list of near-5'-end nucleotide indices
    reverse matches is a *reversed* list of near-3'-end nucleotide indices
    nucleotide fwd_match[n] canonically matches nucleotide rev_matches[n]
    e.g. fwd_matches = 5'-[2, 4, 7,  9]...-3'
         rev_matches = 3'-[22,20,17,15]...-5'
    '''        
    def rev_matches(self):
        l = list(self._rev_matches)
        l.sort(reverse=True)
        return l 
#end class


class SecStructures(object):
    '''
    Given one or two sequence objects compatible with BioPython SeqRecord objects, 
    computes sets of self-dimers, hairpins and (in case of two sequences) cross-dimers.
    To get formatted output just convert to a str: 
    >>> print str(MySecStructs)
    '''
    
    max_dimer_dG   = -3 #kcal/mol #corresponds to conversion degree ~1e-3
    max_hairpin_dG =  3 #kcal/mol #corresponds to conversion degree ~1e-2

    def __init__(self, seq_rec1, seq_rec2=None):
        '''
        seq_rec1 -- Bio.SeqRecord object containing first sequence
        seq_rec2 -- Bio.SeqRecord object containing second sequence. May be None. 
        If so, hairpins and self-dimers for seq_req1 will be calculated.
        Otherwise, only corss-dimers between seq_req1 and seq_req2 will be calculated.
        conversion_threshold -- secondary structures with conversion degree below this value 
        will not be included into the short report. (default 0.01, i.e. 1%)
        '''
        self._seq_rec1 = seq_rec1
        self._seq1 = seq_rec1.seq
        if seq_rec2 and str(seq_rec1.seq) != str(seq_rec2.seq):
            self._seq_rec2 = seq_rec2
            self._seq2 = seq_rec2.seq
        else: 
            self._seq_rec2 = None
            self._seq2 = None
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
                    reactions[d_hash] = Equilibrium.compose_reaction(equilibrium_constant(D.dG, TD_Functions.PCR_T), 
                                                                     str(self._seq1), 
                                                                     None, 
                                                                     '2A')
            if self._seq1_hairpins:
                for h_hash in self._seq1_hairpins:
                    H = self._seq1_hairpins[h_hash]
                    reactions[h_hash] = Equilibrium.compose_reaction(equilibrium_constant(H.dG, TD_Functions.PCR_T), 
                                                                     str(self._seq1), 
                                                                     None, 
                                                                     'A')
        else: #crossdimers
            if self._cross_dimers:
                for d_hash in self._cross_dimers:
                    D = self._cross_dimers[d_hash]
                    reactions[d_hash] = Equilibrium.compose_reaction(equilibrium_constant(D.dG, TD_Functions.PCR_T), 
                                                                     str(self._seq1), 
                                                                     str(self._seq2), 
                                                                     'AB')
        return reactions
    #end def
    
    
    def set_conversion_degree(self, s_hash, r):
        if not self._seq_rec2:
            if   s_hash in self._seq1_dimers:
                self._seq1_dimers[s_hash].conversion_degree   = r
            elif s_hash in self._seq1_hairpins:
                self._seq1_hairpins[s_hash].conversion_degree = r
            else: return False
        else:
            if s_hash in self._cross_dimers:
                self._cross_dimers[s_hash].conversion_degree  = r
            else: return False
        return True
    #end def
    
    
    def max_conversion_degree(self):
        if self._seq2 and self._cross_dimers:
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
        if self._seq2 and self._cross_dimers:
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
        if self._seq2 or not self._seq1_hairpins: 
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
                         tuple(struct.fwd_matches()+struct.rev_matches())))
        else:
            return hash((str(self._seq1), str(self._seq2), 
                         tuple(struct.fwd_matches()+struct.rev_matches())))
    #end def
    
    
    @classmethod
    def _remove_loops(cls, dimer):
        """remove terminal internal loops from a dimer one by one in all 
        possible combinations to produce a bunch of new structures"""
        dimers  = [dimer]
        fwd_matches = dimer.fwd_matches()
        rev_matches = dimer.rev_matches()
        #scan for loops from both ends
        left_loops  = []
        right_loops = []
        left_margin  = (set(),set())
        right_margin = (set(),set())
        for i,j in zip(range(len(fwd_matches)-1), range(len(fwd_matches)-1, 0, -1)):
            #left
            left_margin[0].add(fwd_matches[i])
            left_margin[1].add(rev_matches[i])
            if fwd_matches[i+1]-fwd_matches[i] > 1:
                left_loops.append((left_margin, fwd_matches[i+1]-fwd_matches[i]-1))
                left_margin  = (set(),set())
            #right
            #matches and gaps
            right_margin[0].add(fwd_matches[j])
            right_margin[1].add(rev_matches[j])
            if fwd_matches[j]-fwd_matches[j-1] > 1:
                right_loops.append((right_margin, fwd_matches[j]-fwd_matches[j-1]-1))
                right_margin = (set(),set())
        #of no loops have been found, return original structure
        if not left_loops: return dimers
        #else, remove boundary loops one by one
        l_dimer = dimer
        for l in range(len(left_loops)):
            r_dimer = l_dimer
            for r in range(len(right_loops)-l):
                r_dimer = Dimer(r_dimer[0]-right_loops[r][0][0],
                               r_dimer[1]-right_loops[r][0][1])
                dimers.append(r_dimer)
            l_dimer = Dimer(l_dimer[0]-left_loops[l][0][0],
                           l_dimer[1]-left_loops[l][0][1])
            dimers.append(l_dimer)
        return dimers
    #end def
    
    
    def _find_dimers(self, seq1, seq2):
        """Search for all dimers"""
        all_structures = {'dimers': dict(), 'hairpins': dict()}
        self_dimers = str(seq1) == str(seq2)
        dimers   = all_structures['dimers']
        hairpins = all_structures['hairpins']
        fwd_str  = str(seq1)
        rev_str  = str(seq2.reverse_complement())
        fwd_len  = len(fwd_str)
        rev_len  = len(rev_str)
        full_len = fwd_len+rev_len-1
        for front in range(full_len):
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
            for fwd,rev in zip(range(fwd_left,fwd_right),range(rev_left,rev_right)):
                if(fwd_str[fwd] == rev_str[rev]):
                    dimer.add(fwd, rev)
            #if found some
            if dimer[0]:
                #add a copy to the hairpins dict to search for hairpin structures later
                if self_dimers and dimer not in hairpins: 
                    hairpins[self._structure_hash(dimer)] = dimer
                #find the most stable dimer
                min_dG = self.max_dimer_dG
                min_dG_dimer = None
                new_dimers = self._remove_loops(dimer)
                for dimer in new_dimers:
                    #calculate dG for every dimer
                    dG = dimer_dG_corrected(dimer, seq1, seq2)
                    dimer.dG = dG
                    #find the most stable
                    if min_dG > dG:
                        min_dG = dG
                        min_dG_dimer = dimer
                #if there is a stable dimer, append it
                if min_dG_dimer:
                    d_hash = self._structure_hash(min_dG_dimer)
                    if d_hash not in dimers:
                        dimers[d_hash] = min_dG_dimer
                        #check for 3' dimer
                        if max(min_dG_dimer.fwd_matches()) > fwd_len-4 \
                        or min(min_dG_dimer.rev_matches()) < 3:
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
            fwd_matches = dimer.fwd_matches()
            rev_matches = dimer.fwd_matches()
            while (rev_matches and fwd_matches)       and \
                  (min(rev_matches) < max(fwd_matches) or \
                   min(rev_matches) - max(fwd_matches) < 4): #at least 3 bases separation
                fwd_matches.remove(max(fwd_matches))
                rev_matches.remove(min(rev_matches))
            if fwd_matches:
                #find the most stable
                min_dG = self.max_hairpin_dG
                min_dG_hairpin = None
                new_dimers = self._remove_loops(Dimer(fwd_matches, set(seq_len-1-m for m in rev_matches)))
                for new_dimer in new_dimers:
                    hairpin = Hairpin(new_dimer[0], set(seq_len-1-m for m in new_dimer[1]))
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
                        if max(min_dG_hairpin.rev_matches()) > seq_len-4:
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
        dimer_string = ''
        fwd_matches = dimer.fwd_matches()
        rev_matches = dimer.rev_matches()
        fwd_matches.sort()
        rev_matches.sort()
        #construct matches string
        matches = ''
        for m in range(len(fwd_matches)):
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
        dimer_string += "%(spacer1)s5' %(seq1)s 3'\n" \
            % {'spacer1':spacer1,
               'seq1'   :str(seq1)}
        dimer_string += "%(spacerM)s%(matches)s\n"    \
            % {'spacerM':spacerM,
               'matches':matches}
        dimer_string += "%(spacer2)s3' %(seq2)s 5'\n" \
            % {'spacer2':spacer2,
               'seq2'   :str(seq2[::-1])}
        #dG
        if dimer.dG != None:
            dimer_string += 'dG(%.1fC) = %.2f kcal/mol\n' % (TD_Functions.PCR_T, dimer.dG)
        #conversion degree
        if dimer.conversion_degree != None:
            dimer_string += 'conversion degree = %.4f%%\n' % (dimer.conversion_degree*100)
        #check for 3' dimer
        if fwd_matches[-1] > len(seq1)-4 or \
           rev_matches[0]  < 3:
            dimer_string += '(3\'-dimer)\n'
        dimer_string += '\n'
        return dimer_string
    #end def


    def _format_dimers_header(self, dimers, seq1, seq2):
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
        dimers_string += self._format_dimers_header(dimers, seq1, seq2)
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
        dimers_string += self._format_dimers_header(dimers, seq1, seq2)
        #print dimers
        for dimer in dimers:
            dimers_string += self._format_dimer(dimer, seq1, seq2)
        return dimers_string
    #end def


    @classmethod
    def _format_hairpin(cls, hairpin, seq):
        hairpin_string = ''
        fwd_matches = hairpin.fwd_matches()
        rev_matches = hairpin.rev_matches()
        loop        = (min(rev_matches) - max(fwd_matches)) - 1
        fwd_loop    = loop/2
        rev_loop    = loop - fwd_loop
        fwd_tail    = min(fwd_matches) 
        rev_tail    = len(seq) - max(rev_matches) - 1
        break_pos   = max(fwd_matches) + fwd_loop + 1
        #construct matches string
        matches = ''
        for m in range(len(fwd_matches)):
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
        loop1   = '-'*(rev_loop-fwd_loop)
        loopM   = ' '*rev_loop    
        #construct hairpins string
        hairpin_string += "%(spacer1)s5' %(fwd)s%(loop1)s\n" \
            % {'spacer1':spacer1,
               'loop1'  :loop1,
               'fwd'    :str(seq[:break_pos])}
        hairpin_string += "%(spacerM)s%(matches)s%(loopM)s]\n"  \
            % {'spacerM':spacerM,
               'loopM'  :loopM,
               'matches':matches}
        hairpin_string += "%(spacer2)s3' %(rev)s\n" \
            % {'spacer2':spacer2,
               'rev'    :str((seq[break_pos:])[::-1])}
        #dG
        if hairpin.dG != None:
            hairpin_string += 'dG(%.1fC) = %.2f kcal/mol\n' % (TD_Functions.PCR_T, hairpin.dG)
        #conversion degree
        if hairpin.conversion_degree != None:
            hairpin_string += 'conversion degree = %.4f%%\n' % (hairpin.conversion_degree*100)
        #check for 3' hairpin
        if rev_matches[0] > len(seq)-4:
            hairpin_string += '(3\'-hairpin)\n'
        hairpin_string += '\n'
        return hairpin_string
    #end def


    def _format_hairpins_header(self, hairpins, seq):
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
        hairpins_string += self._format_hairpins_header(hairpins, seq)
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
        hairpins_string += self._format_hairpins_header(hairpins, seq)
        #print hairpins
        for hairpin in hairpins:
            hairpins_string += self._format_hairpin(hairpin, seq)
        return hairpins_string
    #end def
#end class


class AllSecStructures(object):
    '''
    Given two lists of non-degenerate primers calculate all possible secondary 
    structures (hairpins, dimers, cross-dimers) and equilibrium state between them.
    '''
    def __init__(self, fwd_primers, rev_primers):
        #initial check
        if not fwd_primers and not rev_primers:
            raise ValueError('AllSecStructures:__init__: at least one list of primers should be provided.')
        self._fwd_primers = fwd_primers
        self._rev_primers = rev_primers
        #all primers list
        self._all_primers = list()
        if fwd_primers: self._all_primers += fwd_primers
        if rev_primers: self._all_primers += rev_primers
        #reactions
        self._reactions = dict()
        #secondary structures and reactions of their formation
        self._fwd_self = []
        self._rev_self = []
        if self._fwd_primers:
            #selfdimers and hairpins
            self._fwd_self = [SecStructures(p) for p in self._fwd_primers]
            for struct in self._fwd_self:
                self._reactions.update(struct.compose_reactions())
        if self._rev_primers:
            #selfdimers and hairpins
            self._rev_self = [SecStructures(p) for p in self._rev_primers]
            for struct in self._rev_self:
                self._reactions.update(struct.compose_reactions())
        #crossdimers
        self._cross = []
        for i in range(len(self._all_primers)):
            for j in range(i+1, len(self._all_primers)):
                self._cross.append(SecStructures(self._all_primers[i],
                                                 self._all_primers[j]))
        for struct in self._cross:
                self._reactions.update(struct.compose_reactions())
        #concentrations
        self._concentrations = dict().fromkeys([str(p.seq) for p in self._all_primers], C_Prim_M())
        #equilibrium system
        self._equilibrium = None
        self._equilibrium_concentrations = None
    #end def
    
    #property functions
    def reactions(self): return self._reactions
    def concentrations(self): return self._concentrations
    def equilibrium_concentrations(self): return self._equilibrium_concentrations
    
    
    def calculate_equilibrium(self):
        if not self._reactions: return False
        #calculate equilibrium
        self._equilibrium = Equilibrium(self._reactions, self._concentrations)
        sol = self._equilibrium.calculate()[0]
        #calculate equilibrium primer concentrations
        self._equilibrium_concentrations = dict()
        for primer in self._concentrations:
            _C, _None = self._equilibrium.reactants_consumtion(primer)
            self._equilibrium_concentrations[primer] = self._concentrations[primer] - _C
        #set conversion degrees
        for r_hash in sol:
            for struct in (self._fwd_self+self._rev_self+self._cross):
                if struct.set_conversion_degree(r_hash, sol[r_hash]): break
    #end def
    
    
    @classmethod
    def _print_structures(cls, structures, full=True):
        structures_string = ''
        for struct in structures:
            if not struct: continue
            if full: structures_string += struct.formatFull()
            else: structures_string += struct.formatShort()
        return structures_string
    #end def
    
    def print_structures(self, full=True):
        structures_string = ''
        if self._fwd_self:
            fwd_structs = self._print_structures(self._fwd_self, full)
            if fwd_structs:
                structures_string += hr(' secondary structures of forward primers ', '#')
                structures_string += fwd_structs 
        if self._rev_self:
            rev_structs = self._print_structures(self._rev_self, full)
            if rev_structs:
                structures_string += hr(' secondary structures of reverse primers ', '#')
                structures_string += rev_structs 
        if self._cross:
            cross_structs = self._print_structures(self._cross, full)
            if cross_structs:
                structures_string += hr(' cross-dimers ', '#')
                structures_string += cross_structs
        if structures_string:
            structures_string = hr(' stable secondary structures ', symbol='#') + structures_string
        return structures_string
    #end def
    
    def print_structures_short(self):
        return self.print_structures(full=False)
#end class


#tests
if __name__ == '__main__':
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    TD_Functions.PCR_T = 55
    
    fwd_primers = [SeqRecord(Seq('ATATTCTACGACGGCTATCC').reverse_complement()), SeqRecord(Seq('ATATGCTACGACGGCTATCC').reverse_complement())]
    rev_primers = [SeqRecord(Seq('CAAGGGTTAGAAGCGGAAG').complement())]
    all_structs = AllSecStructures(fwd_primers, rev_primers)

    print all_structs.concentrations
    print all_structs.reactions
    
    all_structs.calculate_equilibrium()
    print all_structs.equilibrium_concentrations
    print ''
    print all_structs.print_structures()