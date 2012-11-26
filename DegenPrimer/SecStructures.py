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
from TD_Functions import dimer_dG_corrected, hairpin_dG_corrected
from StringTools import hr

    
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
        self.dG = 0
    
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

    def __init__(self, seq_rec1, seq_rec2=None, dG_threshold=0):
        """Constructor"""
        self._seq_rec1 = seq_rec1
        self._seq1 = seq_rec1.seq
        if seq_rec2 and str(seq_rec1.seq) != str(seq_rec2.seq):
            self._seq_rec2 = seq_rec2
            self._seq2 = seq_rec2.seq
        else: 
            self._seq_rec2 = None
            self._seq2 = None
        self._dG_threshold = dG_threshold
        self._recalculate()
    #end def
    
    def belowThreshold(self):
        return self.dimerMin_dG() and self.dimerMin_dG() <= self._dG_threshold \
            or self.dimerMin_3prim_dG() and self.dimerMin_3prim_dG() <= self._dG_threshold+1 \
            or self.hairpinMin_dG() and self.hairpinMin_dG() <= self._dG_threshold+2 \
            or self.hairpinMin_3prim_dG() and self.hairpinMin_3prim_dG() <= self._dG_threshold+3
    #end def
    
    def dimerMin_dG(self):
        if self._seq2 and self._cross_dimers:
            return self._cross_dimers[0].dG
        elif self._seq1_dimers:
            return self._seq1_dimers[0].dG
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
        else: return self._seq1_hairpins[0].dG
    #end def
    
    def hairpinMin_3prim_dG(self):
        if self._min_3prim_hairpin_dG:
            return self._min_3prim_hairpin_dG
        else: return None
    #end def
    
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
            string += format_dimers(self._cross_dimers, self._seq1, self._seq2)
        else:
            string += hr(' %s: %s ' % (self._seq_rec1.id, str(self._seq1)), symbol='=')
            string += hr(' self-dimers ')
            string += format_dimers(self._seq1_dimers, self._seq1, self._seq1)
            string += hr(' hairpins ')
            string += format_hairpins(self._seq1_hairpins, self._seq1)
        string += hr('',symbol='=')
        return string
    
    def formatFull(self):
        return self._format_str(full=True)
    
    def formatShort(self):
        return self._format_str(full=False)
    
    def __str__(self):
        return self._format_str()
    
    def _recalculate(self):
        """For single sequence calculate self-dimers and hairpins.
        For two sequences calculate only cross-dimers"""
        self._3prim_dimers         = 0
        self._min_3prim_dimer_dG   = 0
        self._min_3prim_dimer      = None
        self._3prim_hairpins       = 0
        self._min_3prim_hairpin_dG = 0
        self._min_3prim_hairpin    = None 
        if self._seq2:
            self._seq1_dimers    = None
            self._seq1_hairpins  = None
            self._seq2_dimers    = None
            self._seq2_hairpins  = None
            self._cross_dimers   = self._find_dimers(self._seq1, self._seq2)['dimers']
        else:
            all_structures       = self._find_dimers(self._seq1, self._seq1)
            self._seq1_dimers    = all_structures['dimers']
            self._seq1_hairpins  = self._find_hairpins(all_structures['hairpins'], self._seq1)
            self._seq2_dimers    = None
            self._seq2_hairpins  = None
            self._cross_dimers   = None
    #end def        
    
    def _remove_loops(self, dimer):
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
        all_structures = {'dimers': [], 'hairpins': []}
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
                if self_dimers and dimer not in hairpins: 
                    hairpins.append(dimer)
                min_dG = 0
                min_dG_structure = None
                new_structures = self._remove_loops(dimer)
                for struct in new_structures:
                    #calculate dG for every dimer
                    dG = dimer_dG_corrected(struct, seq1, seq2)
                    struct.dG = dG
                    #find the most stable
                    if min_dG > dG:
                        min_dG = dG
                        min_dG_structure = struct
                #if there is a stable dimer, append it
                if min_dG_structure and min_dG_structure not in dimers:
                    dimers.append(min_dG_structure)
                    #check for 3' dimer
                    if max(min_dG_structure[0]) > len(seq1)-4 \
                    or min(min_dG_structure[1]) < 3:
                        self._3prim_dimers += 1
                        if self._min_3prim_dimer_dG > min_dG:
                            self._min_3prim_dimer_dG = min_dG
                            self._min_3prim_dimer    = min_dG_structure 
        #sort by dG
        dimers.sort(key=lambda x: x.dG)
        return all_structures
    #end def

    def _find_hairpins(self, dimers, seq):
        """Find stable hairpins in given set of SYMMETRIC SELF-dimers"""
        seq_len = len(seq)
        hairpins   = []
        if not dimers: return hairpins
        for dimer in dimers:
            fwd_matches = dimer.fwd_matches()
            rev_matches = dimer.fwd_matches()
            while (rev_matches and fwd_matches)       and \
                  (min(rev_matches) < max(fwd_matches) or \
                   min(rev_matches) - max(fwd_matches) < 4): #at least 3 bases separation
                fwd_matches.remove(max(fwd_matches))
                rev_matches.remove(min(rev_matches))
            if fwd_matches:
                min_dG = 0
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
                if min_dG_hairpin and min_dG_hairpin not in hairpins:
                    hairpins.append(min_dG_hairpin) 
                    #check for 3' structure
                    if min(min_dG_hairpin[0]) < 3:
                        self._3prim_hairpins += 1
                        if self._min_3prim_hairpin_dG > min_dG:
                            self._min_3prim_hairpin_dG = min_dG
                            self._min_3prim_hairpin    = min_dG_hairpin 
        hairpins.sort(key=lambda x: x.dG)
        return hairpins
    #end def

    def _format_dimer(self, dimer, seq1, seq2):
        structure_string = ''
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
        structure_string += "%(spacer1)s5' %(seq1)s 3'\n" \
            % {'spacer1':spacer1,
               'seq1'   :str(seq1)}
        structure_string += "%(spacerM)s%(matches)s\n"    \
            % {'spacerM':spacerM,
               'matches':matches}
        structure_string += "%(spacer2)s3' %(seq2)s 5'\n" \
            % {'spacer2':spacer2,
               'seq2'   :str(seq2[::-1])}
        #dG
        structure_string += 'dG(%.1fC) = %.2f kcal/mol\n' % (TD_Functions.PCR_T, dimer.dG)
        #check for 3' dimer
        if fwd_matches[-1] > len(seq1)-4 or \
           rev_matches[0]  < 3:
            structure_string += '(3\'-dimer)\n'
        structure_string += '\n'
        return structure_string
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
        header_string += 'Most stable dimer: %.2f kcal/mol\n' % dimers[0].dG
        if self._3prim_dimers: 
            header_string += 'Most stable 3\' dimer: %.2f kcal/mol\n' % self._min_3prim_dimer_dG
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
        if self.dimerMin_dG() and self.dimerMin_dG() < self._dG_threshold:
            if dimers[0] != self._min_3prim_dimer: #if it is the 3' the next block will print it
                dimers_string += self._format_dimer(dimers[0], seq1, seq2)
        #print the most stable 3'-dimer
        if self.dimerMin_3prim_dG() and self.dimerMin_3prim_dG() < self._dG_threshold+1:
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

    def _format_hairpin(self, hairpin, seq):
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
        hairpin_string += 'dG(%.1fC) = %.2f kcal/mol\n' % (TD_Functions.PCR_T, hairpin.dG)
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
        header_string += 'Most stable hairpin: %.2f kcal/mol\n' % hairpins[0].dG
        if self._3prim_hairpins: 
            header_string += 'Most stable 3\' hairpin: %.2f kcal/mol\n' % self._min_3prim_hairpin_dG
        header_string += '\n'
        return header_string
    #end def

    def _format_hairpins_short(self, hairpins, seq):
        hairpins_string  = ''
        if not hairpins: return 'No hairpins found\n\n'
        #print header
        hairpins_string += self._format_hairpins_header(hairpins, seq)
        if self.hairpinMin_dG() and self.hairpinMin_dG() < self._dG_threshold+2:
            if hairpins[0] != self._min_3prim_hairpin: #if it is the 3' the next block will print it
                hairpins_string += self._format_hairpin(hairpins[0], seq)
        #print the most stable hairpin
        if self.hairpinMin_3prim_dG() and self.hairpinMin_3prim_dG() < self._dG_threshold+3:
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