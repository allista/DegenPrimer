# Copyright (C) 2012 Allis Tauri <allista@gmail.com>
# 
# indicator_gddccontrol is free software: you can redistribute it and/or modify it
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
from copy import deepcopy
from TD_Functions import dimer_dG, hairpin_dG

def hr(string, width=80, symbol='-'):
        left_hr  = (width - len(string))/2
        right_hr = (width - len(string) - left_hr)
        return symbol*left_hr + string + symbol*right_hr + '\n\n'
#end def
    
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
            return self._cross_dimers[0][2][0]
        elif self._seq1_dimers:
            return self._seq1_dimers[0][2][0]
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
        else: return self._seq1_hairpins[0][2][0]
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
    
    def _remove_loops(self, structure):
        """remove terminal internal loops one by one in all possible combinations
        to produce a bunch of new structures"""
        structures  = [structure]
        fwd_matches = list(structure[0])
        rev_matches = list(structure[1])
        fwd_matches.sort()
        rev_matches.sort()
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
        if not left_loops: return structures
        #else, remove boundary loops one by one
        l_structure = structure
        for l in range(len(left_loops)):
            r_structure = l_structure
            for r in range(len(right_loops)-l):
                r_structure = (r_structure[0]-right_loops[r][0][0],
                               r_structure[1]-right_loops[r][0][1],[])
                structures.append(r_structure)
            l_structure = (l_structure[0]-left_loops[l][0][0],
                           l_structure[1]-left_loops[l][0][1],[])
            structures.append(l_structure)
        return structures
    #end def
    
    def _find_dimers(self, seq1, seq2):
        """Search for all secondary structures"""
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
            structure = (set(), set(), [])
            for fwd,rev in zip(range(fwd_left,fwd_right),range(rev_left,rev_right)):
                if(fwd_str[fwd] == rev_str[rev]):
                    structure[0].add(fwd)
                    structure[1].add(rev)
            #if found some
            if structure[0]:
                if self_dimers and structure not in hairpins: 
                    hairpins.append(structure)
                min_dG = 0
                min_dG_structure = None
                new_structures = self._remove_loops(structure)
                for struct in new_structures:
                    #calculate dG for every structure
                    dG = dimer_dG(struct, seq1, seq2)
                    struct[2].append(dG)
                    #find the most stable
                    if min_dG > dG:
                        min_dG = dG
                        min_dG_structure = struct
                #if there is a stable structure, append it
                if min_dG_structure and min_dG_structure not in dimers:
                    dimers.append(min_dG_structure)
                    #check for 3' structure
                    if max(min_dG_structure[0]) > len(seq1)-4 \
                    or min(min_dG_structure[1]) < 3:
                        self._3prim_dimers += 1
                        if self._min_3prim_dimer_dG > min_dG:
                            self._min_3prim_dimer_dG = min_dG
                            self._min_3prim_dimer    = min_dG_structure 
        #sort by dG
        dimers.sort(key=lambda x: x[2][0])
        return all_structures
    #end def

    def _find_hairpins(self, dimers, seq):
        """Find stable hairpins in given set of SYMMETRIC SELF-dimers"""
        seq_len = len(seq)
        hairpins   = []
        if not dimers: return hairpins
        for dimer in dimers:
            fwd_matches = deepcopy(dimer[0])
            rev_matches = deepcopy(dimer[0])
            while (rev_matches and fwd_matches)       and \
                  (min(rev_matches) < max(fwd_matches) or \
                   min(rev_matches) - max(fwd_matches) < 4): #at least 3 bases separation
                fwd_matches.remove(max(fwd_matches))
                rev_matches.remove(min(rev_matches))
            if fwd_matches:
                min_dG = 0
                min_dG_hairpin = None
                new_dimers = self._remove_loops((fwd_matches, set(seq_len-1-m for m in rev_matches), []))
                for new_dimer in new_dimers:
                    hairpin = (deepcopy(new_dimer[0]), set(seq_len-1-m for m in new_dimer[1]), [])
                    dG = hairpin_dG(hairpin, seq)
                    hairpin[2].append(dG)
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
        hairpins.sort(key=lambda x: x[2][0])
        return hairpins
    #end def

    def _format_dimer(self, structure, seq1, seq2):
        structure_string = ''
        fwd_matches = list(structure[0])
        rev_matches = list(structure[1])
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
        structure_string += 'dG(37C) = %.2f kcal/mol\n' % structure[2][0]
        #check for 3' dimer
        if fwd_matches[-1] > len(seq1)-4 or \
           rev_matches[0]  < 3:
            structure_string += '(3\'-dimer)\n'
        structure_string += '\n'
        return structure_string
    #end def

    def _format_dimers_header(self, structures, seq1, seq2):
        header_string  = ''
        if not structures: return 'No dimers found.\n\n'
        #dimers count
        header_string += 'Dimers count: %d\n' % len(structures)
        #3'-count
        if self._3prim_dimers: 
            header_string += '3\'-dimers count: %d\n' % self._3prim_dimers
        #minimum dG
        header_string += 'Most stable: %.2f kcal/mol\n' % structures[0][2][0]
        if self._3prim_dimers: 
            header_string += 'Most stable 3\' dimer: %.2f kcal/mol\n' % self._min_3prim_dimer_dG
        header_string += '\n'
        return header_string
    #end def

    def _format_dimers_short(self, structures, seq1, seq2):
        """short plain text representation of dimers"""
        structures_string  = ''
        if not structures: return 'No dimers found.\n\n'
        #print header
        structures_string += self._format_dimers_header(structures, seq1, seq2)
        #print the most stable dimer
        if self.dimerMin_dG() and self.dimerMin_dG() < self._dG_threshold:
            if structures[0] != self._min_3prim_dimer: #if it is the 3' the next block will print it
                structures_string += self._format_dimer(structures[0], seq1, seq2)
        #print the most stable 3'-dimer
        if self.dimerMin_3prim_dG() and self.dimerMin_3prim_dG() < self._dG_threshold+1:
            structures_string += self._format_dimer(self._min_3prim_dimer, seq1, seq2)
        return structures_string
    #end def

    def _format_dimers(self, structures, seq1, seq2):
        """plain text representation of dimers"""
        structures_string  = ''
        if not structures: return 'No dimers found.\n\n'
        #print header
        structures_string += self._format_dimers_header(structures, seq1, seq2)
        #print dimers
        for struct in structures:
            structures_string += self._format_dimer(struct, seq1, seq2)
        return structures_string
    #end def

    def _format_hairpin(self, hairpin, seq):
        hairpin_string = ''
        fwd_matches = list(hairpin[0])
        rev_matches = list(hairpin[1])
        fwd_matches.sort()
        rev_matches.sort()
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
        hairpin_string += 'dG(37C) = %.2f kcal/mol\n' % hairpin[2][0]
        #check for 3' hairpin
        if rev_matches[-1] > len(seq)-4:
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
        header_string += 'Most stable: %.2f kcal/mol\n' % hairpins[0][2][0]
        if self._3prim_hairpins: 
            header_string += 'Most stable 3\' hairpin: %.2f kcal/mol\n' % self._min_3prim_hairpin_dG
        header_string += '\n'
        return header_string
    #end def

    def _format_hairpins_short(self, hairpins, seq):
        hairpins_string  = ''
        if not hairpins: return 'No hairpins found\n\n'
        #print header
        if self.hairpinMin_dG() and self.hairpinMin_dG() < self._dG_threshold+2:
            if hairpins[0] != self._min_3prim_hairpin: #if it is the 3' the next block will print it
                hairpins_string += self._format_hairpins_header(hairpins, seq)
        #print the most stable hairpin
        if self.hairpinMin_3prim_dG() and self.hairpinMin_3prim_dG() < self._dG_threshold+3:
            hairpins_string += self._format_hairpin(hairpins[0], seq)
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