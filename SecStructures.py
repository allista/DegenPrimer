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

def hr(string, width=80):
        left_hr  = (width - len(string))/2
        right_hr = (width - len(string) - left_hr)
        return '-'*left_hr + string + '-'*right_hr + '\n\n'
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

    def __init__(self, seq_rec1, seq_rec2=None, seed_len=1):
        """Constructor"""
        self._seq_rec1 = seq_rec1
        self._seq1 = seq_rec1.seq
        if seq_rec2 and str(seq_rec1.seq) != str(seq_rec2.seq):
            self._seq_rec2 = seq_rec2
            self._seq2 = seq_rec2.seq
        else: 
            self._seq_rec2 = None
            self._seq2 = None
        self.setSeedLength(seed_len)
    #end def
    
    def isSingle(self):
        if self._seq2: return False
        return True
    #end def
    
    def setSeedLength(self, seed_len):
        if seed_len > len(self._seq1):
            seed_len = len(self._seq1)
        if self._seq2 and seed_len > len(self._seq2):
            seed_len = len(self._seq2)
        self._seed_length = seed_len
        self._recalculate()
    #end def

    def __str__(self):
        """For single sequence output self-dimers and hairpins.
        For two sequences output only cross-dimers"""
        string  = '\n'
        if self._seq2:
            string += hr(' %s vs %s cross-dimers ' % (self._seq_rec1.id, self._seq_rec2.id))
            string += self._format_dimers(self._cross_dimers, self._seq1, self._seq2)
        else:
            string += hr(' %s: %s ' % (self._seq_rec1.id, str(self._seq1)))
            string += hr(' self-dimers ')
            string += self._format_dimers(self._seq1_dimers, self._seq1, self._seq1)
            string += hr(' hairpins ')
            string += self._format_hairpins(self._seq1_hairpins, self._seq1)
        string += '-'*80+'\n\n'
        return string
    #end def
    
    def formatFull(self):
        """For each sequence output self-dimers and hairpins, then output cross-dimers"""
        string  = '\n'
        string += hr(' %s: %s ' % (self._seq_rec1.id, str(self._seq1)))
        string += hr(' self-dimers ')
        string += self._format_dimers(self._seq1_dimers, self._seq1, self._seq1)
        string += hr(' hairpins ')
        string += self._format_hairpins(self._seq1_hairpins, self._seq1)
        if self._seq2:
            string += '\n'
            string += hr(' %s: %s ' % (self._seq_rec2.id, str(self._seq2)))
            string += hr(' self-dimers ')
            string += self._format_dimers(self._seq2_dimers, self._seq2, self._seq2)
            string += hr(' hairpins ')
            string += self._format_hairpins(self._seq2_hairpins, self._seq2)
            string += '\n'
            string += hr(' %s vs %s cross-dimers ' % (self._seq_rec1.id, self._seq_rec2.id))
            string += self._format_dimers(self._cross_dimers, self._seq1, self._seq2)
        string += '-'*80+'\n\n'
        return string
    #end def
    
    def formatShort(self):
        """For single sequence output COUNTS and min dG of self-dimers and hairpins.
        For two sequences output the same only for cross-dimers"""
        string  = '\n'
        if self._seq2:
            string += hr(' %s vs %s cross-dimers ' % (self._seq_rec1.id, self._seq_rec2.id))
            string += self._format_dimers_short(self._cross_dimers, self._seq1, self._seq2)
        else:
            string += hr(' %s: %s ' % (self._seq_rec1.id, str(self._seq1)))
            string += hr(' self-dimers ')
            string += self._format_dimers_short(self._seq1_dimers, self._seq1, self._seq1)
            string += hr(' hairpins ')
            string += self._format_hairpins_short(self._seq1_hairpins, self._seq1)
        string += '-'*80+'\n\n'
        return string
    #end def
    
    def _recalculate(self):
        """For single sequence calculate self-dimers and hairpins.
        For two sequences calculate only cross-dimers"""
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
    
    def recalculateFull(self):
        """For each sequence calculate self-dimers and hairpins, then calculate cross-dimers"""
        all_structures       = self._find_dimers(self._seq1, self._seq1)
        self._seq1_dimers    = all_structures['dimers']
        self._seq1_hairpins  = self._find_hairpins(all_structures['hairpins'], self._seq1) 
        if self._seq2:
            all_structures       = self._find_dimers(self._seq2, self._seq2)
            self._seq2_dimers    = all_structures['dimers']
            self._seq2_hairpins  = self._find_hairpins(all_structures['hairpins'], self._seq2)
            self._cross_dimers   = self._find_dimers(self._seq1, self._seq2)['dimers']
        else:
            self._seq2_dimers    = None
            self._seq2_hairpins  = None
            self._cross_dimers   = None
    #end def
    
    def _find_matches(self, l, seq):
        matches = []
        position = seq.find(l)
        last = 0
        while position != -1:
            matches.append(last + position)
            last += position + 1
            position = seq[last:].find(l)
        return matches
    #end def

    def _find_seeds(self, seq_str1, seq_str2):
        """Search for all occurrences of all seq_str1 substrings of seed_length in seq_str2"""
        seeds = []
        scanned_patches = []
        for i in range(len(seq_str1) - self._seed_length):
            patch = seq_str1[i:i + self._seed_length]
            if patch not in scanned_patches:
                scanned_patches.append(patch)
                matches = self._find_matches(patch, seq_str2)
                if matches:
                    for position in matches:
                        seeds.append((i, position))
        return seeds
    #end def
    
    def _remove_loops(self, structure):
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
    
    def _find_dimers_by_slide(self, seq1, seq2):
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
        #sort by dG
        dimers.sort(key=lambda x: x[2][0])
        return all_structures
    #end def

    def _find_dimers_by_seed(self, seq1, seq2):
        """Search for all secondary structures containing at least one complementary 
        patch of given length"""
        structures, unique_structures, = [], []
        seq_str1 = str(seq1)
        seq_str2 = str(seq2.reverse_complement())
        seeds = self._find_seeds(seq_str1, seq_str2)
        if not seeds: return unique_structures
        for seed in seeds:
            pos1, pos2 = seed
            matches = (set(range(pos1, pos1+self._seed_length)), 
                       set(range(pos2, pos2+self._seed_length)),
                       [])
            #moving upstream
            up_pos1 = pos1-1
            up_pos2 = pos2-1
            while up_pos1 >= 0 and up_pos2 >= 0:
                if seq_str1[up_pos1] == seq_str2[up_pos2]:
                    matches[0].add(up_pos1)
                    matches[1].add(up_pos2)
                up_pos1 -= 1
                up_pos2 -= 1
            #moving downstream
            ds_pos1 = pos1+self._seed_length
            ds_pos2 = pos2+self._seed_length
            while ds_pos1 < len(seq_str1) and ds_pos2 < len(seq_str2):
                if seq_str1[ds_pos1] == seq_str2[ds_pos2]:
                    matches[0].add(ds_pos1)
                    matches[1].add(ds_pos2)
                ds_pos1 += 1
                ds_pos2 += 1
            #append structure
            structures.append(matches)
        #report only unique structures
        for struc in structures:
            if struc not in unique_structures:
                unique_structures.append(struc)
        #calculate dG
        for struc in unique_structures:
            struc[2].append(dimer_dG(struc, seq1, seq2))
        #sort structures by dG
        unique_structures.sort(key=lambda x: x[2][0])
        return unique_structures
    #end def
    
    def _find_dimers(self, seq1, seq2):
        if self._seed_length > 1:
            return self._find_dimers_by_seed(seq1, seq2)
        else: return self._find_dimers_by_slide(seq1, seq2)

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
        #count stable dimers
        num_stable  = 0
        three_prim  = 0
        tp_stable   = 0
        for struct in structures: 
            if struct[2][0] < 0:
                num_stable += 1
                #find 3' dimers
                fwd_matches = list(struct[0])
                rev_matches = list(struct[1])
                fwd_matches.sort()
                rev_matches.sort()
                if fwd_matches[-1] > len(seq1)-4 or \
                   rev_matches[0]  < 3:
                    three_prim += 1
                    if tp_stable > struct[2][0]:
                        tp_stable = struct[2][0]
        if not num_stable: return 'No dimers found.\n\n'
        #dimers count
        header_string += 'Dimers count: %d\n' % num_stable
        #3'-count
        if three_prim: header_string += '3\'-dimers count: %d\n' % three_prim
        #minimum dG
        header_string += 'Most stable: %.2f kcal/mol\n' % structures[0][2][0]
        if three_prim: header_string += 'Most stable 3\' dimer: %.2f kcal/mol\n' % tp_stable
        header_string += '\n'
        return header_string
    #end def

    def _format_dimers_short(self, structures, seq1, seq2):
        """short plain text representation of dimers containing at least one complementary 
        patch spanning seed_length"""
        structures_string  = ''
        if not structures: return 'No dimers found.\n\n'
        #print header
        structures_string += self._format_dimers_header(structures, seq1, seq2)
        #if there are stable dimers, print out the most stable
        if structures[0][2][0] < 0:
            structures_string += self._format_dimer(structures[0], seq1, seq2)
        return structures_string
    #end def

    def _format_dimers(self, structures, seq1, seq2):
        """plain text representation of dimers containing at least one complementary 
        patch spanning seed_length"""
        structures_string  = ''
        if not structures: return 'No dimers found.\n\n'
        #print header
        structures_string += self._format_dimers_header(structures, seq1, seq2)
        #print dimers
        for struct in structures:
            #only print stable dimers
            #if struct[2][0] >= 0: continue
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
        three_prim = 0
        tp_stable  = 0
        for hairpin in hairpins:
            rev_matches = list(hairpin[1])
            rev_matches.sort()
            #check for 3' hairpin
            if rev_matches[-1] > len(seq)-4:
                three_prim += 1
                if tp_stable > hairpin[2][0]:
                    tp_stable = hairpin[2][0]
        #hairpins count
        header_string += 'Hairpins count: %d\n' % len(hairpins)
        if three_prim: header_string += '3\'-hairpins count: %d\n' % three_prim
        #minimum dG
        header_string += 'Most stable: %.2f kcal/mol\n' % hairpins[0][2][0]
        if three_prim: header_string += 'Most stable 3\' hairpin: %.2f kcal/mol\n' % tp_stable
        header_string += '\n'
        return header_string
    #end def

    def _format_hairpins_short(self, hairpins, seq):
        hairpins_string  = ''
        if not hairpins: return 'No hairpins found\n\n'
        #print header
        hairpins_string += self._format_hairpins_header(hairpins, seq)
        #print the most stable hairpin
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

if __name__ == '__main__':
    print all_combinations([(1,5),(2,6),(3,7),(4,8)])