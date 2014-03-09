# coding=utf-8

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
Created on Jun 24, 2012

@author: Allis Tauri <allista@gmail.com>

All calculations are based on:
﻿
1) SantaLucia, J., & Hicks, D. (2004). 
The thermodynamics of DNA structural motifs. Annual review of biophysics and 
biomolecular structure, 33, 415-40. doi:10.1146/annurev.biophys.32.110601.141800

2) ﻿von Ahsen, N., Wittwer, C. T., & Schütz, E. (2001). Oligonucleotide melting 
temperatures under PCR conditions: nearest-neighbor corrections for Mg(2+), 
deoxynucleotide triphosphate, and dimethyl sulfoxide concentrations with 
comparison to alternative empirical formulas. Clinical chemistry, 47(11), 1956-61.
'''

from math import sqrt, log, exp
from UnifiedNN import UnifiedNN
import StringTools
###############################################################################


class PCR_Parameters(object):
    def __init__(self):
        #standard PCR conditions
        self.Mg    = 1.5e-3 #M
        self.Na    = 50e-3  #M
        self.dNTP  = 0.1e-3 #M
        self.DNA   = 5e-9   #M; DNA template concentration
        self.DMSO  = 0.0    #%; v/v-percent concentration of DMSO
        self.PCR_T = 60.0   #C; temperature at which PCR is conducted
    #end def
    
    def set(self, params):
        for name in vars(self).keys():
            if name in params:
                setattr(self, name, params[name])
#end class    


#Unified Nearest Neighbour model initialization
NN = UnifiedNN()
if not NN:
    raise Exception('TD_Functions: Unable to initialize UnifiedNN.')

#PCR conditions
PCR_P = PCR_Parameters()


#pseudo-concentration of Na equivalent to current concentration of Mg2+
def C_Na_eq():
    ''''divalent cation correction:
    all concentrations should be in mM (Ahsen et al., 2001)'''
    return (PCR_P.Na*1e3 + 120*sqrt((PCR_P.Mg - PCR_P.dNTP)*1e3))*1e-3
#end def


def primer_template_Tr(sequence, concentration, conversion_degree):
    '''Calculate temperature for primer-template annealing equilibrium 
    with a given conversin_degree using two-state equilibrium model and 
    the Nearest Neighbor thermodynamic tables (SantaLucia & Hicks, 2004).
    Note, that two-state equilibrium model used here is based on assumption, that 
    primer sequence is not self-complementary.'''
    #value constraints
    if not( 0 < conversion_degree < 1):
        raise ValueError('TD_Functions.primer_template_Tr: equilibrium ratio should be in the (0;1) interval.')
    #definitions
    seq_str = str(sequence)
    rev_str = str(sequence.complement())
    dH, dS = 0, 0
    #concentrations
    P   = concentration
    D   = PCR_P.DNA
    DUP = conversion_degree*min(P,D)
    #equilibrium constant 
    K   = DUP/((P-DUP)*(D-DUP))
    #initial corrections
    dH += NN.pair_dPar_37('ini', 'ini', 'dH')
    dS += NN.pair_dPar_37('ini', 'ini', 'dS')
    #test for AT terminals
    if seq_str[0] == 'A' or seq_str[0] == 'T':
        dH += NN.pair_dPar_37('ter', 'ter', 'dH')
        dS += NN.pair_dPar_37('ter', 'ter', 'dS')
    if seq_str[-1] == 'A' or seq_str[-1] == 'T':
        dH += NN.pair_dPar_37('ter', 'ter', 'dH')
        dS += NN.pair_dPar_37('ter', 'ter', 'dS')
    #stacking interactions
    for n in xrange(len(seq_str)-1):
        pair    = seq_str[n:n+2]
        reverse = rev_str[n:n+2]
        dH += NN.pair_dPar_37(pair, reverse, 'dH')
        dS += NN.pair_dPar_37(pair, reverse, 'dS')
    #salt concentration correction
    dS = dS + NN.dS_Na_coefficient * len(seq_str) * log(C_Na_eq()) #C_Na qM
    #final temperature calculation
    return NN.K0 + dH * 1000/(dS - NN.R * log(K)) - NN.T_DMSP_coefficient * PCR_P.DMSO #DMSO correction from [2]
#end def

def primer_template_Tm(sequence, concentration): 
    return primer_template_Tr(sequence, concentration, 0.5)


def format_concentration(concentration):
    return StringTools.format_quantity(concentration, 'M') 
    
    
def format_PCR_conditions(primers, polymerase=None):
    conditions = [['C(Na)',   ]+format_concentration(PCR_P.Na).split(),
                  ['C(Mg)',   ]+format_concentration(PCR_P.Mg).split(),
                  ['C(dNTP)', ]+format_concentration(PCR_P.dNTP).split(),
                  ['C(DNA)',  ]+format_concentration(PCR_P.DNA).split(),]
    if primers:
        for primer in primers:
            conditions.append(['C(%s)' % primer.id,]+format_concentration(primer.total_concentration).split())
    conditions.append(['C(DMSO)', '%.1f' % PCR_P.DMSO, '%'])
    if polymerase:
        conditions.append(['C(Poly)', ]+StringTools.format_quantity(polymerase*1e-6, 'u/ul').split())
    conditions.append(['T', '%.1f' % PCR_P.PCR_T, 'C'])
    return StringTools.print_table(conditions, delimiter='')
#end_def


def dimer_dG_corrected(dimer, seq1, seq2):
    '''calculate dG of a dimer corrected to current PCR conditions including 
    salt concentrations and temperature.
    dimer -- an instance of Dimer class which represents a dimer structure.
    seq1, seq2 -- sequences constituting a dimer, given in 5'->3' orientation.'''
    fwd_matches = dimer.fwd_matches
    #e.g. 5'-(2 ,3 ,4 ,8 ,9 )-3'
    rev_matches = dimer.rev_matches
    #e.g. 3-'(13,14,15,19,20)-5'
    seq_str = str(seq1)
    seq_len = len(seq_str)
    rev_str = str(seq2)[::-1]
    rev_len = len(rev_str)
    #initial dG of annealing
    dG = 0
    dH = NN.pair_dPar_37('ini', 'ini', 'dH')
    dS = NN.pair_dPar_37('ini', 'ini', 'dS')
    #check for 'left' dangling end
    if   fwd_matches[0] == 0 and rev_matches[0] > 0: #3' dangling
        pair    = rev_str[rev_matches[0]]+rev_str[rev_matches[0]-1]
        reverse = seq_str[0]+'-' 
        dH +=  NN.pair_dPar_37(pair, reverse, 'dH')
        dS +=  NN.pair_dPar_37(pair, reverse, 'dS')
    elif rev_matches[0] == 0 and fwd_matches[0] > 0: #5' dangling
        pair    = seq_str[fwd_matches[0]-1]+seq_str[fwd_matches[0]]
        reverse = '-'+rev_str[0] 
        dH +=  NN.pair_dPar_37(pair, reverse, 'dH')
        dS +=  NN.pair_dPar_37(pair, reverse, 'dS')
    #check for 'left' terminal mismatch
    elif fwd_matches[0] > 0 and rev_matches[0] > 0:
        dG += NN.Terminal_mismatch_mean
    #check for 'left' terminal AT
    elif fwd_matches[0] == 0 and rev_matches[0] == 0:
        if seq_str[0] == 'A' or seq_str[0] == 'T':
            dH += NN.pair_dPar_37('ter', 'ter', 'dH')
            dS += NN.pair_dPar_37('ter', 'ter', 'dS')
    #check for 'right' dangling end
    if   fwd_matches[-1] == seq_len-1 and rev_matches[-1] < rev_len-1: #5' dangling
        pair    = rev_str[rev_matches[-1]+1]+rev_str[rev_matches[-1]]
        reverse = '-'+seq_str[-1] 
        dH +=  NN.pair_dPar_37(pair, reverse, 'dH')
        dS +=  NN.pair_dPar_37(pair, reverse, 'dS')
    elif rev_matches[-1] == rev_len-1 and fwd_matches[-1] < seq_len-1: #3' dangling
        pair    = seq_str[fwd_matches[-1]]+seq_str[fwd_matches[-1]+1]
        reverse = rev_str[-1]+'-' 
        dH +=  NN.pair_dPar_37(pair, reverse, 'dH')
        dS +=  NN.pair_dPar_37(pair, reverse, 'dS')
    #check for 'right' terminal mismatch
    elif fwd_matches[-1]  < seq_len-1 and rev_matches[0] < rev_len-1:
        dG += NN.Terminal_mismatch_mean
    #check for 'right' terminal AT
    elif fwd_matches[-1] == seq_len-1 and rev_matches[-1] == rev_len-1:
        if seq_str[-1] == 'A' or seq_str[-1] == 'T':
            dH += NN.pair_dPar_37('ter', 'ter', 'dH')
            dS += NN.pair_dPar_37('ter', 'ter', 'dS')
    #stacking and mismatches
    for i in xrange(len(fwd_matches)-1):
        f_match = fwd_matches[i]
        f_next  = fwd_matches[i+1]
        r_match = rev_matches[i]
        r_next  = rev_matches[i+1]
        #if either || or |x| or |xx|
        if f_next-f_match < 4:
            pair    = seq_str[f_match:f_match+2]
            reverse = rev_str[r_match:r_match+2]
            #dH and salt-corrected dS
            dH +=  NN.pair_dPar_37(pair, reverse, 'dH')
            dS +=  NN.pair_dPar_37(pair, reverse, 'dS')
            #if ||
            if f_next-f_match == 1: continue
            #if |x| or |xx|
            elif f_next-f_match < 4:
                pair1    = rev_str[r_next-1:r_next+1][::-1]
                reverse1 = seq_str[f_next-1:f_next+1][::-1]
                dH +=  NN.pair_dPar_37(pair1, reverse1, 'dH')
                dS +=  NN.pair_dPar_37(pair1, reverse1, 'dS')
                continue
        #loop
        else:
            loop_seq = seq_str[f_match:f_next+1]
            dS += NN.loop_dS_37(loop_seq, 'internal')
    #dS salt correction
    dS += NN.dS_Na_coefficient * (fwd_matches[-1]-fwd_matches[0]) * log(C_Na_eq())
    return dG + dH - NN.temp_K(PCR_P.PCR_T)*dS/1000.0
#end def


def hairpin_dG_corrected(hairpin, seq):
    '''calculate dG of a hairpin corrected to current PCR conditions including 
    salt concentrations and temperature.
    hairpin -- an instance of Hairpin class which represents a hairpin structure.
    seq -- sequence which folds into the hairpin, given in 5'->3' orientation.'''
    fwd_matches = hairpin.fwd_matches
    #e.g. 5'-(2 ,3 ,4 ,8 ,9 )...-3'
    rev_matches = hairpin.rev_matches
    #e.g  3'-(24,23,22,18,17)...-5'
    seq_str = str(seq)
    seq_len = len(seq_str)
    #initial dG of annealing
    dG = 0
    dH = NN.pair_dPar_37('ini', 'ini', 'dH')
    dS = NN.pair_dPar_37('ini', 'ini', 'dS')
    #check for 'left' dangling end
    if   fwd_matches[0] == 0 and rev_matches[0] < seq_len-1: #3'-end
        pair    = seq_str[rev_matches[0]]+seq_str[rev_matches[0]+1]
        reverse = seq_str[0] + '-'
        dH +=  NN.pair_dPar_37(pair, reverse, 'dH')
        dS +=  NN.pair_dPar_37(pair, reverse, 'dS')
    elif fwd_matches[0] > 0 and rev_matches[0] == seq_len-1: #5'-end
        pair    = seq_str[fwd_matches[0]-1]+seq_str[fwd_matches[0]]
        reverse = '-' + seq_str[-1]
        dH +=  NN.pair_dPar_37(pair, reverse, 'dH')
        dS +=  NN.pair_dPar_37(pair, reverse, 'dS')
    #check for 'left' terminal mismatch
    elif fwd_matches[0] > 0 and rev_matches[0] < seq_len-1:
        dG += NN.Terminal_mismatch_mean
    #check for 'left' terminal AT
    elif fwd_matches[0] == 0 and rev_matches[0] == seq_len-1:
        if seq_str[0] == 'A' or seq_str[0] == 'T':
            dH += NN.pair_dPar_37('ter', 'ter', 'dH')
            dS += NN.pair_dPar_37('ter', 'ter', 'dS')
    #stacking and mismatches
    for i in xrange(len(fwd_matches)-1):
        f_match = fwd_matches[i]
        f_next  = fwd_matches[i+1]
        r_match = rev_matches[i]
        r_next  = rev_matches[i+1]
        #if either || or |x| or |xx|
        if f_next-f_match < 4:
            pair    = seq_str[f_match:f_match+2]
            reverse = seq_str[r_match-1:r_match+1][::-1]
            #dH and salt-corrected dS
            dH +=  NN.pair_dPar_37(pair, reverse, 'dH')
            dS +=  NN.pair_dPar_37(pair, reverse, 'dS')
            #if ||
            if f_next-f_match == 1: continue
            #if |x| or |xx|
            elif f_next-f_match < 4:
                pair1    = seq_str[r_next:r_next+2]
                reverse1 = seq_str[f_next-1:f_next+1][::-1]
                dH +=  NN.pair_dPar_37(pair1, reverse1, 'dH')
                dS +=  NN.pair_dPar_37(pair1, reverse1, 'dS')
                continue
        #loop
        else:
            loop_seq = seq_str[f_match:f_next+1]
            dS += NN.loop_dS_37(loop_seq, 'internal')
    #hairpin loop
    hp_str = seq_str[fwd_matches[-1]:rev_matches[-1]+1]
    dH += NN.loop_dH_37(hp_str, 'hairpin')
    dS += NN.loop_dS_37(hp_str, 'hairpin')
    #dS salt correction
    dS += NN.dS_Na_coefficient * (fwd_matches[-1]-fwd_matches[0]) * log(C_Na_eq())
    return dG + dH - NN.temp_K(PCR_P.PCR_T)*dS/1000.0
#end def


def equilibrium_constant(dG_T, T):
    '''calculate equilibrium constant of annealing at a given temperature, 
    given standard dG at this temperature'''
    return exp(-1000*dG_T/(NN.R*NN.temp_K(T))) #annealing equilibrium constant
#end def


def primer_DNA_conversion_degree(primer_concentration, K):
    '''calculate conversion degree of Primer/DNA dimerisation
    given primer_concentration and equilibrium constant'''
    P = primer_concentration
    D = PCR_P.DNA
    #quadratic equation with respect to DUP = r*min(P,D), 
    #where 'r' is a conversion degree
    _b   = (K*P+K*D+1) #MINUS b; always positive
    disc = _b*_b - 4*K*(K*P*D) #this should always be >= 0 given non-negative K, P and D
    DUP  = (_b-sqrt(disc))/(2*K) #take the smallest positive root
    return DUP/min(P,D)
#end def


#tests
if __name__ == '__main__':
    from Bio.Seq import Seq
    from SecStructures import Dimer
    seq1 = Seq('ATATTCTACGACGGCTATCC')
    seq2 = Seq('ATATTCTACGACGGCTATCC').reverse_complement()
    seq3 = Seq('GAACGCAAAGATTGGGAAC')
    seq4 = Seq('GAACGCAAAGATTCTGAAC').reverse_complement()
    dimer12 = Dimer.from_sequences(seq1, seq2)
    dimer34 = Dimer.from_sequences(seq3, seq4)
    
    PCR_P.Na     = 50.0e-3
    PCR_P.Mg     = 3.0e-3 
    PCR_P.dNTP   = 0.1e-3
    PCR_P.DNA    = 1.0e-9
    PCR_P.Primer = 0.43e-6
    PCR_P.PCR_T    = 65
    
    print vars(PCR_P)
    PCR_P.set({'PCR_T': 30})
    print vars(PCR_P)
    
    print seq1, seq2
    print dimer_dG_corrected(dimer12, seq1, seq2)
    print ''
    print seq3, seq4
    print dimer_dG_corrected(dimer34, seq3, seq4)