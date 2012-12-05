# coding=utf-8

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
from StringTools import print_exception
try:
    from Bio.SeqFeature import SeqFeature, FeatureLocation
except ImportError, e:
    print 'The BioPython must be installed in your system.'
    raise
###############################################################################

#standard PCR conditions
C_Mg   = 1.5  #mM
C_Na   = 50   #mM; should be above 0.05M and below 1.1M
C_dNTP = 0    #mM
C_DNA  = 50   #nM; DNA template concentration
C_Prim = 0.1  #uM; Primer concentration
C_DMSO = 0    #percent
PCR_T  = 37   #C;  temperature at which PCR is conducted. Default is 37 for dG tables contain values of standard Gibbs energy at 37C

#Unified Nearest Neighbour model initialization
NN = UnifiedNN()
if not NN:
    raise Exception('TD_Functions: Unable to initialize UnifiedNN.')

#pseudo-concentration of Na equivalent to current concentration of Mg2+
def C_Na_eq():
    """divalent cation correction (Ahsen et al., 2001)"""
    return C_Na + 120*sqrt(C_Mg - C_dNTP)
#end def

#concentrations of primer and DNA in mols
def C_Prim_M(): return C_Prim*1e-6
def C_DNA_M(): return C_DNA*1e-9


def NN_Tr(seq, r):
    '''Calculate temperature for primer-template association equilibrium 
    with 'r' ratio using two-state equilibrium model and the Nearest Neighbor \
    TD tables and from the paper of SantaLucia & Hicks (2004).
    Note, that two-state equilibrium model used here is based on assumption, that 
    primer sequence is not self-complementary.'''
    #value constraints
    if r >=1 or r <=0:
        raise ValueError('TD_Functions.NN_Tr: equilibrium ratio should be in the (0;1) interval.')
    #definitions
    seq_str = str(seq)
    rev_str = str(seq.complement())
    dH, dS = 0, 0
    #concentrations
    P   = C_Prim_M()
    D   = C_DNA_M()
    DUP = r*min(P,D)
    #equilibrium constant 
    K   = DUP/((P-DUP)*(D-DUP))
    #initial corrections
    dH += NN.pair_dH_37('ini', 'ini')
    dS += NN.pair_dS_37('ini', 'ini')
    #test for AT terminals
    if seq_str[0] == 'A' or seq_str[0] == 'T':
        dH += NN.pair_dH_37('ter', 'ter')
        dS += NN.pair_dS_37('ter', 'ter')
    if seq_str[-1] == 'A' or seq_str[-1] == 'T':
        dH += NN.pair_dH_37('ter', 'ter')
        dS += NN.pair_dS_37('ter', 'ter')
    #stacking interactions
    for n in range(len(seq_str)-1):
        pair    = seq_str[n:n+2]
        reverse = rev_str[n:n+2]
        dH += NN.pair_dH_37(pair, reverse)
        dS += NN.pair_dS_37(pair, reverse)
    #salt concentration correction
    dS = dS + NN.dS_Na_coefficient * len(seq_str) * log(C_Na_eq()*1e-3) #C_Na mM
    #final temperature calculation
    return NN.K0 + dH * 1000/(dS - NN.R * log(K)) - NN.T_DMSP_coefficient * C_DMSO #DMSO correction from [2]
#end def

def NN_Tm(seq): return NN_Tr(seq, 0.5)


def source_feature(seq_rec):
    feature = None
    for f in seq_rec.features:
        if f.type == 'source':
            feature = f
            break
    if not feature:
        feature = SeqFeature(FeatureLocation(0,len(seq_rec.seq)),
                             type = 'source')
        seq_rec.features.append(feature)
    return feature
#end def


def format_PCR_conditions():
    conc_str  = ''
    spacer = max(len(str(C_Na)), 
                 len(str(C_Mg)), 
                 len(str(C_dNTP)),
                 len(str(C_DNA)),
                 len(str(C_Prim)),
                 len(str(C_DMSO)),
                 len(str(PCR_T)))
    conc_str += 'C(Na)     = ' + str(C_Na)  + ' '*(spacer-len(str(C_Na)))   +' mM\n'
    conc_str += 'C(Mg)     = ' + str(C_Mg)  + ' '*(spacer-len(str(C_Mg)))   +' mM\n'
    conc_str += 'C(dNTP)   = ' + str(C_dNTP)+ ' '*(spacer-len(str(C_dNTP))) +' mM\n'
    conc_str += 'C(DNA)    = ' + str(C_DNA) + ' '*(spacer-len(str(C_DNA)))  +' nM\n'
    conc_str += 'C(Primer) = ' + str(C_Prim)+ ' '*(spacer-len(str(C_Prim))) +' uM\n'
    conc_str += 'C(DMSO)   = ' + str(C_DMSO)+ ' '*(spacer-len(str(C_DMSO))) +' %\n'
    conc_str += 'T         = ' + str(PCR_T) + ' '*(spacer-len(str(PCR_T)))  +' C\n'
    return conc_str
#end_def


def add_PCR_conditions(feature):
    try:
        feature.qualifiers['C_Na']      = str(C_Na)+  ' mM'
        feature.qualifiers['C_Mg']      = str(C_Mg)+  ' mM'
        feature.qualifiers['C_dNTP']    = str(C_dNTP)+' mM'
        feature.qualifiers['C_DNA']     = str(C_DNA)+ ' nM'
        feature.qualifiers['C_Primer']  = str(C_Prim)+' uM'
        feature.qualifiers['C_DMSO']    = str(C_DMSO)+' %'
        feature.qualifiers['T']         = str(PCR_T)+' C'
    except Exception, e:
        print 'add_PCR_conditions:'
        print_exception(e)
#end def


def calculate_Tr(seq_rec, r):
    primer_Tr = NN_Tr(seq_rec.seq, r)
    feature = source_feature(seq_rec)
    add_PCR_conditions(feature)
    feature.qualifiers['T-'+str(r)] = str(primer_Tr)
    return primer_Tr
#end def


def calculate_Tm(seq_rec):
    primer_Tm = NN_Tm(seq_rec.seq)
    feature = source_feature(seq_rec)
    add_PCR_conditions(feature)
    feature.qualifiers['Tm'] = str(primer_Tm)
    return primer_Tm
#end def


def dimer_dG(dimer, seq1, seq2):
    '''Calculate 'standard' dG at 37C of dimer annealing process.
    dimer -- an instance of Dimer class which represents a dimer structure.
    seq1, seq2 -- sequences constituting a dimer, given in 5'->3' orientation.'''
    fwd_matches = dimer.fwd_matches()
    #e.g. 5'-(2 ,3 ,4 ,8 ,9 )-3'
    rev_matches = dimer.rev_matches()
    #e.g. 3-'(13,14,15,19,20)-5'
    seq_str = str(seq1)
    seq_len = len(seq_str)
    rev_str = str(seq2[::-1])
    rev_len = len(rev_str)
    #initial dG of annealing
    dG = NN.pair_dG_37('ini', 'ini')
    #check for 'left' dangling end
    if   fwd_matches[0] == 0 and rev_matches[0] > 0: #3' dangling
        pair    = rev_str[rev_matches[0]]+rev_str[rev_matches[0]-1]
        reverse = seq_str[0]+'-' 
        dG +=  NN.pair_dG_37(pair, reverse)
    elif rev_matches[0] == 0 and fwd_matches[0] > 0: #5' dangling
        pair    = seq_str[fwd_matches[0]-1]+seq_str[fwd_matches[0]]
        reverse = '-'+rev_str[0] 
        dG +=  NN.pair_dG_37(pair, reverse)
    #check for 'left' terminal mismatch
    elif fwd_matches[0] > 0 and rev_matches[0] > 0:
        dG += NN.Terminal_mismatch_mean
    #check for 'left' terminal AT
    elif fwd_matches[0] == 0 and rev_matches[0] == 0:
        if seq_str[0] == 'A' or seq_str[0] == 'T':
            dG += NN.pair_dG_37('ter', 'ter')
    #check for 'right' dangling end
    if   fwd_matches[-1] == seq_len-1 and rev_matches[-1] < rev_len-1: #5' dangling
        pair    = rev_str[rev_matches[-1]+1]+rev_str[rev_matches[-1]]
        reverse = '-'+seq_str[-1] 
        dG +=  NN.pair_dG_37(pair, reverse)
    elif rev_matches[-1] == rev_len-1 and fwd_matches[-1] < seq_len-1: #3' dangling
        pair    = seq_str[fwd_matches[-1]]+seq_str[fwd_matches[-1]+1]
        reverse = rev_str[-1]+'-' 
        dG +=  NN.pair_dG_37(pair, reverse)
    #check for 'right' terminal mismatch
    elif fwd_matches[-1]  < seq_len-1 and rev_matches[0] < rev_len-1:
        dG += NN.Terminal_mismatch_mean
    #check for 'right' terminal AT
    elif fwd_matches[-1] == seq_len-1 and rev_matches[-1] == rev_len-1:
        if seq_str[-1] == 'A' or seq_str[-1] == 'T':
            dG += NN.pair_dG_37('ter', 'ter')
    #stacking and mismatches
    for i in range(len(fwd_matches)-1):
        f_match = fwd_matches[i]
        f_next  = fwd_matches[i+1]
        r_match = rev_matches[i]
        r_next  = rev_matches[i+1]
        #if either || or |x| or |xx|
        if f_next-f_match < 4:
            pair    = seq_str[f_match:f_match+2]
            reverse = rev_str[r_match:r_match+2]
            #salt-corrected dG
            dG += NN.pair_dG_37(pair, reverse)
            #if ||
            if f_next-f_match == 1: continue
            #if |x| or |xx|
            elif f_next-f_match < 4:
                pair1    = rev_str[r_next-1:r_next+1][::-1]
                reverse1 = seq_str[f_next-1:f_next+1][::-1]
                dG += NN.pair_dG_37(pair1, reverse1)
                continue
        #loop
        else:
            dG += NN.internal_loop_dG_37(f_next-f_match-1)
    #dG salt correction
    dG_Na   = NN.dG_Na_coefficient_oligo * (fwd_matches[-1]-fwd_matches[0]) * log(C_Na_eq()*1e-3)
    return dG + dG_Na
#end def


def dimer_dG_corrected(dimer, seq1, seq2):
    '''calculate dG of a dimer corrected to current PCR conditions including 
    salt concentrations and temperature. Dummy for now'''
    fwd_matches = dimer.fwd_matches()
    #e.g. 5'-(2 ,3 ,4 ,8 ,9 )-3'
    rev_matches = dimer.rev_matches()
    #e.g. 3-'(13,14,15,19,20)-5'
    seq_str = str(seq1)
    seq_len = len(seq_str)
    rev_str = str(seq2[::-1])
    rev_len = len(rev_str)
    #initial dG of annealing
    dG = 0
    dH = NN.pair_dH_37('ini', 'ini')
    dS = NN.pair_dS_37('ini', 'ini')
    #check for 'left' dangling end
    if   fwd_matches[0] == 0 and rev_matches[0] > 0: #3' dangling
        pair    = rev_str[rev_matches[0]]+rev_str[rev_matches[0]-1]
        reverse = seq_str[0]+'-' 
        dH +=  NN.pair_dH_37(pair, reverse)
        dS +=  NN.pair_dS_37(pair, reverse)
    elif rev_matches[0] == 0 and fwd_matches[0] > 0: #5' dangling
        pair    = seq_str[fwd_matches[0]-1]+seq_str[fwd_matches[0]]
        reverse = '-'+rev_str[0] 
        dH +=  NN.pair_dH_37(pair, reverse)
        dS +=  NN.pair_dS_37(pair, reverse)
    #check for 'left' terminal mismatch
    elif fwd_matches[0] > 0 and rev_matches[0] > 0:
        dG += NN.Terminal_mismatch_mean
    #check for 'left' terminal AT
    elif fwd_matches[0] == 0 and rev_matches[0] == 0:
        if seq_str[0] == 'A' or seq_str[0] == 'T':
            dH += NN.pair_dH_37('ter', 'ter')
            dS += NN.pair_dS_37('ter', 'ter')
    #check for 'right' dangling end
    if   fwd_matches[-1] == seq_len-1 and rev_matches[-1] < rev_len-1: #5' dangling
        pair    = rev_str[rev_matches[-1]+1]+rev_str[rev_matches[-1]]
        reverse = '-'+seq_str[-1] 
        dH +=  NN.pair_dH_37(pair, reverse)
        dS +=  NN.pair_dS_37(pair, reverse)
    elif rev_matches[-1] == rev_len-1 and fwd_matches[-1] < seq_len-1: #3' dangling
        pair    = seq_str[fwd_matches[-1]]+seq_str[fwd_matches[-1]+1]
        reverse = rev_str[-1]+'-' 
        dH +=  NN.pair_dH_37(pair, reverse)
        dS +=  NN.pair_dS_37(pair, reverse)
    #check for 'right' terminal mismatch
    elif fwd_matches[-1]  < seq_len-1 and rev_matches[0] < rev_len-1:
        dG += NN.Terminal_mismatch_mean
    #check for 'right' terminal AT
    elif fwd_matches[-1] == seq_len-1 and rev_matches[-1] == rev_len-1:
        if seq_str[-1] == 'A' or seq_str[-1] == 'T':
            dH += NN.pair_dH_37('ter', 'ter')
            dS += NN.pair_dS_37('ter', 'ter')
    #stacking and mismatches
    for i in range(len(fwd_matches)-1):
        f_match = fwd_matches[i]
        f_next  = fwd_matches[i+1]
        r_match = rev_matches[i]
        r_next  = rev_matches[i+1]
        #if either || or |x| or |xx|
        if f_next-f_match < 4:
            pair    = seq_str[f_match:f_match+2]
            reverse = rev_str[r_match:r_match+2]
            #dH and salt-corrected dS
            dH +=  NN.pair_dH_37(pair, reverse)
            dS +=  NN.pair_dS_37(pair, reverse)
            #if ||
            if f_next-f_match == 1: continue
            #if |x| or |xx|
            elif f_next-f_match < 4:
                pair1    = rev_str[r_next-1:r_next+1][::-1]
                reverse1 = seq_str[f_next-1:f_next+1][::-1]
                dH +=  NN.pair_dH_37(pair1, reverse1)
                dS +=  NN.pair_dS_37(pair1, reverse1)
                continue
        #loop
        else:
            loop_seq = seq_str[f_match:f_next+1]
            dS += NN.loop_dS_37(loop_seq, 'internal')
    #dS salt correction
    dS += NN.dS_Na_coefficient * (fwd_matches[-1]-fwd_matches[0]) * log(C_Na_eq()*1e-3)
    return dG + dH - NN.temp_K(PCR_T)*dS/1000.0
#end def


def hairpin_dG(hairpin, seq):
    '''Calculate 'standard' dG at 37C of hairpin annealing process.
    hairpin -- an instance of Hairpin class which represents a hairpin structure.
    seq -- sequence which folds into the hairpin, given in 5'->3' orientation.'''
    fwd_matches = hairpin.fwd_matches()
    #e.g. 5'-(2 ,3 ,4 ,8 ,9 )...-3'
    rev_matches = hairpin.rev_matches()
    #e.g  3'-(24,23,22,18,17)...-5'
    seq_str = str(seq)
    seq_len = len(seq_str)
    #initial dG of annealing
    dG = NN.pair_dG_37('ini', 'ini')
    #check for 'left' dangling end
    if   fwd_matches[0] == 0 and rev_matches[0] < seq_len-1: #3'-end
        pair    = seq_str[rev_matches[0]]+seq_str[rev_matches[0]+1]
        reverse = seq_str[0] + '-'
        dG +=  NN.pair_dG_37(pair, reverse)
    elif fwd_matches[0] > 0 and rev_matches[0] == seq_len-1: #5'-end
        pair    = seq_str[fwd_matches[0]-1]+seq_str[fwd_matches[0]]
        reverse = '-' + seq_str[-1]
        dG +=  NN.pair_dG_37(pair, reverse)
    #check for 'left' terminal mismatch
    elif fwd_matches[0] > 0 and rev_matches[0] < seq_len-1:
        dG += NN.Terminal_mismatch_mean
    #check for 'left' terminal AT
    elif fwd_matches[0] == 0 and rev_matches[0] == seq_len-1:
        if seq_str[0] == 'A' or seq_str[0] == 'T':
            dG += NN.pair_dG_37('ter', 'ter')
    #stacking and mismatches
    for i in range(len(fwd_matches)-1):
        f_match = fwd_matches[i]
        f_next  = fwd_matches[i+1]
        r_match = rev_matches[i]
        r_next  = rev_matches[i+1]
        #if either || or |x| or |xx|
        if f_next-f_match < 4:
            pair    = seq_str[f_match:f_match+2]
            reverse = seq_str[r_match-1:r_match+1][::-1]
            #salt-corrected dG
            dG += NN.pair_dG_37(pair, reverse)
            #if ||
            if f_next-f_match == 1: continue
            #if |x| or |xx|
            elif f_next-f_match < 4:
                pair1    = seq_str[r_next:r_next+2]
                reverse1 = seq_str[f_next-1:f_next+1][::-1]
                dG += NN.pair_dG_37(pair1, reverse1)
        #internal loop
        else:
            dG += NN.internal_loop_dG_37(f_next-f_match-1)
    #hairpin loop
    hp_str = seq_str[fwd_matches[-1]:rev_matches[-1]+1]
    dG += NN.hairpin_loop_dG_37(hp_str)
    #dG salt correction
    dG_Na   = NN.dG_Na_coefficient_oligo * (fwd_matches[-1]-fwd_matches[0]) * log(C_Na_eq()*1e-3)
    return dG + dG_Na
#end def


def hairpin_dG_corrected(hairpin, seq):
    '''calculate dG of a hairpin corrected to current PCR conditions including 
    salt concentrations and temperature.'''
    fwd_matches = hairpin.fwd_matches()
    #e.g. 5'-(2 ,3 ,4 ,8 ,9 )...-3'
    rev_matches = hairpin.rev_matches()
    #e.g  3'-(24,23,22,18,17)...-5'
    seq_str = str(seq)
    seq_len = len(seq_str)
    #initial dG of annealing
    dG = 0
    dH = NN.pair_dH_37('ini', 'ini')
    dS = NN.pair_dS_37('ini', 'ini')
    #check for 'left' dangling end
    if   fwd_matches[0] == 0 and rev_matches[0] < seq_len-1: #3'-end
        pair    = seq_str[rev_matches[0]]+seq_str[rev_matches[0]+1]
        reverse = seq_str[0] + '-'
        dH +=  NN.pair_dH_37(pair, reverse)
        dS +=  NN.pair_dS_37(pair, reverse)
    elif fwd_matches[0] > 0 and rev_matches[0] == seq_len-1: #5'-end
        pair    = seq_str[fwd_matches[0]-1]+seq_str[fwd_matches[0]]
        reverse = '-' + seq_str[-1]
        dH +=  NN.pair_dH_37(pair, reverse)
        dS +=  NN.pair_dS_37(pair, reverse)
    #check for 'left' terminal mismatch
    elif fwd_matches[0] > 0 and rev_matches[0] < seq_len-1:
        dG += NN.Terminal_mismatch_mean
    #check for 'left' terminal AT
    elif fwd_matches[0] == 0 and rev_matches[0] == seq_len-1:
        if seq_str[0] == 'A' or seq_str[0] == 'T':
            dH += NN.pair_dH_37('ter', 'ter')
            dS += NN.pair_dS_37('ter', 'ter')
    #stacking and mismatches
    for i in range(len(fwd_matches)-1):
        f_match = fwd_matches[i]
        f_next  = fwd_matches[i+1]
        r_match = rev_matches[i]
        r_next  = rev_matches[i+1]
        #if either || or |x| or |xx|
        if f_next-f_match < 4:
            pair    = seq_str[f_match:f_match+2]
            reverse = seq_str[r_match-1:r_match+1][::-1]
            #dH and salt-corrected dS
            dH +=  NN.pair_dH_37(pair, reverse)
            dS +=  NN.pair_dS_37(pair, reverse)
            #if ||
            if f_next-f_match == 1: continue
            #if |x| or |xx|
            elif f_next-f_match < 4:
                pair1    = seq_str[r_next:r_next+2]
                reverse1 = seq_str[f_next-1:f_next+1][::-1]
                dH +=  NN.pair_dH_37(pair1, reverse1)
                dS +=  NN.pair_dS_37(pair1, reverse1)
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
    dS += NN.dS_Na_coefficient * (fwd_matches[-1]-fwd_matches[0]) * log(C_Na_eq()*1e-3)
    return dG + dH - NN.temp_K(PCR_T)*dS/1000.0
#end def


def equilibrium_constant(dG_T, T):
    '''calculate equilibrium constant of at a given temperature, 
    given standard dG at this temperature'''
    return exp(-1000*dG_T/(NN.R*NN.temp_K(T))) #annealing equilibrium constant


def primer_DNA_conversion_degree(dG_T, T):
    '''calculate conversion degree of Primer/DNA dimerisation
    given standard dG(kcal/mol) of annealing at T(C) temperature'''
    K = equilibrium_constant(dG_T, T)
    P = C_Prim_M()
    D = C_DNA_M()
    #quadratic equation with respect to DUP = r*min(P,D), 
    #where 'r' is a conversion degree
    _b   = (K*P+K*D+1) #MINUS b; always positive
    disc = _b*_b - 4*K*(K*P*D) #this should always be >= 0 given non-negative K, P and D
    DUP  = (_b-sqrt(disc))/(2*K) #take the smallest positive root
    return DUP/min(P,D)
#end def


def dimer_conversion_degree(dG_T, T):
    '''calculate conversion degree of Primer self-dimerisation
    given standard dG(kcal/mol) of annealing at T(C) temperature'''
    K = equilibrium_constant(dG_T, T)
    P = C_Prim_M()
    #quadratic equation with respect to DUP = r*P, 
    #where 'r' is a conversion degree
    _b   = (2*K*P+1) #MINUS b; always positive
    disc = _b*_b - 4*K*(K*P*P) #this should always be >= 0 given non-negative K and P
    DUP  = (_b-sqrt(disc))/(2*K) #take the smallest positive root
    return DUP/P
#end def


def hairpin_conversion_degree(dG_T, T):
    '''calculate conversion degree of hairpin formation on Primer
    given standard dG(kcal/mol) of annealing at T(C) temperature'''
    K = equilibrium_constant(dG_T, T)
    return K/(1+K)
#end def


#tests
if __name__ == '__main__':
    from Bio.Seq import Seq
    from OligoFunctions import compose_dimer
    seq1 = Seq('ATATTCTACGACGGCTATCC').reverse_complement()
    seq2 = Seq('ATATTCTACAACGGCTATCC') #
    seq3 = Seq('ATATTCTACAACGGCTATCC')
    seq4 = Seq('CCTATCGGCAACATCTTATA'[::-1])
    dim1 = compose_dimer(seq1, seq2) #Dimer((0,1,2,3), (16,17,18,19))
    dim2 = compose_dimer(seq3, seq4)
    
    C_Na     = 50.0
    C_Mg     = 3.0 
    C_dNTP   = 0.3 
    C_DNA    = 50.0
    C_Primer = 0.43
    PCR_T = 55
    
    _dG = 3
    const = equilibrium_constant(_dG, PCR_T)
    print 'Equilibrium constant:     ', const
    print 'Hairpin conversion degree:', hairpin_conversion_degree(_dG, PCR_T)*100
    print 'Dimer conversion degree:  ', dimer_conversion_degree(_dG, PCR_T)*100
    print ''
    
    print seq1, seq2
    print dimer_dG(dim1, seq1, seq2)
    print dimer_dG_corrected(dim1, seq1, seq2)
    print ''
    print seq3, seq4
    print dimer_dG(dim2, seq3, seq4)
    print dimer_dG_corrected(dim2, seq3, seq4)
    print ''
#    for i in range(-200,-30):
#        dG = i/10.0
#        print 'dG: %f' % dG
#        print 'K: %e'  % equilibrium_constant(dG, PCR_T)
#        print 'r: %f%%'  % (dimer_conversion_degree(dG, PCR_T)*100)
#        print ''
        
#    print
#    Pa = 0.8
#    Pb = 0.2
#    PA = 2*Pa
#    PB = 2*Pb
#    PaPbA = Pa*Pb 
#    PaPbB = Pa*Pb
#    for i in range(3,30):
#        PaPbA_cur = PaPbA 
#        PaPbA += Pa*PB + PaPbB
#        PaPbB += Pb*PA + PaPbA_cur
#        PA += Pa
#        PB += Pb
#        print 'N: %d; PA: %f; PaPbA: %f; SUM: %f, exp(Pa*Pb)+PA: %f; %% %f' % (i, PA, PaPbA, PA+PaPbA, 4*(Pa*Pb)*2**(i-2)+Pa*Pb*(i-1), 4*(Pa*Pb)*2**(i-2)/(PA+PaPbA)*100)
#        print 'N: %d; PB: %f; PaPbB: %f; SUM: %f, exp(Pa*Pb)+PB: %f; %% %f' % (i, PB, PaPbB, PB+PaPbB, 4*(Pa*Pb)*2**(i-2)+Pa*Pb*(i-1), 4*(Pa*Pb)*2**(i-2)/(PB+PaPbB)*100)
#        print 'A-B%%: %f' % ((PA+PaPbA-(PB+PaPbB))/(PA+PaPbA)*100)
#        print ''
        
        
