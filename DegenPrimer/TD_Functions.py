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
from UnifiedNN import *
from StringTools import print_exception
try:
    from Bio.SeqFeature import SeqFeature, FeatureLocation
except Exception, e:
    print_exception(e)
    raise ImportError('The BioPython must be installed in your system.')


#utility functions
def print_exception(e):
    print "Exception occurred: " + str(type(e)) + " : " + e.__str__()
###############################################################################

#standard PCR conditions
C_Mg   = 1.5  #mM
C_Na   = 50   #mM; should be above 0.05M and below 1.1M
C_dNTP = 0    #mM
C_DNA  = 50   #nM; DNA template concentration
C_Prim = 0.1  #uM; Primer concentration
C_DMSO = 0    #percent
PCR_T  = 37   #C;  temperature at which PCR is conducted. Default is 37 for dG tables contain values of standard Gibbs energy at 37C

def C_Na_eq():
    """divalent cation correction (Ahsen et al., 2001)"""
    global C_Na, C_Mg, C_dNTP
    return C_Na + 120*sqrt(C_Mg - C_dNTP)
#end def


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
    global C_Prim, C_DNA, C_DMSO, R, K0, Sym_Correction
    seq_str = str(seq)
    rev_com = str(seq.reverse_complement())
    seq_len = len(seq)
    dH, dS = 0, 0
    #concentrations
    P   = C_Prim*1e-6
    D   = C_DNA *1e-9
    DUP = r*min(P,D)
    #equilibrium constant 
    K   = DUP/((P-DUP)*(D-DUP))
    #initial corrections
    dH += delta_H('ini', 'ini')
    dS += delta_S('ini', 'ini')
    #test for AT terminals
    if seq_str[0] == 'A' or seq_str[0] == 'T':
        dH += delta_H('ter', 'ter')
        dS += delta_S('ter', 'ter')
    if seq_str[-1] == 'A' or seq_str[-1] == 'T':
        dH += delta_H('ter', 'ter')
        dS += delta_S('ter', 'ter')
    #stacking interactions
    for n in range(len(seq_str)-1):
        NN  = seq_str[n:n+2]
        RC  = rev_com[seq_len-n-2:seq_len-n]
        dH += delta_H(NN, RC)
        dS += delta_S(NN, RC)
    #salt concentration correction
    dS = dS + dS_Na_coefficient * len(seq_str) * log(C_Na_eq()*1e-3) #C_Na mM
    #final temperature calculation
    return dH * 1000/(dS - R * log(K)) + K0 - 0.75 * C_DMSO #DMSO correction from [2]
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
                 len(str(C_DMSO)))
    conc_str += 'C(Na)     = ' + str(C_Na)  + ' '*(spacer-len(str(C_Na)))   +' mM\n'
    conc_str += 'C(Mg)     = ' + str(C_Mg)  + ' '*(spacer-len(str(C_Mg)))   +' mM\n'
    conc_str += 'C(dNTP)   = ' + str(C_dNTP)+ ' '*(spacer-len(str(C_dNTP))) +' mM\n'
    conc_str += 'C(DNA)    = ' + str(C_DNA) + ' '*(spacer-len(str(C_DNA)))  +' nM\n'
    conc_str += 'C(Primer) = ' + str(C_Prim)+ ' '*(spacer-len(str(C_Prim))) +' uM\n'
    conc_str += 'C(DMSO)   = ' + str(C_DMSO)+ ' '*(spacer-len(str(C_DMSO))) +' %\n'
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
    fwd_matches = dimer.fwd_matches()
    #e.g. 5'-(2 ,3 ,4 ,8 ,9 )-3'
    rev_matches = dimer.rev_matches()
    #e.g. 3-'(13,14,15,19,20)-5'
    seq_str = str(seq1)
    seq_len = len(seq_str)
    rev_str = str(seq2[::-1])
    rev_len = len(rev_str)
    dG_Na   = dG_Na_coefficient_oligo * 1 * log(C_Na_eq()*1e-3)
    dG = delta_G('ini', 'ini')
    #check for 'left' dangling end
    if   fwd_matches[0] == 0 and rev_matches[0] > 0: #3' dangling
        dG += DanglingNN[rev_str[rev_matches[0]]+'X'][rev_str[rev_matches[0]-1]]
    elif rev_matches[0] == 0 and fwd_matches[0] > 0: #5' dangling
        dG += DanglingNN['X'+seq_str[fwd_matches[0]]][seq_str[fwd_matches[0]-1]]
    #check for 'left' terminal mismatch
    elif fwd_matches[0] > 0 and rev_matches[0] > 0:
        dG += Terminal_mismatch_mean
    #check for 'left' terminal AT
    elif fwd_matches[0] == 0 and rev_matches[0] == 0:
        if seq_str[0] == 'A' or seq_str[0] == 'T':
            dG += delta_G('ter', 'ter')
    #check for 'right' dangling end
    if   fwd_matches[-1] == seq_len-1 and rev_matches[-1] < rev_len-1: #5' dangling
        dG += DanglingNN['X'+rev_str[rev_matches[-1]]][rev_str[rev_matches[-1]+1]]
    elif rev_matches[-1] == rev_len-1 and fwd_matches[-1] < seq_len-1: #3' dangling
        dG += DanglingNN[seq_str[fwd_matches[-1]]+'X'][seq_str[fwd_matches[-1]+1]]
    #check for 'right' terminal mismatch
    elif fwd_matches[-1]  < seq_len-1 and rev_matches[0] < rev_len-1:
        dG += Terminal_mismatch_mean
    #check for 'right' terminal AT
    elif fwd_matches[-1] == seq_len-1 and rev_matches[-1] == rev_len-1:
        if seq_str[-1] == 'A' or seq_str[-1] == 'T':
            dG += delta_G('ter', 'ter')
    #stacking and mismatches
    for i in range(len(fwd_matches)-1):
        f_match = fwd_matches[i]
        f_next  = fwd_matches[i+1]
        r_match = rev_matches[i]
        r_next  = rev_matches[i+1]
        #if either || or |x| or |xx|
        if f_next-f_match < 4:
            NN  = seq_str[f_match:f_match+2]
            RV  = rev_str[r_match:r_match+2]
            #salt-corrected dG
            dG += MismatchNN[NN][RV] + dG_Na
            #if ||
            if f_next-f_match == 1: continue
            #if |x| or |xx|
            elif f_next-f_match < 4:
                NN1 = rev_str[r_next-1:r_next+1][::-1]
                RV1 = seq_str[f_next-1:f_next+1][::-1]
                dG += MismatchNN[NN1][RV1] + dG_Na
                continue
        #loop
        elif f_next-f_match < 31:
            dG += loop_dG(f_next-f_match-1, 'I') + 2*Terminal_mismatch_mean
        else: pass
    return dG
#end def


def hairpin_dG(hairpin, seq):
    fwd_matches = hairpin.fwd_matches()
    #e.g. 5'-(2 ,3 ,4 ,8 ,9 )...-3'
    rev_matches = hairpin.rev_matches()
    #e.g  3'-(24,23,22,18,17)...-5'
    seq_str = str(seq)
    seq_len = len(seq_str)
    dG_Na   = dG_Na_coefficient_oligo * 1 * log(C_Na_eq()*1e-3)
    dG = delta_G('ini', 'ini')
    #check for 'left' dangling end
    if   fwd_matches[0] == 0 and rev_matches[0] < seq_len-1:
        dG += DanglingNN['X'+seq_str[rev_matches[0]]][seq_str[rev_matches[0]+1]]
    elif fwd_matches[0] > 0 and rev_matches[0] == seq_len-1:
        dG += DanglingNN['X'+seq_str[fwd_matches[0]]][seq_str[fwd_matches[0]-1]]
    #check for 'left' terminal mismatch
    elif fwd_matches[0] > 0 and rev_matches[0] < seq_len-1:
        dG += Terminal_mismatch_mean
    #check for 'left' terminal AT
    elif fwd_matches[0] == 0 and rev_matches[0] == seq_len-1:
        if seq_str[0] == 'A' or seq_str[0] == 'T':
            dG += delta_G('ter', 'ter')
    #stacking and mismatches
    for i in range(len(fwd_matches)-1):
        f_match = fwd_matches[i]
        f_next  = fwd_matches[i+1]
        r_match = rev_matches[i]
        r_next  = rev_matches[i+1]
        #if either || or |x| or |xx|
        if f_next-f_match < 4:
            NN  = seq_str[f_match:f_match+2]
            RV  = seq_str[r_match-1:r_match+1][::-1]
            #salt-corrected dG
            dG += MismatchNN[NN][RV] + dG_Na
            #if ||
            if f_next-f_match == 1: continue
            #if |x| or |xx|
            elif f_next-f_match < 4:
                NN1 = seq_str[r_next:r_next+2]
                RV1 = seq_str[f_next-1:f_next+1][::-1]
                dG += MismatchNN[NN1][RV1] + dG_Na
                continue
        #internal loop
        elif f_next-f_match < 31:
            dG += loop_dG(f_next-f_match-1, 'I') + 2*Terminal_mismatch_mean
        else: pass
    #hairpin loop
    hp_len = rev_matches[-1]-fwd_matches[-1]-1
    dG += loop_dG(hp_len, 'H')
    #3-4 loop
    if hp_len < 5:
        hp_str = seq_str[fwd_matches[-1]:rev_matches[-1]+1]
        if hp_str in Tri_Tetra_Loops:
            dG += Tri_Tetra_Loops[hp_str]
        if hp_len == 3:
            if seq_str[fwd_matches[-1]] == 'A' or seq_str[fwd_matches[-1]] == 'T':
                dG += 0.5 #kcal/mol; AT-closing penalty
        elif hp_len == 4:
            dG += Terminal_mismatch_mean
    else: dG += Terminal_mismatch_mean
    return dG
#end def


def equilibrium_constant(dG_T, T):
    '''calculate equilibrium constant of the annealing reaction
    at a given temperature, given standard dG at this temperature'''
    global R, K0
    return exp(-1000*dG_T/(R*(T-K0))) #annealing equilibrium constant


def conversion_degree(dG_T, T):
    '''calculate conversion degree at equilibrium
    given standard dG(kcal/mol) of annealing at T(C) temperature'''
    global C_Prim, C_DNA 
    K = equilibrium_constant(dG_T, T)
    P    = C_Prim*1e-6 #M
    D    = C_DNA *1e-9 #M
    #quadratic equation with respect to DUP = r*min(P,D), 
    #where 'r' is a conversion degree
    _b    = (K*P+K*D+1) #MINUS b; always positive
    disc = _b*_b - 4*K*(K*P*D) #this should always be >= 0 given non-negative K, P and D
    DUP  = (_b-sqrt(disc))/(2*K) #take the smallest positive root
    return DUP/min(P,D)
#end def