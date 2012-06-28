# coding=utf-8

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

from math import sqrt, log
from Bio.SeqFeature import SeqFeature, FeatureLocation
from UnifiedNN import *

#utility functions
def print_exception(e):
    print "Exception occurred: " + str(type(e)) + " : " + e.__str__()
###############################################################################

#standard PCR conditions
C_Mg   = 1.5  #mM
C_Na   = 50   #mM; should be above 0.05M and below 1.1M
C_dNTP = 0    #mM
C_DNA  = 50   #nM; DNA template concentration
C_Prim = 0.11 #uM; Primer concentration

def C_Na_eq():
    """divalent cation correction (Ahsen et al., 2001)"""
    global C_Na, C_Mg, C_dNTP
    return C_Na + 120*sqrt(C_Mg - C_dNTP)
#end def


def NN_Tr(seq, r):
    """Calculate temperature for equilibrium with 'r' ratio 
    using Nearest Neighbor TD tables and two-state equilibrium
    equations from the paper of SantaLucia & Hicks (2004)"""
    global C_DNA, R, K0
    seq_str = str(seq)
    rev_com = str(seq.reverse_complement())
    seq_len = len(seq)
    dH, dS, K = 0, 0, 0
    #initial corrections
    dH += delta_H('ini', 'ini')
    dS += delta_S('ini', 'ini')
    #test for self-complementarity
    if seq_str == rev_com:
        dH += delta_H('sym', 'sym')
        dS += delta_S('sym', 'sym')
        #concentrations
        #C_Prim uM; C_DNA nM; DNA 2 strands
        A = C_Prim*1e-6
        B = 2*C_DNA*1e-9 #both strands
        if A*(1-r) <= B*r: 
            raise ValueError('For self-complementary oligonucleotides: for '
                             'equilibrium at %d%% ssPrimer concentration '
                             'should be grater than %f of target dsDNA concentration' % \
                             (100*r, r/(1-r)))
        #equilibrium constant
        #primer binds to itself as well as to both strands of target DNA
        #thus effective primer concentration and consequently Tm are decreased 
        K = r**2/(2*(1-r)) * A/(A*(1-r)-B*r)**2
    else:
        #concentrations
        A = C_Prim*1e-6
        B = C_DNA*1e-9
        if A <= B*r: 
            raise ValueError('For equilibrium at %d%% ssPrimer concentration '
                             'should be grater than %f of target dsDNA concentration' % \
                             (100*r, r))
        #equilibrium constant 
        #primer binds
        K  = r/((1-r)*(A-B*r))
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
    return dH * 1000/(dS - R * log(K)) + K0
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
    conc_str += 'C(Na)     = ' + str(C_Na)+  ' mM\n'
    conc_str += 'C(Mg)     = ' + str(C_Mg)+  ' mM\n'
    conc_str += 'C(dNTP)   = ' + str(C_dNTP)+' mM\n'
    conc_str += 'C(DNA)    = ' + str(C_DNA)+ ' nM\n'
    conc_str += 'C(Primer) = ' + str(C_Prim)+' uM\n'
    return conc_str
#end_def


def add_PCR_conditions(feature):
    try:
        feature.qualifiers['C_Na']      = str(C_Na)+  ' mM'
        feature.qualifiers['C_Mg']      = str(C_Mg)+  ' mM'
        feature.qualifiers['C_dNTP']    = str(C_dNTP)+' mM'
        feature.qualifiers['C_DNA']     = str(C_DNA)+ ' nM'
        feature.qualifiers['C_Primer']  = str(C_Prim)+' uM'
    except Exception, e:
        print_exception(e)
#end def


def calculate_Tr(seq_rec, r):
    try:
        primer_Tr = NN_Tr(seq_rec.seq, r)
        feature = source_feature(seq_rec)
        add_PCR_conditions(feature)
        feature.qualifiers['T-'+str(r)] = str(primer_Tr)
        return str(primer_Tr)
    except Exception, e:
        print_exception(e)
        return "unknown"
#end def


def calculate_Tm(seq_rec):
    try:
        primer_Tm = NN_Tm(seq_rec.seq)
        feature = source_feature(seq_rec)
        add_PCR_conditions(feature)
        feature.qualifiers['Tm'] = str(primer_Tm)
        return str(primer_Tm)
    except Exception, e:
        print_exception(e)
        return "unknown"
#end def


def dimer_dG(dimer, seq1, seq2):
    fwd_matches = list(dimer[0])
    fwd_matches.sort()
    #e.g. (2 ,3 ,4 ,8 ,9 )
    rev_matches = list(dimer[1])
    rev_matches.sort()
    #e.g. (13,14,15,19,20)
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
    fwd_matches = list(hairpin[0])
    fwd_matches.sort()
    #e.g. (2 ,3 ,4 ,8 ,9 )
    rev_matches = list(hairpin[1])
    rev_matches.sort(reverse=True)
    #e.g  (24,23,22,18,17)
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