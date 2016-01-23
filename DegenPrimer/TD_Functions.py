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

from contextlib import contextmanager
from copy import deepcopy
from math import sqrt, log, exp

from .PCR_Parameters import PCR_Parameters
from .UnifiedNN import UnifiedNN
from BioUtils.Tools import Text


###############################################################################
#Unified Nearest Neighbour model and PCR parameters initialization
NN = UnifiedNN()
if not NN: raise RuntimeError('TD_Functions: Unable to initialize UnifiedNN.')
PCR_P  = PCR_Parameters()
_PCR_P = None

@contextmanager
def AcquireParameters():
    global _PCR_P, PCR_P
    assert _PCR_P is None, 'TD_Functions: PCR parameters are already acquired.'
    try: 
        _PCR_P, PCR_P = PCR_P, deepcopy(PCR_P)
        yield
    finally: 
        PCR_P, _PCR_P = _PCR_P, None
        
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
    #concentrations
    P   = concentration
    D   = PCR_P.DNA
    DUP = conversion_degree*min(P,D)
    #equilibrium constant 
    K   = DUP/((P-DUP)*(D-DUP))
    #initial corrections
    dP = NN.ini
    dH = dP.dH; dS = dP.dS
    #test for AT terminals
    if seq_str[0]  in 'AT':
        dP  = NN.ter
        dH += dP.dH; dS += dP.dS
    if seq_str[-1] in 'AT':
        dP  = NN.ter
        dH += dP.dH; dS += dP.dS
    #stacking interactions
    for n in xrange(len(seq_str)-1):
        pair    = seq_str[n:n+2]
        reverse = rev_str[n:n+2]
        dP  = NN.int_dPar_37(pair, reverse)
        dH += dP.dH; dS += dP.dS
    #salt concentration correction
    dS = dS + NN.dS_Na_coefficient * len(seq_str) * PCR_P.Na_eq_log #C_Na qM
    #final temperature calculation
    return NN.K0 + dH * 1000/(dS - NN.R * log(K)) - NN.T_DMSO_coefficient * PCR_P.DMSO #DMSO correction from [2]
#end def

def primer_template_Tm(sequence, concentration): 
    return primer_template_Tr(sequence, concentration, 0.5)


def format_concentration(concentration):
    return Text.format_quantity(concentration, 'M') 
    
    
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
        conditions.append(['C(Poly)', ]+Text.format_quantity(polymerase*1e-6, 'u/ul').split())
    conditions.append(['T', '%.1f' % PCR_P.PCR_T, 'C'])
    return Text.print_table(conditions, delimiter='')
#end_def


def dimer_dG_corrected(dimer, seq1, seq2):
    '''calculate dG of a dimer corrected to current PCR conditions including 
    salt concentrations and temperature.
    dimer -- an instance of Dimer class which represents a dimer structure.
    seq1, seq2 -- sequences constituting a dimer, given in 5'->3' orientation.'''
    fwd_matches = dimer.fwd_matches #e.g. 5'-(2 ,3 ,4 ,8 ,9 )-3'
    offset  = dimer.offset
    seq_str = str(seq1)
    seq_len = len(seq_str)
    rev_str = str(seq2)[::-1]
    rev_len = len(rev_str)
    #initial dG of annealing
    dG = 0
    dP = NN.ini
    dH = dP.dH
    dS = dP.dS
    #check for 'left' dangling end
    fwd0 = fwd_matches[0] 
    if fwd0 == 0 and offset < 0: #3' dangling
        rev0    = fwd0-offset
        pair    = rev_str[rev0]+rev_str[rev0-1]
        reverse = seq_str[0]+'-' 
        dP  = NN.end_dPar_37(pair, reverse)
        dH += dP.dH; dS += dP.dS
    #check for 'left' terminal AT
    elif fwd0 == 0 and offset == 0:
        if seq_str[0] in 'AT':
            dP  = NN.ter
            dH += dP.dH; dS += dP.dS
    elif fwd0 == offset: #5' dangling
        pair    = seq_str[fwd0-1]+seq_str[fwd0]
        reverse = '-'+rev_str[0]
        dP  = NN.end_dPar_37(pair, reverse)
        dH += dP.dH; dS += dP.dS 
    #check for 'left' terminal mismatch
    elif fwd0 > offset:
        dG += NN.Terminal_mismatch_mean
    #check for 'right' dangling end
    fwd_1 = fwd_matches[-1]
    rev_1 = fwd_1-offset
    if fwd_1 == seq_len-1 and rev_1 < rev_len-1: #5' dangling
        pair    = rev_str[rev_1+1]+rev_str[rev_1]
        reverse = '-'+seq_str[-1] 
        dP  = NN.end_dPar_37(pair, reverse)
        dH += dP.dH; dS += dP.dS
    elif rev_1 == rev_len-1 and fwd_1 < seq_len-1: #3' dangling
        pair    = seq_str[fwd_1]+seq_str[fwd_1+1]
        reverse = rev_str[-1]+'-' 
        dP  = NN.end_dPar_37(pair, reverse)
        dH += dP.dH; dS += dP.dS
    #check for 'right' terminal mismatch
    elif fwd_1 < seq_len-1 and rev_1 < rev_len-1:
        dG += NN.Terminal_mismatch_mean
    #check for 'right' terminal AT
    elif fwd_1 == seq_len-1 and rev_1 == rev_len-1:
        if seq_str[-1] in 'AT':
            dP  = NN.ter
            dH += dP.dH; dS += dP.dS
    #stacking and mismatches
    for i in xrange(len(fwd_matches)-1):
        f_match = fwd_matches[i]
        f_next  = fwd_matches[i+1]
        r_match = f_match-offset
        r_next  = f_next-offset
        gap     = f_next-f_match
        #if either || or |x| or |xx|
        if gap < 4:
            pair    = seq_str[f_match:f_match+2]
            reverse = rev_str[r_match:r_match+2]
            #dH and salt-corrected dS
            dP  = NN.int_dPar_37(pair, reverse)
            dH += dP.dH; dS += dP.dS
            #if |x| or |xx|
            if gap > 1:
                pair1    = rev_str[r_next-1:r_next+1][::-1]
                reverse1 = seq_str[f_next-1:f_next+1][::-1]
                dP  = NN.int_dPar_37(pair1, reverse1)
                dH += dP.dH; dS += dP.dS
        #loop
        else:
            loop_seq = seq_str[f_match:f_next+1]
            dS += NN.internal_loop_dS_37(loop_seq)
    #dS salt correction
    dS += NN.dS_Na_coefficient * (fwd_matches[-1]-fwd_matches[0]) * PCR_P.Na_eq_log
    return dG + dH - (PCR_P.PCR_T-NN.K0)*dS/1000.0
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
    dP = NN.ini
    dH = dP.dH
    dS = dP.dS
    #check for 'left' dangling end
    if   fwd_matches[0] == 0 and rev_matches[0] < seq_len-1: #3'-end
        pair    = seq_str[rev_matches[0]]+seq_str[rev_matches[0]+1]
        reverse = seq_str[0] + '-'
        dP  = NN.end_dPar_37(pair, reverse)
        dH += dP.dH; dS +=  dP.dS
    elif fwd_matches[0] > 0 and rev_matches[0] == seq_len-1: #5'-end
        pair    = seq_str[fwd_matches[0]-1]+seq_str[fwd_matches[0]]
        reverse = '-' + seq_str[-1]
        dP  = NN.end_dPar_37(pair, reverse)
        dH += dP.dH; dS += dP.dS
    #check for 'left' terminal mismatch
    elif fwd_matches[0] > 0 and rev_matches[0] < seq_len-1:
        dG += NN.Terminal_mismatch_mean
    #check for 'left' terminal AT
    elif fwd_matches[0] == 0 and rev_matches[0] == seq_len-1:
        if seq_str[0] in 'AT':
            dP  = NN.ter
            dH += dP.dH; dS += dP.dS
    #stacking and mismatches
    for i in xrange(len(fwd_matches)-1):
        f_match = fwd_matches[i]
        f_next  = fwd_matches[i+1]
        r_match = rev_matches[i]
        r_next  = rev_matches[i+1]
        gap     = f_next-f_match
        #if either || or |x| or |xx|
        if gap < 4:
            pair    = seq_str[f_match:f_match+2]
            reverse = seq_str[r_match-1:r_match+1][::-1]
            #dH and salt-corrected dS
            dP  = NN.int_dPar_37(pair, reverse)
            dH += dP.dH; dS += dP.dS
            #if |x| or |xx|
            if gap > 1:
                pair1    = seq_str[r_next:r_next+2]
                reverse1 = seq_str[f_next-1:f_next+1][::-1]
                dP  = NN.int_dPar_37(pair1, reverse1)
                dH += dP.dH; dS += dP.dS
        #loop
        else:
            loop_seq = seq_str[f_match:f_next+1]
            dS += NN.internal_loop_dS_37(loop_seq)
    #hairpin loop
    hp_str = seq_str[fwd_matches[-1]:rev_matches[-1]+1]
    dH += NN.hairpin_loop_dH_37(hp_str)
    dS += NN.hairpin_loop_dS_37(hp_str)
    #dS salt correction
    dS += NN.dS_Na_coefficient * (fwd_matches[-1]-fwd_matches[0]) * PCR_P.Na_eq_log
    return dG + dH - (PCR_P.PCR_T-NN.K0)*dS/1000.0
#end def


def equilibrium_constant(dG_T, T):
    '''calculate equilibrium constant of annealing at a given temperature, 
    given standard dG at this temperature'''
    return exp(-1000*dG_T/(NN.R*(T-NN.K0))) #annealing equilibrium constant
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