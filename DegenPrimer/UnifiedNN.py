# coding=utf-8
#
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

All calculations and data are based on:
﻿SantaLucia, J., & Hicks, D. (2004). 
The thermodynamics of DNA structural motifs. Annual review of biophysics and 
biomolecular structure, 33, 415-40. doi:10.1146/annurev.biophys.32.110601.141800
'''


from math import log
from Tri_Tetra_Loops import Tri_Tetra_Loops


#constants
R                 =  1.9872 #Universal gas constant cal/(K*mol)
K0                = -273.15 #Absolute temperature zero in degree Celsius
dG_Na_coefficient_oligo = -0.114 #kcal/mol for oligomers with length =< 16
dG_Na_coefficient_poly  = -0.175 #kcal/mol for longer polymers
dS_Na_coefficient = +0.368 #e.u.
Loop_coefficient  =  2.44
#The NN stabilities at 37◦ C range from −1.23 to −0.21 kcal/mol 
#for CG/GA and AC/TC, respectively
Terminal_mismatch_mean  = (-1.23 + -0.21)/2


UnifiedNN = {
            #'NN' :{'dH': kcal/mol, 'dS': e.u., 'dG'(37C): kcal/mol},
             'AA' :{'dH':  -7.6, 'dS': -21.3, 'dG': -1.00}, #TT
             'AT' :{'dH':  -7.2, 'dS': -20.4, 'dG': -0.88}, #TA
             'TA' :{'dH':  -7.2, 'dS': -21.3, 'dG': -0.58}, #AT
             'CA' :{'dH':  -8.5, 'dS': -22.7, 'dG': -1.45},
             'GT' :{'dH':  -8.4, 'dS': -22.4, 'dG': -1.44},
             'CT' :{'dH':  -7.8, 'dS': -21.0, 'dG': -1.28},
             'GA' :{'dH':  -8.2, 'dS': -22.2, 'dG': -1.30},
             'CG' :{'dH': -10.6, 'dS': -27.2, 'dG': -2.17},
             'GC' :{'dH':  -9.8, 'dS': -24.4, 'dG': -2.24},
             'GG' :{'dH':  -8.0, 'dS': -19.9, 'dG': -1.84},
             
             'ini':{'dH':  +0.2, 'dS':  -5.7, 'dG': +1.96}, #initiation penalty
             'ter':{'dH':  +2.2, 'dS':  +6.9, 'dG': +0.05}, #terminal A-T pair penalty
             'sym':{'dH':   0.0, 'dS':  -1.4, 'dG': +0.43}, #self-complementary sequence penalty
}


MismatchNN = {
             #dG(37C) for sequence with single mismatches
             'GA' :{'CA':  0.17, 'CC':  0.81, 'CG': -0.25, 'CT':  UnifiedNN['CT']['dG']},
             'GC' :{'CA':  0.47, 'CC':  0.79, 'CG':  UnifiedNN['CG']['dG'], 'CT':  0.62},
             'GG' :{'CA': -0.52, 'CC':  UnifiedNN['GG']['dG'], 'CG': -1.11, 'CT':  0.08},
             'GT' :{'CA':  UnifiedNN['CA']['dG'], 'CC':  0.98, 'CG': -0.59, 'CT':  0.45},
             
             'CA' :{'GA':  0.43, 'GC':  0.75, 'GG':  0.03, 'GT':  UnifiedNN['GT']['dG']},
             'CC' :{'GA':  0.79, 'GC':  0.70, 'GG':  UnifiedNN['GG']['dG'], 'GT':  0.62},
             'CG' :{'GA':  0.11, 'GC':  UnifiedNN['GC']['dG'], 'GG': -0.11, 'GT': -0.47},
             'CT' :{'GA':  UnifiedNN['GA']['dG'], 'GC':  0.40, 'GG': -0.32, 'GT': -0.12},
             
             'AA' :{'TA':  0.61, 'TC':  0.88, 'TG':  0.14, 'TT':  UnifiedNN['AA']['dG']},
             'AC' :{'TA':  0.77, 'TC':  1.33, 'TG':  UnifiedNN['CA']['dG'], 'TT':  0.64},
             'AG' :{'TA':  0.02, 'TC':  UnifiedNN['GA']['dG'], 'TG': -0.13, 'TT':  0.71},
             'AT' :{'TA':  UnifiedNN['TA']['dG'], 'TC':  0.73, 'TG':  0.07, 'TT':  0.69},
             
             'TA' :{'AA':  0.69, 'AC':  0.92, 'AG':  0.42, 'AT':  UnifiedNN['AT']['dG']},
             'TC' :{'AA':  1.33, 'AC':  1.05, 'AG':  UnifiedNN['CT']['dG'], 'AT':  0.97},
             'TG' :{'AA':  0.74, 'AC':  UnifiedNN['GT']['dG'], 'AG':  0.44, 'AT':  0.43},
             'TT' :{'AA':  UnifiedNN['AA']['dG'], 'AC':  0.75, 'AG':  0.34, 'AT':  0.68},
}


DanglingNN = {
             #dG(37C) kcal/mol for dangling terminals
             #5' dangling ends
             'XA' :{'A': -0.51, 'C': -0.42, 'G': -0.62, 'T': -0.71},
             'XC' :{'A': -0.96, 'C': -0.52, 'G': -0.72, 'T': -0.58},
             'XG' :{'A': -0.58, 'C': -0.34, 'G': -0.56, 'T': -0.61},
             'XT' :{'A': -0.50, 'C': -0.02, 'G':  0.48, 'T': -0.10},
             #3' dangling ends
             'AX' :{'A': -0.12, 'C':  0.28, 'G': -0.01, 'T':  0.13},
             'CX' :{'A': -0.82, 'C': -0.31, 'G': -0.01, 'T': -0.52},
             'GX' :{'A': -0.92, 'C': -0.23, 'G': -0.44, 'T': -0.35},
             'TX' :{'A': -0.48, 'C': -0.19, 'G': -0.50, 'T': -0.29},
}


LoopNN = {
         #dG(37C) kcal/mol increment for loops of different size
         3  : {'I': 3.2, 'H': 3.5},
         4  : {'I': 3.6, 'H': 3.5},
         5  : {'I': 4.0, 'H': 3.3},
         6  : {'I': 4.4, 'H': 4.0},
         7  : {'I': 4.6, 'H': 4.2},
         8  : {'I': 4.8, 'H': 4.3},
         9  : {'I': 4.9, 'H': 4.5},
         10 : {'I': 4.9, 'H': 4.6},
         12 : {'I': 5.2, 'H': 5.0},
         14 : {'I': 5.4, 'H': 5.1},
         16 : {'I': 5.6, 'H': 5.3},
         18 : {'I': 5.8, 'H': 5.5},
         20 : {'I': 5.9, 'H': 5.7},
         25 : {'I': 6.3, 'H': 6.1},
         30 : {'I': 6.6, 'H': 6.3},
}


def delta_Par(seq, rev_comp, par):
    if seq in UnifiedNN:
        return UnifiedNN[seq][par]
    elif rev_comp in UnifiedNN:
        return UnifiedNN[rev_comp][par]
#end def

def delta_G(seq, rev_comp): return delta_Par(seq, rev_comp, 'dG')
def delta_H(seq, rev_comp): return delta_Par(seq, rev_comp, 'dH')
def delta_S(seq, rev_comp): return delta_Par(seq, rev_comp, 'dS')


def loop_dG(length, loop_type):
    if length > 30: raise ValueError('Loop length should not exceed 30')
    if length in LoopNN:
        return LoopNN[length][loop_type]
    else:
        exp_len = length
        while exp_len not in LoopNN: exp_len += 1
        return LoopNN[exp_len][loop_type] + Loop_coefficient * R * 310.15/1000 * log(float(length)/exp_len)
#end def