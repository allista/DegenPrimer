# coding=utf-8
#
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
Created on 2016-01-14

@author: Allis Tauri <allista@gmail.com>
'''

import DegenPrimer.TD_Functions as tdf
from DegenPrimer.SecStructures import Duplex

if __name__ == '__main__':
    tdf.PCR_P.Na = 50.0e-3
    tdf.PCR_P.Mg = 3.0e-3
    tdf.PCR_P.dNTP = 0.15e-6
    tdf.PCR_P.DNA = 1.0e-9
    tdf.PCR_P.DMSO = 0.0
    tdf.PCR_P.PCR_T = 60.0
    du = Duplex('GAACGCAAAGATCGGGAAC', 'CTTGCGTTTCTAACCCTTG'[::-1])
    print du
    
