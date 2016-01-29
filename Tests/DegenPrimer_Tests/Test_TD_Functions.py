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

def test():
    from copy import deepcopy
    print tdf.PCR_P
    tdf.PCR_P.set({'PCR_T': 30})
    print tdf.PCR_P

    pcrp = deepcopy(tdf.PCR_P)
    print type(tdf.PCR_P), type(pcrp)
    tdf.PCR_P.set({'PCR_T': 60})
    print tdf.PCR_P
    print pcrp