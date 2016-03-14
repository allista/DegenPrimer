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

def test():
    import cProfile
    import DegenPrimer.TD_Functions as tdf
    from DegenPrimer.SecStructures import Duplex, reverse_complement, Dimer

    tdf.PCR_P.Na = 50.0e-3
    tdf.PCR_P.Mg = 3.0e-3
    tdf.PCR_P.dNTP = 0.15e-6
    tdf.PCR_P.DNA = 1.0e-9
    tdf.PCR_P.DMSO = 0.0
    tdf.PCR_P.PCR_T = 60.0
    with tdf.AcquireParameters():
        du = Duplex('AGAGAACGCAAAGATCGGGAAC', 'CTTGCGTTTCTAACCCTTG'[::-1], dimer=Dimer((3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21), 3))
        print du
        cProfile.runctx('for x in xrange(100000): du.print_most_stable()', 
                        globals(), locals(), 'Duplex.print_stable.profile')
    
#    seq = 'ATGCGTCACTACCAGT'*10000
#    cProfile.runctx('''for x in xrange(100): 
#    reverse_complement(seq)''', 
#                    globals(), locals(), 'reverse_complement.profile')

test()