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

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from BioUtils.Tools.Multiprocessing import MPMain
from DegenPrimer.Primer import Primer
from DegenPrimer import TD_Functions
from DegenPrimer.AllSecStructures import AllSecStructures
from DegenPrimer.WorkCounter import WorkCounter

class AllSecStructurs_Test(MPMain):

    def test(self): return self()
    
    def _main(self):
        TD_Functions.PCR_P.PCR_T = 53
        fwd_primer  = Primer(SeqRecord(Seq('ATARTCTYCGAMGGCTATCC').reverse_complement(),id='primer1'), 0.9e-6, True)
        rev_primer  = Primer(SeqRecord(Seq('NAAGGYTTAGAKGCGGAAG').complement(),id='primer2'), 0.9e-6, True)
        print repr(fwd_primer)
        print repr(rev_primer)
        
        all_structs = AllSecStructures(self.abort_event, 'test', [fwd_primer, rev_primer])
        all_structs.find_structures()
    
        print all_structs.concentrations()
        print all_structs.reactions()
    
        c = WorkCounter()
        all_structs.calculate_equilibrium(c)
        print c.percent()
        
        print all_structs.equilibrium_concentrations()
        print ''
        print all_structs.print_structures()
        
if __name__ == '__main__':
    AllSecStructurs_Test()