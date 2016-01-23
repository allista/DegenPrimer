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

from BioUtils.Tools.Multiprocessing import MPMain

class Equilibrium_Test(MPMain):
    def test(self): self()
    
    def _main(self):
        from DegenPrimer.WorkCounter import WorkCounter
        from DegenPrimer.Equilibrium import Equilibrium, EquilibriumSolver, EquilibriumBase
        import cProfile
        import shelve
        
        db = shelve.open('../data/reactions.shelve', 'r', protocol=-1)
        reactions = db['reactions']
        concentrations = db['concentrations']
        
        print len(reactions)
        print len(concentrations)
        
        cProfile.runctx('''for i in xrange(1):
        eq = EquilibriumSolver(self.abort_event, reactions, concentrations, 1e-10)
        eq.calculate()''', 
        globals(), locals(),
        'EquilibriumSolver.profile')
    
        print 'Done'
        
if __name__ == '__main__':
    Equilibrium_Test()