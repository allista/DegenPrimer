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

from time import sleep
_pid = -1
abort_event = None
def sig_handler(signal, frame):
    if _pid != os.getpid(): return
    abort_event.set()
    sleep(0.1)
#end def


if __name__ == '__main__':
    from DegenPrimer.WorkCounter import WorkCounter
    from DegenPrimer.Equilibrium import Equilibrium
    #setup signal handler
    import signal
    signal.signal(signal.SIGINT,  sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)
    signal.signal(signal.SIGQUIT, sig_handler)
    
    import cProfile, os, sys
    from multiprocessing import Event
    import shelve
    
    abort_event = Event()
    _pid = os.getpid()
#    os.chdir('../')
    
    db = shelve.open('../results/Equilibrium-test.new', 'r', protocol=-1)
    reactions = db['reactions']
    concentrations = db['concentrations']
    
    print len(reactions)
    print len(concentrations)
    
    eq = Equilibrium(abort_event, reactions, concentrations, 1e-10)
    eq.calculate(WorkCounter())
    print eq.objective_value
    print eq.solution
    
    sys.exit(0)
    
    cProfile.run('''
for i in xrange(size):
    eq_system = Equilibrium(reactions, C_dict, 1e-10)
    solutions.append(eq_system.calculate())
''',
'equilibrium.profile')

    print 'Done'