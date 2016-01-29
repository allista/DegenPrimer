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
    from multiprocessing import Manager
    from threading import Lock
    from BioUtils.Tools import WaitingThread
    
    from DegenPrimer.Primer import Primer
    from DegenPrimer.SeqUtils import load_sequence
    from DegenPrimer.WorkCounter import WorkCounterManager
    from DegenPrimer.PCR_Optimizer import PCR_Optimizer
    from DegenPrimer import TD_Functions as tdf
    
    mgr = Manager()
    abort_event = mgr.Event()
    
    tdf.PCR_P.PCR_T = 53
    tdf.PCR_P.Mg    = 3e-3
    tdf.PCR_P.dNTP  = 300e-6
    tdf.PCR_P.DNA   = 1e-10
    fwd_primer = Primer(load_sequence('ATATTCTACRACGGCTATCC', 'fwd_test', 'fwd_test'), 0.43e-6, True)
    rev_primer = Primer(load_sequence('GAASGCRAAKATYGGGAAC', 'rev_test', 'rev_test'), 0.43e-6, True)
    optimizer  = PCR_Optimizer(abort_event,
                               100, 5,
                               max_mismatches=5,
                               job_id='test-job', 
                               primers=[fwd_primer,
                                        rev_primer], 
                               min_amplicon=50, 
                               max_amplicon=2000, 
                               polymerase=40000, 
                               with_exonuclease=False, 
                               num_cycles=30,
                               side_reactions=None, 
                               side_concentrations=None,
                               include_side_annealings=False)
    cmgr = WorkCounterManager()
    cmgr.start()
    counter = cmgr.WorkCounter()
    
    plock = Lock()
    
    job = WaitingThread(plock, 1, target=optimizer.optimize_PCR_parameters, 
                        name='optimize PCR',
                        args=('../data/ThGa.fa', 
                              ({'start':47920, 'end':49321},), 
                              (
                               {'name':'PCR_T',
                                'min':50, 'ini':60, 'max':72},
#                               {'name':'dNTP',
#                                'min':100e-6, 'ini':200e-6, 'max':900e-6},
                               ),
                             ))
    job.start(); print ''
    job.join()

    optimizer.write_report()