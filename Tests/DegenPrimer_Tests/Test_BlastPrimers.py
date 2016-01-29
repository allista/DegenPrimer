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

if __name__ == '__main__':
    import sys, os
    from time import sleep
    from multiprocessing import Manager
    from DegenPrimer.Primer import Primer
    from DegenPrimer.SeqUtils import load_sequence
    from DegenPrimer.WorkCounter import WorkCounterManager
    from DegenPrimer.WaitingThread import WaitingThread
    from DegenPrimer.TD_Functions import PCR_P
    from DegenPrimer.BlastPrimers import BlastPrimers
    from threading import Lock
    
    os.chdir('../')
    
    mgr = Manager()
    abort_event = mgr.Event()
    
    PCR_P.PCR_T = 53
    PCR_P.Mg    = 3e-3
    PCR_P.dNTP  = 300e-6
    PCR_P.DNA   = 1e-10
    
    fwd_primer = Primer(load_sequence('ATATTCTACRACGGCTATCC', 'F-TGAM_0057-268_d1', 'F-TGAM_0057-268_d1'), 0.43e-6, True)
    rev_primer = Primer(load_sequence('GAASGCRAAKATYGGGAAC', 'R-TGAM_0055-624-d4', 'R-TGAM_0055-624-d4'), 0.43e-6, True)
    
    blastp = BlastPrimers(abort_event,
                          job_id='F-TGAM_0057-268_d1-R-TGAM_0055-624-d4', 
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
    if not blastp.have_saved_results(): sys.exit()

    cmgr = WorkCounterManager()
    cmgr.start()
    
    plock = Lock()    
    
    counter = cmgr.WorkCounter()
    job = WaitingThread(plock, 1, target=blastp.load_and_analyze, 
                        name='load_and_analyze', kwargs=dict(counter=counter))
    job.start(); print ''
    while job.is_alive():
        if counter.changed_percent():
            with plock: print counter #.rstr()
        sleep(0.1)
    job.join()
    with plock: print counter #.rstr()
    
#    blastp.write_reports()
