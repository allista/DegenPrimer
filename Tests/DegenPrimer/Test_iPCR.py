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

class _Main(MPMain):
    def _main(self):
        import time
        import DegenPrimer.TD_Functions as tdf
        from multiprocessing import Manager
        from DegenPrimer.Primer import Primer
        from DegenPrimer.SeqUtils import load_sequence
        from DegenPrimer.WorkCounter import WorkCounterManager
        from BioUtils.Tools import WaitingThread
        from DegenPrimer.iPCR import iPCR
        from threading import Lock
        
        mgr = Manager()
        abort_event = mgr.Event()
        
        tdf.PCR_P.PCR_T = 53
        tdf.PCR_P.Mg    = 3e-3
        tdf.PCR_P.dNTP  = 300e-6
        tdf.PCR_P.DNA   = 1e-10
        fwd_primer = Primer(load_sequence('ATATTCTACRACGGCTATCC', 'fwd_test', 'fwd_test'), 0.43e-6, True)
        rev_primer = Primer(load_sequence('GAASGCRAAKATYGGGAAC', 'rev_test', 'rev_test'), 0.43e-6, True)
        ipcr = iPCR(abort_event,
                    max_mismatches=6,
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
                    include_side_annealings=True)
        cmgr = WorkCounterManager()
        cmgr.start()
        counter = cmgr.WorkCounter()
        
        plock = Lock()
            
        job = WaitingThread(plock, 1, target=ipcr.simulate_PCR, 
                            name='simulate_PCR',
                            args=(counter, 
                                  (
                                   '../data/ThGa.fa', #single sequence
    #                               '../data/Ch5_gnm.fa', '../data/ThDS1.fa', '../data/ThES1.fa', #long sequences 
                                   '../data/ThDS1-FC.fa', '../data/ThDS1-850b-product.fa', #short sequences
                                  ),))
        job.start(); print ''
        while job.is_alive():
            if counter.changed_percent():
                with plock: print counter
            time.sleep(1)
        job.join()
        with plock: print counter
        
        ipcr.write_report()

if __name__ == '__main__':
    import sys
    main = _Main()
    sys.exit(main())
    
############################## statistics ######################################
#TD_Functions.PCR_P.PCR_T = 53
#TD_Functions.PCR_P.Mg    = 3e-3
#TD_Functions.PCR_P.dNTP  = 300e-6
#TD_Functions.PCR_P.DNA   = 1e-10
#(abort_event,
#max_mismatches=7,
#job_id='test-job', 
#primers=[#fwd_primer,
#         rev_primer], 
#min_amplicon=50, 
#max_amplicon=2000, 
#polymerase=40000, 
#with_exonuclease=False, 
#num_cycles=30,
#side_reactions=None, 
#side_concentrations=None,
#include_side_annealings=False)

#Results use: 44.114922Mb
#PCR Simulation: searching for possible PCR products in Thermococcus gammatolerans EJ3...
#Results use: 29.034309Mb
#PCR Simulation: searching for possible PCR products in TBCH5...
#Results use: 41.029671Mb
#PCR Simulation: searching for possible PCR products in DS1 complete temp...
#Results use: 21.963631Mb
#PCR Simulation: searching for possible PCR products in Thermococcus sp. ES1...