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
Created on Mar 15, 2013

@author: Allis Tauri <allista@gmail.com>
'''


import numpy as np
from scipy.optimize import fmin_l_bfgs_b
from StringTools import wrap_text
from BioUtils.Tools.Multiprocessing import MultiprocessingBase
from AllSecStructures import AllSecStructures
from PCR_ProductsFinder import PPFManager
from SinglePCR import SinglePCR
from Product import Region
from WorkCounter import WorkCounter
from BioUtils.Tools.Output import OutIntercepter
from BioUtils.Tools.tmpStorage import register_tmp_file, cleanup_file
from iPCR_Base import iPCR_Base
import TD_Functions as tdf


#manager to isolate forking point in a low-memory process
from PCR_Mixture import ShelvedMixture
from BioUtils.Tools.Multiprocessing import Parallelizer
from BioUtils.Tools.UMP import FuncManager, at_manager

class compute_objective_value(object):
    def __init__(self, pcr, product_bounds, purity):
        self._pcr            = pcr
        self._product_bounds = product_bounds
        self._purity         = purity
    #end def
    
    @staticmethod
    def _abs_prod_C(pcr_results, bounds):
        r_id, _ov, products, _r_end = pcr_results
        lookup_region = Region(r_id, bounds['start'], bounds['end'], forward=True)
        total_concentration   = 0
        product_concentration = -1
        for _p_id, product in products.iteritems():
            total_concentration += product.quantity
            if product == lookup_region:
                product_concentration = product.quantity
        if total_concentration == 0 or product_concentration < 0: return 0, 0
        return (product_concentration, 
                product_concentration/total_concentration)
    #end def
    
    def _compute_prod_C(self, pcr_results):
        result = 0
        for bounds in self._product_bounds:
            prod_C  = self._abs_prod_C(pcr_results, bounds)
            result += -1e9*prod_C[0]*(prod_C[1]**self._purity)
        return result
    #end def
        
    def __call__(self, mixture_path):
        pcr_results = self._pcr(ShelvedMixture(mixture_path))
        if pcr_results is None: return 1
        return self._compute_prod_C(pcr_results)
#end class

def _parallel_OV(abort_event, OVCs, mixture_paths):
    parallelize = at_manager(Parallelizer, 'parallelize_both')
    return parallelize(abort_event, False, 0.1, OVCs, mixture_paths)
#end def

OVC_Manager   = FuncManager('OVC_Manager', (_parallel_OV,))
SearchManager = FuncManager('SearchManager', 
                            (at_manager(PPFManager, 'matches_to_mixture'), 
                             at_manager(PPFManager, 'find_matches')))


class PCR_Optimizer(MultiprocessingBase, iPCR_Base):
    '''
    Repeatedly perform PCR simulation varying parameters 
    to find optimal conditions
    '''

    _PCR_report_suffix = 'PCR-optimization'
    _epsilon = 1e-8 #dx for finite-difference approximation of the gradient
    
    def __init__(self, 
                  abort_event, max_steps, purity, 
                  *args, **kwargs):
        MultiprocessingBase.__init__(self, abort_event)
        iPCR_Base.__init__(self, abort_event, *args, **kwargs)
        #PCR simulations
        self._PCR_ProductsFinder = self._new_PCR_ProductsFinder()
        with OutIntercepter(): #to suppress output from SearchManager
            self._searcher       = SearchManager()
            self._searcher.start()
        self._ovcM               = OVC_Manager()
        self._ovcM.start()
        self._last_OV            = None
        self._last_mixture       = None
        self._last_parameters    = None
        #target
        self._sec_structures     = None
        self._matches_list       = None
        self._product_bounds     = None
        self._seq_name           = None
        #parameters
        self._num_params         = None
        self._iscales            = None
        self._iparams            = None
        self._ibounds            = None
        self._param_dict         = None
        self._relative           = None 
        #solution and flags
        self._max_steps          = max_steps
        self._purity             = purity
        self._exit_status        = None
        self._optimized          = False
        self._optimum            = None 
        self._named_optimum      = None
    #end def
    
    
    def _name_parameters(self, parameters):
        named_parameters = dict()
        for name, i in self._param_dict.items():
            named_parameters[name] = parameters[i]*self._iscales[i]+self._iscales[i]/2
        return named_parameters
    #end def
    
    def _set_parameters(self, parameters):
        #make parameter dictionary
        named_parameters = self._name_parameters(parameters)
        #set thermodynamic parameters
        tdf.PCR_P.set(named_parameters)
        #set primer concentrations
        for primer in self._primers:
            if primer.id in named_parameters:
                primer.total_concentration = named_parameters[primer.id]
        #change polymerase concentration
        if 'polymerase' in named_parameters:
            self._polymerase = named_parameters['polymerase']
    #end def
    
    
    def _update_side_reactions(self, parameters):
        sec_structs = AllSecStructures(self._abort_event, self._job_id, self._primers)
        sec_structs.find_structures()
        self._side_concentrations = sec_structs.concentrations()
        self._side_reactions      = sec_structs.reactions()
    #end def
    
    def _prepare_simulation(self, parameters):
        self._set_parameters(parameters)
        #recalculate PCR mixture
        mixture = self._searcher.matches_to_mixture(WorkCounter(), 
                                                    self._seq_name, 
                                                    self._matches_list,
                                                    self._PCR_ProductsFinder)
        if mixture is None: return None, None
        register_tmp_file(mixture)
        self._update_side_reactions(parameters)
        pcr = SinglePCR(self._abort_event, 
                        self._primers_concentrations, 
                        self._elongation_time, 
                        self._polymerase, 
                        self._with_exonuclease, 
                        self._num_cycles, 
                        self._side_reactions, self._side_concentrations)
        ovc = compute_objective_value(pcr, self._product_bounds, self._purity)
        return ovc, mixture
    #end def
    
    def _objective_function(self, parameters):
        self._last_parameters = parameters
        ovc, mixture = self._prepare_simulation(parameters)
        if ovc is None: return -1e32 #FIXME: end search process if aborted
        if self._last_mixture: cleanup_file(self._last_mixture)
        self._last_mixture = mixture
        self._last_OV = ovc(mixture)
        return self._last_OV
    #end def
    
    def _objective_function_grad(self, parameters):
        '''Should always be called right after _objective_function call'''
        assert (parameters == self._last_parameters).all(), \
        'Grad argument differs from last_parameters'
        p0   = parameters
        f0   = self._last_OV
        grad = np.zeros((self._num_params,), float)
        ei   = np.zeros((self._num_params,), float)
        ovcs = []; mixtures = []
        for p in xrange(self._num_params): #each step is highly parallelized
            ei[p] = self._epsilon; pi = p0+ei
            ovc, mixture = self._prepare_simulation(pi)
            if ovc is None: return grad #FIXME: need to handle aborts somehow
            ovcs.append(ovc); mixtures.append(mixture)
            ei[p] = 0.0
        fi = self._ovcM.parallel_OV(self._abort_event, ovcs, mixtures)
        if fi is None: return grad #FIXME: need to handle aborts somehow
        map(cleanup_file, mixtures)
        fi = np.asarray(fi, float)
        dx = np.fromfunction(lambda x: self._epsilon, (self._num_params,))
        return (fi-f0)/dx
    #end def
    
    
    def _prepare_parameters(self, product_bounds, parameters):
        self._product_bounds = product_bounds
        self._iscales        = []
        self._iparams        = []
        self._ibounds        = []
        self._param_dict     = dict()
        for par in parameters:
            self._param_dict[par['name']] = len(self._iparams)
            self._iscales.append(float(par['max']-par['min']))
            self._iparams.append((par['ini']-self._iscales[-1]/2)/self._iscales[-1])
            self._ibounds.append(((par['min']-self._iscales[-1]/2)/self._iscales[-1], 
                                  (par['max']-self._iscales[-1]/2)/self._iscales[-1]))
        self._num_params = len(self._iparams)
    #end def

    def _find_matches(self, seq_file, seq_id=None):
        #try to connect to a database
        if not self._try_connect_db((seq_file,)): return False
        #get names of the sequences
        self._seq_names = self._get_names((seq_id,) if seq_id else None)
        if not self._seq_names:
            self._seq_db.close() 
            return False
        self._seq_name = self._seq_names.values()[0]
        #get sequences
        templates = self._seq_db.get_seqs(self._seq_names.keys())
        self._seq_db.close()
        #find match positions
        self._matches_list = self._searcher.find_matches(WorkCounter(), 
                                                         self._seq_name, templates[0][2], 
                                                         self._max_mismatches, 
                                                         self._PCR_ProductsFinder)
        return self._matches_list is not None
    #end def

    def optimize_PCR_parameters(self, seq_file, product_bounds, 
                                   parameters, seq_id=None):
        self._optimized = False
        #prepare parameters and matchers
        if not (seq_file and product_bounds and parameters): return False
        if not self._find_matches(seq_file, seq_id): return False
        self._prepare_parameters(product_bounds, parameters)
        #optimize parameters
        print ''
        optimum, _unused, info = fmin_l_bfgs_b(self._objective_function, 
                                               self._iparams,
                                               self._objective_function_grad, 
                                               approx_grad=False, 
                                               bounds=self._ibounds, factr=1e12,
                                               maxfun=self._max_steps, disp=1)
        if self._last_mixture: cleanup_file(self._last_mixture)
        #parse results
        if info['warnflag'] == 2:
            self._exit_status = info['task']
        elif info['warnflag'] == 1:
            self._exit_status = ('Number of simulations has '
                                 'exceeded defined limit: %d') % info['funcalls']
        else: self._optimized = True
        #save results
        self._optimum       = optimum
        self._named_optimum = self._name_parameters(self._optimum)
        #simulate PCR at optimum
        ovc, mixture = self._prepare_simulation(self._optimum)
        if ovc is None: return False
        self._PCR_Simulation = self._new_PCR_Simulation()
        self._PCR_Simulation.add_mixture(self._seq_name, mixture)
        self._PCR_Simulation.run(WorkCounter())
        self._have_results = bool(self._PCR_Simulation)
        return self._have_results
    #end def
    simulate_PCR = optimize_PCR_parameters


    def _format_header(self):
        header = ''
        if self._optimized:
            header += 'Optimization was successful.\n'
            header += wrap_text('Results of PCR simulation with optimized values of the parameters '
                                'are given below.\n\n')
        else:
            header += 'Optimization was unsuccessful: %s\n' % self._exit_status
            header += wrap_text('Results of PCR simulation with the values of the parameters '
                                'from the last iteration of the optimization are given below. These '
                                'values, nevertheless, usually give better simulation results than '
                                'the initial ones.\n\n')
        header += iPCR_Base._format_header(self)
        header += self._PCR_Simulation.format_report_header()
        return header
    #end def
    
    def _format_report_body(self):
        body  = ''
        body += self._PCR_Simulation.format_quantity_explanation()
        hit   = self._PCR_Simulation.hits()[0]
        body += self._PCR_Simulation.per_hit_header(hit)
        body += self._PCR_Simulation.per_hit_histogram(hit)
        body += '\n'
        body += self._PCR_Simulation.per_hit_electrophoresis(hit)
        return body
    #end def
#end class



if __name__ == '__main__':
    import sys, time
    from multiprocessing import Manager
    from DegenPrimer.Primer import Primer
    from SeqUtils import load_sequence
    from WorkCounter import WorkCounterManager
    from WaitingThread import WaitingThread
    from threading import Lock
    
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
                        args=('../ThGa.fa', 
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