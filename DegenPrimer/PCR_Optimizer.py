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

import TD_Functions
from scipy.optimize import fmin_l_bfgs_b
from StringTools import wrap_text
from iPCR_Base import iPCR_Base
from MultiprocessingBase import MultiprocessingBase
from AllSecStructures import AllSecStructures


class PCR_Optimizer(MultiprocessingBase, iPCR_Base):
    '''
    Repeatedly perform PCR simulation varying parameters 
    to find optimal conditions
    '''

    _PCR_report_suffix = 'PCR-optimization'
    
    def __init__(self, abort_event, max_steps, purity, *args, **kwargs):
        MultiprocessingBase.__init__(self, abort_event)
        iPCR_Base.__init__(self, *args, **kwargs)
        #simulators
        self._sec_structures     = None
        self._PCR_Simulations    = []
        #target
        self._product_bounds     = None
        self._seq_name           = None
        #parameters
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
    
    
    @staticmethod
    def _recalculate_annealing(annealing):
        for duplex, _unused in annealing[1]:
            duplex.recalculate()
        return annealing
    #end def
    
    
    def _recalculate_annealings(self):
        for i, annealings in self._seq_annealings.items():
            if self._abort_event.is_set(): return
            fwd_annealings = self.parallelize_work(1, 
                                                   self._recalculate_annealing, 
                                                   annealings[0])
            rev_annealings = self.parallelize_work(1, 
                                                   self._recalculate_annealing, 
                                                   annealings[1])
            self._seq_annealings[i] = (fwd_annealings, rev_annealings)
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
        TD_Functions.PCR_P.set(named_parameters)
        #set primer concentrations
        for primer in self._primers:
            if primer.id in named_parameters:
                primer.total_concentration = named_parameters[primer.id]
        #change polymerase concentration
        if 'polymerase' in named_parameters:
            self._polymerase = named_parameters['polymerase']
    #end def
    
    
    def _objective_function(self, parameters):
        result = 1
        self._set_parameters(parameters)
        #recalculate annealings' dG and K
        self._recalculate_annealings()
        #recalculate secondary structures
        sec_structs = AllSecStructures(self._abort_event, self._job_id, self._primers)
        sec_structs.find_structures()
        self._side_concentrations = sec_structs.concentrations()
        self._side_reactions      = sec_structs.reactions()
        #create new simulator
        _PCR_Sim = self._new_PCR_Simulation()
        self._PCR_Simulations.append(_PCR_Sim)
        if not self._add_annealings_for_seqs(_PCR_Sim): return result
        #simulate PCR
        _PCR_Sim.run()
        if not _PCR_Sim:
            raise RuntimeError('\nPCR Optimizer: unable to simulate PCR during '
                               'optimization\n')
        result = 0
        for bounds in self._product_bounds:
            prod_C = _PCR_Sim.product_concentration(self._seq_name, 
                                                    (bounds['start'],bounds['end']))
            result += -1e9*prod_C[0]*(prod_C[1]**self._purity)
        self._PCR_Simulations.remove(_PCR_Sim)
        del _PCR_Sim
        return result
    #end def


    def optimize_PCR_parameters(self, seq_file, product_bounds, 
                                   parameters, seq_id=None):
        #reset optimized flag
        self._optimized = False
        if not product_bounds or not parameters: return False
        #prepare simulation parameters
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
        #find possible products
        if not self._find_annealings((seq_file,), 
                                     (seq_id,) if seq_id else None): 
            return False
        self._seq_name       = self._seq_names.values()[0]
        #optimize parameters
        print '\n'
        optimum, _unused, info = fmin_l_bfgs_b(self._objective_function, 
                                               self._iparams, approx_grad=True, 
                                               bounds=self._ibounds, factr=1e12,
                                               maxfun=self._max_steps, disp=1)
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
        self._set_parameters(self._optimum)
        self._recalculate_annealings()
        sec_structs = AllSecStructures(self._abort_event, self._job_id, self._primers)
        sec_structs.find_structures()
        self._side_concentrations = sec_structs.concentrations()
        self._side_reactions      = sec_structs.reactions()
        self._PCR_Simulation      = self._new_PCR_Simulation()
        if self._add_annealings_for_seqs(self._PCR_Simulation):
            self._PCR_Simulation.run()
        self._have_results = bool(self._PCR_Simulation)
        return self._have_results
    #end def


    def _format_header(self):
        header = ''
        if self._optimized:
            header += 'Optimization was successful.\n'
            header += ('Results of PCR simulation with optimized values of the parameters '
                       'are given below.\n\n')
        else:
            header += 'Optimization was unsuccessful: %s\n' % self._exit_status
            header += ('Results of PCR simulation with the values of the parameters '
                       'from the last iteration of the optimization are given below. These '
                       'values, nevertheless, usually give better simulation results than '
                       'the initial ones.\n\n')
        header  = wrap_text(header)
        header += self._PCR_Simulation.format_report_header()
        return header
    #end def
    
    
    def _format_report_body(self):
        body = ''
        body += self._PCR_Simulation.format_report_header()
        body += self._PCR_Simulation.format_quantity_explanation()
        hit = self._PCR_Simulation.hits()[0]
        body += self._PCR_Simulation.per_hit_header(hit)
        body += self._PCR_Simulation.per_hit_histogram(hit)
        body += '\n'
        body += self._PCR_Simulation.per_hit_electrophoresis(hit)
        return body
    #end def
#end class



if __name__ == '__main__':
    from Primer import Primer, load_sequence
    import os, cProfile
    os.chdir('../')
    TD_Functions.PCR_P.PCR_T = 53
    TD_Functions.PCR_P.Mg    = 3e-3
    TD_Functions.PCR_P.dNTP  = 200e-6
    TD_Functions.PCR_P.DNA   = 1e-10
    fwd_primer = Primer(load_sequence('ATATTCTACRACGGCTATCC', 'fwd_test', 'fwd_test'), 0.5e-6)
    rev_primer = Primer(load_sequence('GAASGCRAAKATYGGGAAC', 'rev_test', 'rev_test'), 0.5e-6)
    optimizer = PCR_Optimizer(100, 
                              5, 'test-job', 
                              [fwd_primer,
                               rev_primer],
                              50, 
                              1500, 
                              100000,
                              False, 
                              33,
                              include_side_annealings=False)
    optimizer.optimize_PCR_parameters('ThGa.fa', ({'start':47920, 'end':49321},), 
                                      (
                                       {'name':'PCR_T',
                                        'min':50, 'ini':60, 'max':72},
                                       {'name':'dNTP',
                                        'min':100e-6, 'ini':200e-6, 'max':900e-6},
                                       ),
                                      )
    optimizer.write_report()
    pass