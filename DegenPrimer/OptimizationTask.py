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
Created on Mar 27, 2013

@author: Allis Tauri <allista@gmail.com>
'''


from PrimerTaskBase import PrimerTaskBase
from PCR_Optimizer import PCR_Optimizer
import TD_Functions


class OptimizationTask(PrimerTaskBase):
    
    def __init__(self):
        PrimerTaskBase.__init__(self)
        self._optimizer          = None
        self._primers            = []
        self._valid_parameters   = []
    #end def


    @staticmethod
    def check_options(args):
        if not args.optimize: return False
        if not PrimerTaskBase.check_options(args): return False
        has_targets = False
        if args.target_product_list:
            for product in args.target_product_list:
                if product['start'] < product['end']:
                    has_targets = True
                    break
        if not has_targets: return False
        has_parameters = False
        if args.optimization_parameter_list:
            for par in args.optimization_parameter_list:
                if par['enabled']:
                    has_parameters = True
                    break
        if not has_parameters: return False
        return (bool(args.fasta_files) or 
                bool(args.sequence_db) and bool(args.use_sequences))
    #end def                
        
    
    def terminate(self):
        if self._optimizer is None: return
        del self._optimizer
        self._optimizer = None
    #end def
        
    
    def run(self, args):
        self._primers = args.primers
        self._define_valid_parameters()
        #set PCR parameters
        TD_Functions.PCR_P.set(args.options)
        #parse optimization parameters
        parameters = args.optimization_parameter_list
        parameters = [par for par in parameters 
                      if par['enabled'] and 
                      self._parse_parameter_name(par)]
        for par in parameters: self._scale_parameter(par, args)
        #perform optimization
        self._optimizer = PCR_Optimizer(args.max_simulations,
                                        args.product_purity,
                                        args.max_mismatches,
                                        args.job_id, 
                                        self._primers, 
                                        args.min_amplicon, 
                                        args.max_amplicon, 
                                        args.polymerase, 
                                        args.with_exonuclease, 
                                        args.cycles,
                                        None, 
                                        None,
                                        args.analyse_all_annealings)
        seq_file = args.sequence_db or args.fasta_files[0]
        seq_id   = args.use_sequences[0] if args.use_sequences else None
        if self._optimizer.optimize_PCR_parameters(seq_file, 
                                                   args.target_product_list,
                                                   parameters,
                                                   seq_id):
            self._optimizer.write_report()
            for report in self._optimizer.reports():
                args.register_report(**report)
        del self._optimizer
        self._optimizer = None
        return 1
    #end def


    def _define_valid_parameters(self):
        self._valid_parameters.extend(primer.id for primer in self._primers)
        self._valid_parameters.extend(vars(TD_Functions.PCR_P))
        self._valid_parameters.extend(['polymerase'])
    #end def
        
    
    def _parse_parameter_name(self, parameter):
        name = parameter['name']
        if  name in self._valid_parameters:
            return True
        _name = name.replace('-', '_').strip('_') 
        if _name in self._valid_parameters: 
            parameter['name'] = _name
            return True
        _name = name.replace(' ', '_').strip('_')
        if _name in self._valid_parameters: 
            parameter['name'] = _name
            return True
        print '\nOptimizationTask: invalid optimization parameter name: %s' % name
        return False
    #end def
        
    
    def _scale_parameter(self, parameter, args):
        name   = parameter['name']
        option = None
        scale  = 1
        if name in [p.id for p in self._primers]:
            option = args.find_option('primers/primer/concentration')
        elif name == 'polymerase':
            option = args.find_option('iPCR/polymerase')
        else:
            option = args.find_option('PCR/'+name)
        if option is not None and option.scale:
            scale = option.scale
        parameter['min'] = scale*float(parameter['min'])
        parameter['ini'] = scale*float(parameter['ini'])
        parameter['max'] = scale*float(parameter['max'])
    #end def
#end class