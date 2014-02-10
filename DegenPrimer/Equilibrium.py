#!/usr/bin/python
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
Created on Nov 25, 2012

@author: Allis Tauri <allista@gmail.com>
'''
from scipy.optimize import fsolve
from numpy.random import random_sample

class Equilibrium(object):
    '''
    Calculate equilibrium parameters in a system of concurrent reactions.
    The system is defined as a dictionary of dicts: 
    {id: {'K':    equilibrium_constant, 
          'Ai':   first reactant_id, 
          'Bi':   second reactant_id, 
          'type': reaction_type), ...}
    Initial concentrations of reactants are provided as a dict:
    {reactant_id: concentration}
    reaction_type is a string identifier of one of the following types of the reaction:
    'AB'    A + B <--> AB    bimolecular association
    '2A'    2A    <--> A_2   dimerization
    'A'     A     <--> A'    transformation
    '''
    
    @staticmethod
    def compose_reaction(K,     #equilibrium constant of the reaction
                           Ai,    #identifier of the first reactant
                           Bi,    #identifier of the second reactant, if present
                           r_type #reaction type, one of the 'AB', '2A', 'A'
                           ):
        return {'constant'  : K,
                'Ai'        : Ai,
                'Bi'        : Bi,
                'type'      : r_type}
    #end def
    
    
    @classmethod
    def indexed_reaction(cls, reaction):
        return [reaction['constant'],
                reaction['Ai'],
                reaction['Bi'],
                reaction['type']]
    #end def
    

    def __init__(self, reactions, concentrations, precision = 1e-10):
        #system in parameters
        self._precision     = precision
        self.reactions      = reactions
        self.concentrations = concentrations
        #dictionaries for 'navigation' within the system
        reaction_i = 0 #reaction index
        reactant_i = 0 #reactant index
        self._rev_r_dict        = dict() #reaction hash by index
        self._Ri_reactions      = list() #indexed list of reactions in which reactant Ri participates
        self._ireactions        = list() #indexed list of reactions
        self._iconcentrations   = list() #indexed list of concentrations
        self._reactants_dict    = dict() #dict of reactant indices
        self._min_concentration = min(self.concentrations.values())
        for r_hash in self.reactions:
            R  = self.reactions[r_hash]
            #check reactants
            if R['Ai'] != None and R['Ai'] not in self.concentrations:
                raise ValueError(('Equilibrium: reactant %s is not present in '
                                  'the list of concentrations.') % str(R['Ai']))
            if R['Bi'] != None and R['Bi'] not in self.concentrations:
                raise ValueError(('Equilibrium: reactant %s is not present in '
                                  'the list of concentrations.') % str(R['Bi']))
            #add reaction to the list and add it's hash to the reverse dicts
            Ri = self.indexed_reaction(R)
            self._ireactions.append(Ri)
            self._rev_r_dict[reaction_i] = r_hash
            #index reactants
            if R['Ai'] not in self._reactants_dict:
                self._reactants_dict[R['Ai']] = reactant_i
                self._iconcentrations.append(self.concentrations[R['Ai']])
                self._Ri_reactions.append([])
                reactant_i += 1
            Ai = self._reactants_dict[R['Ai']]
            Ri[1] = Ai
            if R['Bi'] != None:
                if R['Bi'] not in self._reactants_dict:
                    self._reactants_dict[R['Bi']] = reactant_i
                    self._iconcentrations.append(self.concentrations[R['Bi']])
                    self._Ri_reactions.append([])
                    reactant_i += 1
                Bi = self._reactants_dict[R['Bi']]
            else: Bi = None
            Ri[2] = Bi
            #add reaction to the list of reactant's reactions
            self._Ri_reactions[Ai].append(reaction_i)
            if Bi != None: self._Ri_reactions[Bi].append(reaction_i)
            #next reaction
            reaction_i += 1
        #solution:
        #indexed list of reactants consumption
        self._consumptions = None
        self._solution     = None #pure solver output
        self.solution      = None #solution mapped to reaction hashes
        self.solution_objective_value = -1
    #end def
    
    
    def get_conversion_degree(self, r_hash):
        if self.solution == None:
            raise Exception('No solution has been found yet. Call calculate() method prior to solution lookup.')
        if not r_hash in self.reactions:
            raise KeyError('There is no reaction in the system with the identifier "%s"' % str(r_hash))
        return self.solution[r_hash]
    #end def
    
    
    def get_product_concentration(self, r_hash):
        if self.solution == None:
            raise Exception('No solution has been found yet. Call calculate() method prior to solution lookup.')
        if not r_hash in self.reactions:
            raise KeyError('There is no reaction in the system with the identifier "%s"' % str(r_hash))
        reaction = self.reactions[r_hash]
        r_type   = reaction['type']
        C_A      = self.concentrations[reaction['Ai']]
        r        = self.solution[r_hash]
        if   r_type == 'AB':
            C_B  = self.concentrations[reaction['Bi']]
            C_AB = C_A if C_A < C_B else C_B
            return C_AB*r
        elif r_type == '2A': return C_A/2.0*r
        elif r_type == 'A' : return C_A*r
    #end def
    
    
    def reactants_consumption(self, Ai, Bi=None):
        if self.solution == None:
            raise Exception('No solution has been found yet. Call calculate() method prior to solution lookup.')
        if not Ai in self._reactants_dict:
            return 0, 0
        if Bi is None or not Bi in self._reactants_dict:
            return (self._consumptions[self._reactants_dict[Ai]], None)
        return (self._consumptions[self._reactants_dict[Ai]], 
                self._consumptions[self._reactants_dict[Bi]])
    #end def
    
    
    def __getitem__(self, r_hash):
        if self.solution == None:
            raise Exception('No solution has been found yet. Call calculate() method prior to solution lookup.')
        if not r_hash in self.reactions:
            raise KeyError('There is no reaction in the system with the identifier "%s"' % str(r_hash))
        r_type = self.reactions[r_hash]['type']
        if   r_type == 'AB': return min(self.concentrations[self.reactions[r_hash]['Ai']],
                                        self.concentrations[self.reactions[r_hash]['Bi']])*self.solution[r_hash]
        elif r_type == '2A': return self.concentrations[self.reactions[r_hash]['Ai']] *self.solution[r_hash]/2.0
        elif r_type == 'A':  return self.concentrations[self.reactions[r_hash]['Ai']] *self.solution[r_hash]
    #end def
    
        
    def _reactants_consumption_factory(self, Ai):
        C_A = self._iconcentrations[Ai]
        _A_list = []
        for ri in self._Ri_reactions[Ai]:
            reaction = self._ireactions[ri]
            r_type   = reaction[3]
            if   r_type == 'AB':
                r_Ai = reaction[1]
                r_Bi = reaction[2]
                if Ai == r_Ai:
                    C_B = self._iconcentrations[r_Bi]
                else:
                    C_B = self._iconcentrations[r_Ai]
                C_AB = C_A if C_A < C_B else C_B
                _A_list.append((C_AB, ri))
            else: _A_list.append((C_A, ri)) 
        #consumption function for indexed reactant Ai
        def _reactant_consumption(r):
            _A = 0
            for _C_max, ri in _A_list:
                _A += _C_max*r[ri]
            return _A
        #end def
        return _reactant_consumption
    
    
    #factory for left side functions of the system's equations
    def _function_factory(self, i):
        reaction = self._ireactions[i]
        r_type   = reaction[3]
        Ai       = reaction[1]
        C_A      = self._iconcentrations[Ai]
        K        = reaction[0]
        if r_type   == 'AB':
            Bi   = reaction[2]
            C_B  = self._iconcentrations[Bi]
            C_AB = C_A if C_A < C_B else C_B
            def func(r, consumptions):
                return (C_AB*r[i] - 
                        (C_A - consumptions[Ai])*
                        (C_B - consumptions[Bi])*K)/self._min_concentration
            return func
        elif r_type == '2A':
            C_A_2 = C_A/2.0
            def func(r, consumptions):
                return (C_A_2*r[i] - 
                        ((C_A - consumptions[Ai])**2)*K)/self._min_concentration
            return func
        elif r_type == 'A':
            C_A_K = C_A*K
            def func(r, consumptions):
                return (C_A*r[i] - 
                        C_A_K + consumptions[Ai]*K)/self._min_concentration
            return func
    #end def
    

    def calculate(self):
        #all left-side functions
        l_funcs  = [self._function_factory(ri) 
                    for ri in xrange(len(self.reactions))]
        #reactant consumption functions
        rc_funcs = [self._reactants_consumption_factory(Ai) 
                     for Ai in xrange(len(self._iconcentrations))]
        #initial evaluation point
        r0 = [1e-6,]*len(self.reactions)
        #system function for fsolve and objective function for fmin_
        def sys_func(r):
            consumptions = [func(r) for func in rc_funcs]
            return tuple(l_func(r, consumptions) for l_func in l_funcs)
        objective_function = lambda(r): sum(sf**2 for sf in sys_func(r))
        #try to find the solution using fsolve
        sol = fsolve(sys_func, r0, full_output=True)
        #if failed for the first time, iterate fsolve while jitter r0 a little
        r0_obj = objective_function(sol[0])
        while r0_obj > self._precision \
        or    min(sol[0]) < 0 \
        or    max(sol[0]) > 1:
            sol = fsolve(sys_func, [ri+(random_sample(1)[0]-ri)*0.3 for ri in r0], 
                         xtol=1e-10, full_output=True)
            if  min(sol[0]) >= 0 and max(sol[0]) <= 1: 
                r0 = sol[0]
                r0_obj = objective_function(sol[0])
        r0 = sol[0]
        #in any case use the best guess for the solution
        self._consumptions = [func(r0) for func in rc_funcs]
        self._solution     = r0
        self.solution      = dict()
        for ri in xrange(len(r0)):
            self.solution[self._rev_r_dict[ri]] = r0[ri]
        self.solution_objective_value = objective_function(r0)
        return self.solution, self.solution_objective_value
    #end def
#end class


#tests
if __name__ == '__main__':
    from scipy.stats import sem, histogram
    import numpy as np
    reactions = {'A1A2': {'constant': 1e15,
                       'Ai': 1,
                       'Bi': 2,
                       'type': 'AB'},
                 '2A1': {'constant': 1e6,
                       'Ai': 1,
                       'Bi': None,
                       'type': '2A'},
                 '2A2': {'constant': 1e9,
                       'Ai': 2,
                       'Bi': None,
                       'type': '2A'},
                 'A1*': {'constant': 1e3,
                       'Ai': 1,
                       'Bi': None,
                       'type': 'A'},
                 'A2*': {'constant': 1e2,
                       'Ai': 2,
                       'Bi': None,
                       'type': 'A'},
                 }


    import cProfile, os
    os.chdir('../')
    ov_sum = 0
    size = 1000
    scale = 1e-2/size
    
    C_A1, C_A2 = 0.6543, 0.0353
    C_dict = {1: (1+C_A1*size)*scale, 2: (1+C_A2*size)*scale}
    solutions = []
#    pow = lambda x: x**2
#    cProfile.run('''
#for i in xrange(1000):
#    abs(-123.464)
#    x = pow(-123.464)''',
#'abs-power.profile')
    cProfile.run('''
for i in xrange(size):
    eq_system = Equilibrium(reactions, C_dict, 1e-10)
    solutions.append(eq_system.calculate())
''',
'equilibrium.profile')
    sols = dict((key, []) for key in solutions[0][0].keys())
    for sol, obj in solutions:
        for s in sol:
            sols[s].append(sol[s])
    for s in sols:
        print s, np.mean(sols[s]), sem(sols[s])
        print histogram(sols[s])
        print '\n'
#    for C_A1, C_A2 in zip(random_sample(size), random_sample(size)):
#        print 'A1: %f, A2: %f' % ((1+C_A1*size)*scale, (1+C_A2*size)*scale)
#        C_dict = {1: (1+C_A1*size)*scale, 2: (1+C_A2*size)*scale}
#        eq_system = Equilibrium(reactions, C_dict)
#        sol, obj = eq_system.calculate()
#        ov_sum += obj
#        print 'objective function value at the solution: %e' % obj
#        print 'solution:', sol
#        for R in reactions:
#            print 'Reaction %s [%s]: %e' % (R, reactions[R]['type'], eq_system[R])
#        print ''

    print 'Done'