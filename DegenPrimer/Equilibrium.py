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
# indicator_gddccontrol is distributed in the hope that it will be useful, but
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

from scipy.optimize import fmin_l_bfgs_b, fsolve
from numpy.random import random_sample

class Equilibrium(object):
    '''
    Calculate equilibrium parameters in a system of concurrent reactions.
    The system is defined as a dictionary of dicts: 
    {id: {'K': equilibrium_constant, 'Ai': Ai, 'Bi': Bi, 'type': reaction_type), ...}
    Ai and Bi are reactant identifiers. 
    All concentrations of Ai/Bi are equal to C_A/Bi initialization parameter
    reaction_type is a string identifier of one of the five types of the reaction:
    'AB'    A + B <--> AB
    '2A'    2A    <--> A_2
    '2B'    2B    <--> B_2
    'A'     A     <--> A'
    'B'     B     <--> B'
    '''

    def __init__(self, C_A, C_B, reactions, precision = 1e-10):
        #system in parameters
        self._precision = precision
        self._A         = C_A
        self._B         = C_B
        self._AB        = min(C_A, C_B)
        self._reactions = reactions
        self._max_K     = max(r['constant'] for r in self._reactions.values())
        #dictionaries for 'navigation' within the system
        ri = 0
        self._fwd_r_dict = dict()
        self._rev_r_dict = dict()
        self._Ai_reactions = dict()
        self._Bi_reactions = dict()
        for r_hash in self._reactions:
            #set up forward and reverse dictionaries between reaction hashes and numbers
            self._fwd_r_dict[r_hash] = ri
            self._rev_r_dict[ri] = r_hash
            #set up dictionaries of reactions in which a particular reactant participates
            R = self._reactions[r_hash]
            if not R['Ai'] in self._Ai_reactions:
                self._Ai_reactions[R['Ai']] = [ri]
            else: self._Ai_reactions[R['Ai']].append(ri)
            if not R['Bi'] in self._Bi_reactions:
                self._Bi_reactions[R['Bi']] = [ri]
            else: self._Bi_reactions[R['Bi']].append(ri)
            #next reaction
            ri += 1
        #solution
        self.solution = None
    #end def
    
    
    def __getitem__(self, r_hash):
        if self.solution == None:
            raise Exception('No solution has been found yet. Call calculate() method prior to solution lookup.')
        if not r_hash in self._fwd_r_dict:
            raise KeyError('There is no reaction in the system with the identifier "%s"' % str(r_hash))
        r_type = self._reactions[r_hash]['type']
        if   r_type == 'AB': return self._AB*self.solution[self._fwd_r_dict[r_hash]]
        elif r_type == '2A': return self._A *self.solution[self._fwd_r_dict[r_hash]]/2.0
        elif r_type == '2B': return self._B *self.solution[self._fwd_r_dict[r_hash]]/2.0
        elif r_type == 'A':  return self._A *self.solution[self._fwd_r_dict[r_hash]]
        elif r_type == 'B':  return self._B *self.solution[self._fwd_r_dict[r_hash]]
    #end def
    
    
    def _reactants_consumption(self, r, Ai, Bi):
        _A = 0
        for ri in self._Ai_reactions[Ai]:
            r_type = self._reactions[self._rev_r_dict[ri]]['type']
            if   r_type == 'AB': _A += self._AB*r[ri]
            elif r_type == '2A': _A += self._A *r[ri]
            elif r_type == 'A':  _A += self._A *r[ri]
        _B = 0
        for ri in self._Bi_reactions[Bi]:
            r_type = self._reactions[self._rev_r_dict[ri]]['type']
            if   r_type == 'AB': _B += self._AB*r[ri]
            elif r_type == '2B': _B += self._B *r[ri]
            elif r_type == 'B':  _B += self._B *r[ri]
        return _A, _B
    #end def
    
    
    #factory for left side functions of the system's equations
    def _function_factory(self, i):
        reaction = self._reactions[self._rev_r_dict[i]]
        if reaction['type'] == 'AB':
            def func(r):
                ABi = self._AB*r[i]
                _A, _B = self._reactants_consumption(r, reaction['Ai'], reaction['Bi'])
                return (ABi - ((self._A - _A)*(self._B - _B))*reaction['constant'])
            return func
        elif reaction['type'] == '2A':
            def func(r):
                Ai2 = self._A*r[i]/2.0
                _A, _B = self._reactants_consumption(r, reaction['Ai'], reaction['Bi'])
                return (Ai2 - ((self._A - _A)**2)*reaction['constant'])
            return func
        elif reaction['type'] == '2B':
            def func(r):
                Bi2 = self._B*r[i]/2.0
                _A, _B = self._reactants_consumption(r, reaction['Ai'], reaction['Bi'])
                return (Bi2 - ((self._B - _B)**2)*reaction['constant'])
            return func
        elif reaction['type'] == 'A':
            def func(r):
                A1i = self._A*r[i]
                _A, _B = self._reactants_consumption(r, reaction['Ai'], reaction['Bi'])
                return (A1i - (self._A - _A)*reaction['constant'])
            return func
        elif reaction['type'] == 'B':
            def func(r):
                B1i = self._B*r[i]
                _A, _B = self._reactants_consumption(r, reaction['Ai'], reaction['Bi'])
                return (B1i - (self._B - _B)*reaction['constant'])
            return func
    #end def

    def _initial_estimation(self, obj_func, r0, bounds):
        '''try to estimate the solution roughly with iterations of increasing precision'''
        factr = 1e15 #initial convergence precision
        while factr > 10:
            sol    = fmin_l_bfgs_b(obj_func, r0, factr = factr, bounds=bounds, approx_grad=True)
            objective_ratio = obj_func(r0)/obj_func(sol[0]) 
            if objective_ratio < 10:
                if objective_ratio > 1:
                    r0 = sol[0] 
                break
            factr *= 1e-5
            r0 = sol[0]
        return r0

    def calculate(self):
        #all left-side functions and initial root estimation
        l_funcs  = []
        r0       = [1e-6,]*len(self._reactions)
        for ri in range(len(self._reactions)):
            l_funcs.append(self._function_factory(ri))
        #system function for fsolve and objective function for fmin_ 
        sys_func = lambda(r): tuple(l_funcs[i](r) for i in range(len(l_funcs)))
        objective_function = lambda(r): sum((l_funcs[i](r))**2 for i in range(len(l_funcs)))
        objective_scale    = 1e12/sum((l_funcs[i](r0))**2 for i in range(len(l_funcs)))
        scaled_objective_function = lambda(r): sum((l_funcs[i](r))**2 for i in range(len(l_funcs)))*objective_scale
        #solution bounds
        bounds = ((0,1),)*len(l_funcs)
        #try to estimate the solution roughly with iterations of increasing precision
        r0 = self._initial_estimation(scaled_objective_function, r0, bounds)
        print '\ninitial estimation'
        print r0
        print sys_func(r0)
        print objective_function(r0)
        print ''
        #try to fine-tune the solution using fsolve
        sol = fsolve(sys_func, r0, full_output=True)
        print '\nfsolve'
        print sol[0]
        print sys_func(sol[0])
        print objective_function(sol[0])
        print ''
        #if failed for the first time, iterate fsolve while jitter r0 a little
        r0_value = objective_function(r0)
        while objective_function(sol[0]) >= r0_value \
        or    objective_function(sol[0]) >  self._precision \
        or    min(sol[0]) < 0 \
        or    max(sol[0]) > 1:
            sol = fsolve(sys_func, [ri+(random_sample(1)[0]-0.5)*1e-1 for ri in r0], xtol=1e-10, full_output=True)
            print sol[0]
            print sys_func(sol[0])
            print ''
        r0 = sol[0]
        #in any case use the best guess for the solution
        self.solution = r0
        self.solution_objective_value = objective_function(r0)
        return self.solution, self.solution_objective_value
    #end def
#end class


#tests
if __name__ == '__main__':
    import sys
    reactions = {'1': {'constant': 1e5,
                       'Ai': 1,
                       'Bi': 1,
                       'type': 'AB'},
                 '2': {'constant': 1e6,
                       'Ai': 1,
                       'Bi': 1,
                       'type': '2A'},
                 '3': {'constant': 1e9,
                       'Ai': 1,
                       'Bi': 1,
                       'type': '2B'},
                 '4': {'constant': 1e3,
                       'Ai': 1,
                       'Bi': 1,
                       'type': 'A'},
                 '5': {'constant': 1e2,
                       'Ai': 1,
                       'Bi': 1,
                       'type': 'B'},
                 }


    ov_sum = 0
    size = 1000
    scale = 1e-2/size
    for C_A, C_B in zip(random_sample(size), random_sample(size)):
        print 'A: %f, B: %f' % ((1+C_A*size)*scale, (1+C_B*size)*scale)
        eq_system = Equilibrium(C_A*1e-2, C_B, reactions)
        r0, obj = eq_system.calculate()
        ov_sum += obj
        print 'objective function value at the solution: %e' % obj
        print 'solution:', r0
        for R in reactions:
            print 'Reaction %s [%s]: %e' % (R, reactions[R]['type'], eq_system[R])
        print ''
    print 'OV_mean: %e' % (ov_sum/size)