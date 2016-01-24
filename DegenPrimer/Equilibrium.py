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

from BioUtils.Tools.AbortableBase import AbortableBase
from BioUtils.Tools.Multiprocessing import MultiprocessingBase
from numpy.random import random
from scipy.optimize import fsolve, newton_krylov
from scipy.sparse import csr_matrix, csgraph
import warnings

import numpy as np

class Reaction(object):
    '''
    Structure like object that defines a reversible mono or bimolecular chemical 
    reaction with four parameters:
    K      -- equilibrium constant of the reaction
    Ai, Bi -- identifiers of reactants
    r_type -- a string identifier of one of the following types of the reaction:
    'AB'    A + B <--> AB    bimolecular association
    '2A'    2A    <--> A_2   dimerization
    'A'     A     <--> A'    transformation
    '''
    
    __slots__ = ['constant', 'Ai', 'Bi', 'type']
    
    def __init__(self, 
                  K,     #equilibrium constant of the reaction
                  Ai,    #identifier of the first reactant
                  Bi,    #identifier of the second reactant, may be None
                  r_type #reaction type, one of the 'AB', '2A', 'A'
                 ):
        self.constant = K
        self.Ai       = Ai
        self.Bi       = Bi
        self.type     = r_type
    #end def
    
    def __str__(self):
        return ('[K=%f; Ai: %s; Bi: %s; type: %s]' %
                (self.constant, str(self.Ai), str(self.Bi), self.type))
    #end def
    
    def __repr__(self): return self.__str__()
#end class


class EquilibriumBase(object):
    '''Base class for equilibrium system solvers'''
    def __init__(self, reactions, concentrations, precision):
        self._precision      = precision
        #indexed reactants and concentrations
        self._reactants_ids  = concentrations.keys() #reactant index to ID 
        self._reactants_idx  = dict((r, i) for i, r in enumerate(self._reactants_ids)) #reactant ID to index
        #output values
        self.consumptions    = None #dict of reactants' consumptions
        self.solution        = None #dict of conversion degrees of reactions at equilibrium
        self.objective_value = -1 #value of the objective function of the system at equilibrium
    #end def
    
    def _indexed_reaction(self, R):
        return (R.constant, 
                self._reactants_idx[R.Ai], 
                None if R.Bi is None else self._reactants_idx[R.Bi], 
                R.type)
    #end def
#end class


class EquilibriumSolver(EquilibriumBase, AbortableBase):
    '''
    Calculate equilibrium parameters in a system of connected concurrent 
    reactions. The system is defined by dictionaries of Reactions and 
    concentrations of reactants.
    '''

    def __init__(self, abort_event, reactions, concentrations, precision):
        EquilibriumBase.__init__(self, reactions, concentrations, precision)
        AbortableBase.__init__(self, abort_event)
        #indexed reactants and concentrations
        self._concentrations = [concentrations[r] for r in self._reactants_ids] #list of reactants' concentrations
        self._min_C          = min(self._concentrations)
        #indexed reactions
        self._reaction_ids   = reactions.keys() #reaction index to ID
        self._reactions      = [self._indexed_reaction(reactions[rid])
                                 for rid in self._reaction_ids] #list of reactions
        self._Ri_reactions   = [[] for _n in self._concentrations]
        for ri, R in enumerate(self._reactions): #partition reactions by reactants
            self._Ri_reactions[R[1]].append((ri, R))
            if R[2] is not None:
                self._Ri_reactions[R[2]].append((ri, R))
    #end def
    
        
    def _reactants_consumption_factory(self, Ai):
        C_A = self._concentrations[Ai]
        reactions = self._Ri_reactions[Ai]
        n_reactions = len(reactions)
        C_max = np.zeros(n_reactions, dtype=float)
        ris = [0]*n_reactions
        for i, (ri, R) in enumerate(reactions):
            if R[3] == 'AB':
                r_Ai = R[1]
                if Ai == r_Ai:
                    C_B = self._concentrations[R[2]]
                else:
                    C_B = self._concentrations[r_Ai]
                C_max[i] = C_A if C_A < C_B else C_B
            else: 
                C_max[i] = C_A
            ris[i] = ri
        C_max = np.array(C_max)
        #consumption function for indexed reactant Ai
        def _reactant_consumption(r):
            return np.sum(C_max*r[ris])
        #end def
        return _reactant_consumption
    
    
    #factory for left side functions of the system's equations
    def _function_factory(self, i):
        reaction = self._reactions[i]
        r_type   = reaction[3]
        Ai       = reaction[1]
        C_A      = self._concentrations[Ai]
        K        = reaction[0]
        if r_type   == 'AB':
            Bi   = reaction[2]
            C_B  = self._concentrations[Bi]
            C_AB = C_A if C_A < C_B else C_B
            def func(r, consumptions):
                return (C_AB*r[i] - 
                        (C_A - consumptions[Ai])*
                        (C_B - consumptions[Bi])*K)
            return func
        elif r_type == '2A':
            C_A_2 = C_A/2.0
            def func(r, consumptions):
                C_A_C = (C_A - consumptions[Ai])
                return (C_A_2*r[i] - C_A_C*C_A_C*K)
            return func
        elif r_type == 'A':
            C_A_K = C_A*K
            def func(r, consumptions):
                return (C_A*r[i] - 
                        C_A_K + consumptions[Ai]*K)
            return func
    #end def
    
    
    @staticmethod
    def _fsolve(sys_func, r0, precision):
        return fsolve(sys_func, r0, xtol=precision)
    
    @staticmethod
    def _nksolve(sys_func, r0, precision):
        return newton_krylov(sys_func, r0, x_tol=precision)
        
    
    def _solve(self, sys_func, obj_func, r0, solver):
        #try to find the solution using selected solver
        sol = solver(sys_func, r0, self._precision)
        #if failed for the first time, iterate solver while jitter r0 a little
        r0_obj = obj_func(sol)
        while r0_obj > self._precision \
        or    np.min(sol) < 0 \
        or    np.max(sol) > 1:
            if self.aborted(): return None
            print '%d._solve.objf = %f > %f: %s' % (id(self), r0_obj, self._precision, r0_obj > self._precision)#test
            sol = solver(sys_func, r0 + (random(r0.shape)-r0)*0.3, self._precision)
            if  np.min(sol) >= 0 and np.max(sol) <= 1: 
                r0 = sol
                r0_obj = obj_func(sol)
        return sol
    #end def
    

    def calculate(self):
        #all left-side functions
        n_reactions = len(self._reactions)
        l_funcs  = [self._function_factory(ri) 
                    for ri in xrange(n_reactions)]
        #reactant consumption functions
        n_reactants = len(self._concentrations)
        rc_funcs = [self._reactants_consumption_factory(Ai) 
                    for Ai in xrange(n_reactants)]
        #initial evaluation point
        r0 = np.repeat(1e-6, n_reactions)
        #system function and objective function for solver
        def sys_func(r):
            consumptions = np.zeros(n_reactants, dtype=float)
            for i in xrange(n_reactants): consumptions[i] = rc_funcs[i](r)
            new_r = np.zeros(n_reactions, dtype=float)
            for i in xrange(n_reactions): new_r[i] = l_funcs[i](r, consumptions)
            return new_r/self._min_C 
        def obj_func(r):
            new_r = sys_func(r)
            return np.sum(new_r*new_r)
        #choose the solver
        if n_reactions < 100: 
            solver = self._fsolve #fsole is fast for small systems
            altsolv = self._nksolve
        else: 
            solver = self._nksolve #but nksolve is MUCH faster for big ones
            altsolv = self._fsolve
        #solve the system
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try: sol = self._solve(sys_func, obj_func, r0, solver)
            except: 
                print '\nFailed to solve with %s, trying with %s' % (solver, altsolv)#test
                try: sol = self._solve(sys_func, obj_func, r0, altsolv)
                except Exception, e:
                    print '\nUnable to calculate equilibrium.'
                    print e
                    sol = None
        if sol is None: return False
        #calculate solution objective function and reactant consumptions
        self.objective_value = obj_func(sol)
        self.solution     = dict((self._reaction_ids[ri], r) 
                                 for ri, r in enumerate(sol))
        self.consumptions = dict((self._reactants_ids[ri], func(sol)) 
                                 for ri, func in enumerate(rc_funcs))
        return True
    #end def
#end class


class Equilibrium(EquilibriumBase, MultiprocessingBase):
    '''
    Calculate equilibrium parameters in arbitrary system of concurrent reactions.
    The system is defined as a dictionary of Reaction objects: 
    {reaction_ID: Reaction}.
    Initial concentrations of reactants are provided as a dictionary:
    {reactant_ID: concentration}
    '''

    def __init__(self, abort_event, reactions, concentrations, precision = 1e-10):
        MultiprocessingBase.__init__(self, abort_event)
        EquilibriumBase.__init__(self, reactions, concentrations, precision)
        #input parameters
        self.reactions         = reactions
        self.concentrations    = concentrations
        #group reactions by their connected graph components 
        self._reactions_groups = self._group_reactions()
    #end def
    
    def __nonzero__(self): return self.solution is not None
    
    def _index_reaction(self, reaction):
        rid, R = reaction
        iR = self._indexed_reaction(R)
        Ai = iR[1]; Bi = iR[2] if iR[2] is not None else Ai
        return rid, iR, (Ai,Bi)
    #end def
    
    def _group_reactions(self):
        n_reactants = len(self.concentrations)
        #convert reactions to indexed form and fill reactions graph and
        ireactions = self.parallelize_work(1, self._index_reaction, 
                                           self.reactions.items())
        if ireactions is None: return []
        inds = np.array([AiBi for _rid, _iR, AiBi in ireactions]).transpose()
        data = np.ones(len(ireactions))
        r_graph = csr_matrix((data, (inds[0],inds[1])), shape=(n_reactants, n_reactants), dtype='int8')
        #lable reactants by connected component of the reactions graph
        n_comps, comps = csgraph.connected_components(r_graph, directed=False)
        del r_graph
        #group reactions and concentrations by connected component of the reactions graph
        if n_comps < 2: return []
        groups = [[dict(), dict()] for _n in xrange(n_comps)]
        for rid, iR, _AiBi in ireactions: 
            groups[comps[iR[1]]][0][rid] = self.reactions[rid]
        for cidx, cid in enumerate(self._reactants_ids):
            groups[comps[cidx]][1][cid]  = self.concentrations[cid]
        #remove empty groups
        groups = [group for group in groups if group[0]]
        return groups
    #end def
    
    
    def get_conversion_degree(self, r_hash):
        if self.solution == None:
            raise Exception('No solution has been found yet. Call calculate() method prior to solution lookup.')
        if r_hash not in self.solution:
            raise KeyError('There is no reaction in the system with the identifier "%s"' % str(r_hash))
        return self.solution[r_hash]
    #end def
    
    
    def get_product_concentration(self, r_hash):
        if self.solution == None:
            raise Exception('No solution has been found yet. Call calculate() method prior to solution lookup.')
        if r_hash not in self.solution:
            raise KeyError('There is no reaction in the system with the identifier "%s"' % str(r_hash))
        reaction = self.reactions[r_hash]
        r_type   = reaction.type
        C_A      = self.concentrations[reaction.Ai]
        r        = self.solution[r_hash]
        if   r_type == 'AB':
            C_B  = self.concentrations[reaction.Bi]
            C_AB = C_A if C_A < C_B else C_B
            return C_AB*r
        elif r_type == '2A': return C_A/2.0*r
        elif r_type == 'A' : return C_A*r
    #end def
    
    
    def reactants_consumption(self, Ai, Bi=None):
        if self.solution == None:
            raise Exception('No solution has been found yet. Call calculate() method prior to solution lookup.')
        if Ai not in self.consumptions:
            return 0, 0
        if Bi is None or Bi not in self.consumptions:
            return self.consumptions[Ai],None
        return self.consumptions[Ai],self.consumptions[Bi]
    #end def
    

    def _calculate(self, reactions_group, precision):
        eq = EquilibriumSolver(self._abort_event, 
                               reactions_group[0], reactions_group[1], precision)
        if eq.calculate():
            return eq.objective_value, eq.solution, eq.consumptions
        return -1, None, None
    #end def

    def calculate(self, counter):
        result = False
        if len(self._reactions_groups) < 2:
            #single process mode
            equilibrium = EquilibriumSolver(self._abort_event, 
                                            self.reactions, 
                                            self.concentrations, 
                                            self._precision)
            if equilibrium.calculate():
                self.objective_value = equilibrium.objective_value
                self.solution        = equilibrium.solution
                self.consumptions    = equilibrium.consumptions
                result = True
            counter.done()
        else:
            #parallel processing
            solutions = self.parallelize_work(1, self._calculate, 
                                               self._reactions_groups, 
                                               self._precision, counter=counter)
            if solutions is None: return
            #parse solutions
            self.objective_value = -1
            self.solution        = dict()
            self.consumptions    = dict()
            for obj_val, sol, cons in solutions:
                if sol is None:
                    self.objective_value = -1
                    self.solution        = None
                    break
                if obj_val > self.objective_value: 
                    self.objective_value = obj_val
                self.solution.update(sol)
                self.consumptions.update(cons)
            result = self.solution is not None
        return result
    #end def
#end class