# coding=utf-8

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
Created on Mar 13, 2014

@author: Allis Tauri <allista@gmail.com>
'''

from math import sqrt
from multiprocessing import Value
from textwrap import dedent


class PCR_P_Meta(type):
    def __new__(cls, name, bases, attrs):
        for prop in attrs.get('parameters', ()):
            attrs['_'+prop] = Value('d', 0)
            pdict = dict(); exec dedent('''
            def _get_%(prop)s(self): return self._%(prop)s.value
            def _set_%(prop)s(self, val): self._%(prop)s.value = val
            %(prop)s = property(_get_%(prop)s, _set_%(prop)s)
            ''') % dict(prop=prop) in pdict
            attrs[prop] = pdict[prop]
        return super(PCR_P_Meta, cls).__new__(cls, name, bases, attrs)
    #end def
#end metaclass

class PCR_Parameters_Base(object):
    def _recalculate(self): pass
    
    parameters = []        
    
    def set(self, params):
        for name in params:
            if name in self.parameters:
                setattr(self, name, params[name])
        self._recalculate() 
    #end def
    
    def __str__(self):
        if not self.parameters: return ''
        s  = '\n'
        max_par_len = max(len(str(par)) for par in self.parameters)
        for name in self.parameters:
            s += '%s%s: %f\n' % (name, ' '*(max_par_len-len(name)), 
                                 getattr(self, name)) 
        return s
    #end def
#end class

class PCR_Parameters_Copy(PCR_Parameters_Base): pass
class PCR_Parameters(PCR_Parameters_Base):
    __metaclass__ = PCR_P_Meta
    parameters    = ['Mg', 'Na', 'dNTP', 'DNA', 'DMSO', 'PCR_T', 'Na_eq']
    
    def __init__(self):
        #standard PCR conditions
        self._Mg.value    = 1.5e-3 #M
        self._Na.value    = 50e-3  #M
        self._dNTP.value  = 0.1e-3 #M
        self._DNA.value   = 5e-9   #M; DNA template concentration
        self._DMSO.value  = 0.0    #%; v/v-percent concentration of DMSO
        self._PCR_T.value = 60.0   #C; temperature at which PCR is conducted
        self._recalculate()
    #end def
    
    def _recalculate(self):
        ''''divalent cation correction:
        all concentrations should be in mM (Ahsen et al., 2001)'''
        self._Na_eq.value = (self._Na.value + 
                             sqrt((self._Mg.value - self._dNTP.value)*14.4))
    #end def
    
    def __deepcopy__(self, memo):
        params = dict(parameters=self.parameters)
        for name in self.parameters:
            param = getattr(self, '_'+name)
            params[name] = param.value
        _copy = PCR_Parameters_Copy()
        _copy.__dict__.update(params)
        return _copy
    #end def
#end class