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
Created on Mar 22, 2013

@author: Allis Tauri <allista@gmail.com>
'''

from argparse import Action
from math import ceil, log

from .StringTools import wrap_text


class OptionBase(object):
    def __init__(self, name, desc):
        self._name     = name
        self._desc     = desc
        self._parent   = None
        self.options   = []
    #end def
    
    
    def fill_optional_attributes(self, **kwargs):
        for _name, _value in kwargs.items():
            if hasattr(self, _name): setattr(self, _name, _value)
    #end def
    
    @property
    def name(self): return self._name
    
    @property
    def full_name(self):
        if self._parent is None: return self._name
        return self._parent.full_name + '/' + self._name
    #end def
    
    
    @property
    def desc(self): return self._desc
    
    @property
    def formatted_desc(self):
        return wrap_text(self.desc.replace('%%', '%')).strip()
    
    @property
    def group(self): return self._name
    
    @property
    def is_compound(self): return bool(self.options)
#end class


class OptionGroup(OptionBase):
    def __init__(self, name, desc, **kwargs):
        OptionBase.__init__(self, name, desc)
        self.is_mutex = False
        #get options
        self.options = kwargs.pop('options', None)
        if not self.options: 
            raise ValueError('OptionGroup should be provided with options.')
        #fill optional attributes
        self.fill_optional_attributes(**kwargs)
        #set parent for options and force keyword attributes
        for opt in self.options: 
            opt.set_parent(self)
            opt.fill_optional_attributes(**kwargs)
    #end def
    
    @property
    def is_compound(self): return True
#end class


class Option(OptionBase):
    MULTY_NARGS = ('+', '*')
    CMP_NARGS   = (1, '?', '+', '*')
    
    def __init__(self, name, desc, nargs, py_type, **kwargs):
        OptionBase.__init__(self, name, desc)
        self._parent    = None
        #set mandatory fields
        self._nargs     = nargs
        self.py_type    = py_type
        #define optional fields
        self.units      = None
        self.default    = None
        self.limits     = None
        self.scale      = None
        self.save       = True
        #for GUI
        self.field_type = py_type.__name__ if py_type is not None else None
        #for argparse
        self.args       = ('--%s' % name.replace('_', '-'),)
        self.metavar    = self.field_type        
        self.required   = False
        self.value_required = True
        #fill optional fields_list
        self.fill_optional_attributes(**kwargs)
        #set proper nargs and metavar
        if self.options:
            self.metavar = 'val'
            for opt in self.options: opt.set_parent(self)
            if self._nargs not in self.CMP_NARGS:
                raise ValueError(('\nOption: %s is a compound option;' 
                                  'nargs should be one of %s') % str(self.CMP_NARGS))
        #set proper default
        if self.is_multi or self.is_poly:
            self.default = []
        elif self.py_type is str:
            self.default = ''
    #end def
    
    
    def set_parent(self, parent):
        self._parent = parent
        if isinstance(self._parent, Option) and self._nargs in self.MULTY_NARGS:
            raise ValueError('Option: suboptions must not have * or + nargs.')
    #end def
    

    @property
    def fields_list(self):
        fields_list = []
        if self.options:
            for opt in self.options:
                fields_list.extend(opt.fields_list)
        else: fields_list.append((self._name, self.py_type))
        return fields_list
    #end def
    
    
    @property
    def dest(self):
        if self.is_poly: return self._name+'_list'
        else: return self._name
    
    
    @property
    def group(self):
        if self._parent is None: return self._name
        return self._parent.group
            
    
    @property
    def nargs(self):
        if self.options:
            return sum(opt.nargs if type(opt.nargs) is int else 1 
                       for opt in self.options)
        else: return self._nargs
    #end def
 
        
    @property
    def is_multi(self): #if a non-compound option takes multiple values: --foo val1 val2...
        if self.options: return False
        if self._nargs in self.MULTY_NARGS: return True
        if type(self._nargs) is int and self._nargs > 1: return True
        return False
    #end def
    
    
    @property
    def is_poly(self): #if a compound option can be specified more than once: --foo vals_a --foo vals_b
        if self.options: 
            return self._nargs in self.MULTY_NARGS
        else: return False
    #end def
    
    
    @property
    def decimals(self):
        if not self.limits: return None
        if float(self.limits[0]) > 0:
            return max(ceil(-1*log(self.limits[0], 10)), 2)
        else: return 2
    #end def
    
    
    @property
    def step(self):
        decimals = self.decimals
        if decimals is None: return 1
        return 1.0/(10**ceil(decimals/2.0))
    #end def
    
    
    @property
    def limits_str(self):
        if self.limits is None: return ''
        return '[%s%s%s]' % (('min: %s' % str(self.limits[0]) 
                              if self.limits[0] is not None else ''),
                             (', ' if self.limits[0] is not None 
                              and self.limits[1] is not None else ''),
                             ('max: %s' % str(self.limits[1])
                              if self.limits[1] is not None else ''))
    #end def
    
    
    @property
    def desc(self):
        desc = ''
        if not self.options:
            desc = self._desc.rstrip('. ')
            post = '' 
            if self.default not in (None, '', []):
                post += ' (default: %s%s)' % (str(self.default),
                                             ' '+self.units if self.units else '')
            elif self.units: 
                post += ' (%s)' % self.units 
            if self.limits:
                post += ' %s' % self.limits_str
            if post: desc += ' ' + post
            desc += '.'
        else: 
            desc  = self._desc.rstrip('. ')+'. '
            desc += ('This is a compound option. '
                     'For it exactly %d values should be provided:') % self.nargs
            for i, opt in enumerate(self.options):
                desc += ' %d) %s' % (i+1, opt.desc)
        return desc.lstrip('. ')
    #end def
        
    
    @property
    def action(self):
        class OptionAction(Action):
            fields_list = []
            poly        = False
        
            def __call__(self, parser, namespace, values, option_string=None):
                if len(self.fields_list) != len(values):
                    raise ValueError(('OptionAction: %s: the number of values is not '
                                      'equal to the number of fields.') % self.__class__.__name__)
                _dict = dict()
                for i, val in enumerate(values):
                    _dict[self.fields_list[i][0]] = self.fields_list[i][1](val)
                if not self.poly:
                    setattr(namespace, self.dest, _dict)
                    return
                _dest = getattr(namespace, self.dest)
                if _dest is None:
                    setattr(namespace, self.dest, [_dict])
                else: _dest.append(_dict)
            #end def
        #end class
        
        OptionAction.__name__    = self._name+'_action'
        OptionAction.poly        = self.is_poly
        OptionAction.fields_list = self.fields_list
        
        return OptionAction
    #end def
#end class