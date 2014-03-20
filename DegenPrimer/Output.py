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
Created on Mar 4, 2014

@author: Allis Tauri <allista@gmail.com>
'''

import sys
import traceback
from multiprocessing.queues import Queue
from StringTools import print_exception


class OutIntercepter(object):
    '''A file-like object which intercepts std-out/err'''
    def __init__(self):
        self._oldout = None
        self._olderr = None
    #end def
    
    def write(self, text): pass
    
    def flush(self): pass
    
    def __enter__(self):
        self._oldout = sys.stdout
        self._olderr = sys.stderr
        sys.stdout = sys.stdout = self
        return self
    #end def
    
    def __exit__(self, _type, _value, _traceback):
        if _type is not None:
            print_exception(_value)
            traceback.print_exception(_type, _value, _traceback, file=self._olderr)
        sys.stdout = self._oldout
        sys.stderr = self._olderr
        return True
    #end def
#end class


class OutQueue(OutIntercepter, Queue):
    '''A file-like object which puts text written into it 
    in a cross-process queue'''
    def __init__(self, maxsize=0):
        OutIntercepter.__init__(self)
        Queue.__init__(self, maxsize)
        self._id = None
    #end def

    def set_id(self, _id):
        self._id = _id
        return self
    #end def

    def put(self, obj, block=True, timeout=None):
        Queue.put(self, (self._id, obj), block=block, timeout=timeout)
        
    def write(self, text): self.put(text)
    
    def __enter__(self):
        OutIntercepter.__enter__(self)
        if self._id is None: self._id = -1
        return self
    #end def
    
    def __exit__(self, _type, _value, _traceback):
        OutIntercepter.__exit__(self, _type, _value, _traceback)
        self._id = None
        return True
    #end def
#end class