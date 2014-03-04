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
from multiprocessing import Queue
from StringTools import print_exception

class OutQueue(object):
    '''A file-like object which puts text written into it 
    in a cross-process queue'''
    def __init__(self, queue=None):
        if queue == None:
            self._queue = Queue()
        else: self._queue = queue
        #saved std-streams
        self._oldout    = None
        self._olderr    = None
    #end def
    
    @property
    def queue(self):
        return self._queue
        
    def write(self, text):
        self._queue.put(text)
    
    def flush(self): pass
    
    
    def __enter__(self):
        self._oldout = sys.stdout
        self._olderr = sys.stderr
        sys.stdout = self
        sys.stderr = self
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