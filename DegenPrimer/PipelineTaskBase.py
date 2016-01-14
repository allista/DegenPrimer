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
Created on Mar 28, 2013

@author: Allis Tauri <allista@gmail.com>
'''

from BioUtils.Tools.AbortableBase import AbortableBase
from abc import ABCMeta, abstractmethod


class PipelineTaskBase(AbortableBase):
    '''Base class for pipeline tasks'''
    __metaclass__ = ABCMeta
    
    @staticmethod
    @abstractmethod
    def check_options(args): pass
    
    def __init__(self, abort_event):
        AbortableBase.__init__(self, abort_event)
        self._newline_last = False
    #end def
    
    def _print(self, text=''):
        if self._newline_last:
            if not text: return
            text = text.lstrip('\n')
        if text == '' or text.endswith('\n'):
            self._newline_last = True
        else: self._newline_last = False
        print text
    #end def
#end class