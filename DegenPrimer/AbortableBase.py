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
Created on Feb 17, 2014

@author: Allis Tauri <allista@gmail.com>
'''

from abc import ABCMeta, abstractmethod

class AbortableBase(object):
    '''Base class for all that need to be cleanly aborted through an event'''
    __metaclass__ = ABCMeta

    def __init__(self, abort_event):
        self._abort_event = abort_event
    #end def
#end class
        