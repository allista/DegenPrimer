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

import  errno
from AbortableBase import AbortableBase
from PipelineTaskBase import PipelineTaskBase

class Pipeline(AbortableBase):
    '''Run registered tasks in a sequence.'''

    def __init__(self, abort_event):
        AbortableBase.__init__(self, abort_event)
        self._tasks = []
    #end def
        
        
    def register_task(self, task, priority=10):
        if not isinstance(task, PipelineTaskBase):
            raise TypeError('Pipeline.register_task: task should be an '
                            'instance of a subclass of PipelineTaskBase.')
        self._tasks.append((task, priority))
        self._tasks.sort(key=lambda x: x[1])
    #end def
    
    
    def run(self, args):
        result = 0
        for task, _priority in self._tasks:
            if self._abort_event.is_set(): break
            if not task.check_options(args): continue
            args.save_configuration()
            print ''
            try: 
                result = task.run(args)
            except KeyboardInterrupt: break
            except IOError, e: 
                if e.errno == errno.EINTR: break
                else: raise
            except: raise
            if result == 2: continue #task successful, try next task
            if result == 1: break    #task successful, stop here
            if result <  0:
                print ('\nPipeline: error in %s. '
                       '\nExit code %d') % (task.__class__.__name__, result)
                break
        return result
#end class        