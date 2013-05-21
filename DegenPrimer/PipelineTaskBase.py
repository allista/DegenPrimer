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


class PipelineTaskBase(object):
    '''Base class for pipeline tasks'''
    
    def __del__(self):
        self.terminate()
    
    @staticmethod
    def check_options(args): pass
    
    def terminate(self): pass
    
    def run(self, args): pass
#end class