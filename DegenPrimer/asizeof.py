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
Created on Feb 7, 2014

@author: Allis Tauri <allista@gmail.com>
'''

from pympler.asizeof import asizeof
from guppy import hpy

heappy = hpy()

def mem_used(obj, name=None, use_heappy=False):
    Mb_used = asizeof(obj)/1024.0/1024.0
    _type   = type(obj).__name__
    if not name: name = 'unknown'
    print '%s "%s" uses %.02fMb\n' % (_type, name, Mb_used)
    if use_heappy: print heappy.iso(obj)+'\n'
#end def

def sum_mem_used(objects, name, use_heappy=False):
    Mb_used = sum(asizeof(obj) for obj in objects)/1024.0/1024.0
    print '"%s" uses %.02fMb\n' % (name, Mb_used)
    if use_heappy: print heappy.iso(objects)+'\n'
#end def