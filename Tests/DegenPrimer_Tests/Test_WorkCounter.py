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
Created on 2016-01-14

@author: Allis Tauri <allista@gmail.com>
'''

from time import sleep
from DegenPrimer.WorkCounter import WorkCounter, WorkCounterManager

def test():
    cmgr = WorkCounterManager(); cmgr.start()
    c = WorkCounter(10)
    c.add_subwork(10, weights=(1.1,2.2,3,4,5,6,7,8,9,10))
    
    for i in xrange(10):
        print c
        c[i].count()
        c.count()
        sleep(1)
    print c