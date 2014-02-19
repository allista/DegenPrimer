#!/usr/bin/python
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
Created on Feb 8, 2014

@author: Allis Tauri <allista@gmail.com>
'''

from multiprocessing.managers import BaseManager
from time import sleep
from random import randint
import gc
#from DegenPrimer.asizeof import heappy

class Foo(object):
    def func(self, i):
        sleep(3)
        print '\nstarting to generate huge_data'
        huge_data = [randint(1, 10000000) for _i in xrange(3000000)]
        print '\nhuge_data was generated'
        return huge_data
    #end def
    
    def __call__(self, i): return self.func(i)
#end class

class MyManager(BaseManager): pass
MyManager.register('Foo', Foo)
MyManager.register('gc_collect', gc.collect)

if __name__ == '__main__':
    
    mgr = MyManager()
    mgr.start()
    for _i in xrange(1):
        mfoo = mgr.Foo() 
        res  = mfoo.func(1)
        del mfoo
        print '\nmfoo should now be deleted...'
        sleep(3)
        del res
        print '\nresults should now be deleted...'
        sleep(3)
        gc.collect()
        mgr.gc_collect()
        print '\nall garbage should be collected'
        sleep(3)
    sleep(5)
    print 'Done.'