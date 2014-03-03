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


from multiprocessing.managers import BaseManager, SyncManager
from multiprocessing import Queue, Process
from Queue import Empty
from time import sleep
from random import randint
import gc, sys
#from DegenPrimer.asizeof import heappy


#import multiprocessing
## We must import this explicitly, it is not imported by the top-level
## multiprocessing module.
#import multiprocessing.pool
#import time
#
#from random import randint
#
#
#class NoDaemonProcess(multiprocessing.Process):
#    # make 'daemon' attribute always return False
#    def _get_daemon(self):
#        return False
#    def _set_daemon(self, value):
#        pass
#    daemon = property(_get_daemon, _set_daemon)
#
## We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
## because the latter is only a wrapper function, not a proper class.
#class MyPool(multiprocessing.pool.Pool):
#    Process = NoDaemonProcess
#
#def sleepawhile(t):
#    print("Sleeping %i seconds..." % t)
#    time.sleep(t)
#    return t
#
#def work(num_procs):
#    print("Creating %i (daemon) workers and jobs in child." % num_procs)
#    pool = multiprocessing.Pool(num_procs)
#
#    result = pool.map(sleepawhile,
#        [randint(1, 5) for x in range(num_procs)])
#
#    # The following is not really needed, since the (daemon) workers of the
#    # child's pool are killed when the child is terminated, but it's good
#    # practice to cleanup after ourselves anyway.
#    pool.close()
#    pool.join()
#    return result
#
#def test():
#    print("Creating 5 (non-daemon) workers and jobs in main process.")
#    pool = MyPool(5)
#
#    result = pool.map(work, [randint(1, 5) for x in range(5)])
#
#    pool.close()
#    pool.join()
#    print(result)



class Foo(object):
    def func(self, i, q):
        sleep(3)
        print '\nstarting to generate huge_data'
        huge_data = [randint(1, 10000000)for _i in xrange(300000)]
        print '\nhuge_data was generated'
        q.put(huge_data)
        return 'asdfasdfasdf'
    #end def
    
    def __call__(self, i, q): return self.func(i, q)
#end class


def _foo(i, q): Foo()(i, q)
class MyManager(BaseManager): pass
MyManager.register('foo', _foo, method_to_typeid=dict(foo='str'))


def bar(): return 'asdfasdfasdf'

def foo(i, q):
    fooM = MyManager()
    fooM.start()
    res = fooM.foo(i, q)
    fooM.shutdown()
    return res
#end def
class FooManager(BaseManager): pass
FooManager.register('foo', foo, method_to_typeid=dict(foo='str'))
FooManager.register('bar', bar)

if __name__ == '__main__':
    
    
    mgr = FooManager()
    mgr.start()
    
    res = mgr.bar()
    print res._getvalue()
    print str(res)[0]
    print type(res)
    
    sys.exit(0)
    
    qm = SyncManager()
    qm.start()
    
    mgr = FooManager()
    mgr.start()
    results = []
    for _i in xrange(1):
        q = qm.Queue()
        res = mgr.foo(3, q)
        print res
        print type(res)
        results.append(q.get())
        print '='*80
        sleep(3)
    sleep(5)
    print 'Done.'