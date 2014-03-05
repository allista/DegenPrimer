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
Created on Mar 4, 2014

@author: Allis Tauri <allista@gmail.com>
'''

from time import time
from collections import Sequence, Sized
from multiprocessing.managers import BaseProxy

class WorkCounter(Sequence):
    '''A counter that reports percentage of !recursive! work done'''

    def __init__(self, work_len=1, parent=None, subcounters=None, auto_subcounters=False):
        self._parent      = parent
        self._subcounters = subcounters if subcounters else []
        self.set_len(work_len, auto_subcounters)
        self._last_ratio  = None
        self._changed     = False
    #end def
    
    
    def __nonzero__(self): return True
    
    def __len__(self): return len(self._subcounters)
    
    def __getitem__(self, index): return self._subcounters[index]
            
    def set_len(self, work_len, auto_subcounters=False):
        self.reset()
        if isinstance(work_len, Sized): 
            self._len = len(work_len)
        else: self._len = work_len
        if auto_subcounters:
            if self._subcounters:
                raise RuntimeError('WorkCounter: you cannot use auto_subcounters '
                                   'if a counter already has some subcounters.')
            for _i in xrange(self._len):
                self._subcounters.append(WorkCounter(parent=self))
    #end def
    
    
    def elapsed(self): return time()-self._time0
    
    def changed(self): return self._changed
    
    
    def count(self, num=1):
        if self._count+num > self._len-len(self._subcounters):
            raise ValueError('WorkCounter: the work counted exceeds total work amount') 
        self._count += num
    #end def
    
    def done(self): self._count = self._len-len(self._subcounters)
    
    def reset(self): 
        for counter in self._subcounters: counter.reset()
        self._count = 0
        self._time0 = time()
    #end def
    
    
    def set_parent(self, parent): 
        self._parent = parent
    
    def add_subcounter(self, counter):
        if len(self._subcounters) == self._len:
            raise ValueError('WorkCounter: all subcounters are already present.') 
        self._subcounters.append(counter)
        counter.set_parent(self)
    #end def
    
    
    def total(self): return self._len

    def progress(self):
        self._time0 = time()
        return sum(c.ratio() for c in self._subcounters) + self._count

    def ratio(self):
        self._time0 = time()
        ratio = (sum(c.ratio() for c in self._subcounters)+self._count)/float(self._len)
        self._changed    = ratio == self._last_ratio
        self._last_ratio = ratio
        return ratio
    #end def
    
    def percent(self): return self.ratio()*100.0
    
    def changed_percent(self):
        percent = self.percent()
        return percent if self._changed else None
    #end def
#end class



#tests
if __name__ == '__main__':
    from time import sleep
    c = WorkCounter(10, auto_subcounters=True)
    print c
    if c: print 'c!'
    else: print 'no c!'
    
    for i in xrange(10):
        print 'Elapsed: %.03fs' % c.elapsed()
        print 'Work done: %.02f%%' % c.percent()
        c[i].count()
        sleep(0.5)
    print 'Elapsed: %.03fs' % c.elapsed()
    print 'Work done: %.02f%%' % c.percent()