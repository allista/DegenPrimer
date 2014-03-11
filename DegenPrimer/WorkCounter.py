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

import warnings
import numpy as np
from scipy import stats
from time import time
from datetime import timedelta
from collections import Sequence, Sized
from abc import ABCMeta


class ReturnProxy(type):
    '''A metaclass to work around http://bugs.python.org/issue20854'''
    def __call__(cls, *args, **kwargs):
        if not kwargs and args and isinstance(args[0], cls):
            return args[0]
        return super(ReturnProxy, cls).__call__(*args, **kwargs)
#end metaclass

class SequenceReturnProxy(ReturnProxy, ABCMeta): pass


class WorkCounter(Sequence):
    '''A counter that reports percentage of !recursive! work done'''
    __metaclass__ = SequenceReturnProxy

    def __init__(self, work=1, weights=None, parent=None, subwork=False):
        self._parent      = parent
        self._subcounters = []
        self.set_work(work, weights, subwork)
        self._last_ratio  = None
        self._changed     = False
        self._history     = []
    #end def
    
    
    def __str__(self):
        s  = '\n'
        s += 'Counter ID: %x\n' % id(self)
        s += 'Work size: %f\n' % self._work
        if self._own_work:
            s += '    Own work: %f\n' % self._own_work
        if self._subwork:
            s += '    Subwork: %f\n' % self._subwork
            s += '    %s\n' % self._weights
        progress = self.progress()
        percent  = progress/float(self._work)*100
        s += 'Progress: %f; %6.2f%%\n' % (progress, percent)
        s += 'Elapsed: %s' % timedelta(seconds=self.elapsed())
        s += '; ETA: %s' % self.ETA()
        return s
    #end def
    
    def rstr(self, level=0):
        '''String representation of a counter and all its' subcounters,
        recursively.'''
        s  = str(self)
        if not self._subcounters: return s
        s += '\n%3d %s' % (level+1, '-'*(80-4*(level+1)))
        for i, c in enumerate(self._subcounters):
            s += c.rstr(level+1).replace('\n', '\n%3d ' % (i+1))
            s += '\n'
        hl = ('.' if level else '=')*(80-4*level)
        s = '\n%s%s%s' % (hl,s,hl) 
        return s
    #end def
    
    def __repr__(self): return self.rstr()
    
    
    def __nonzero__(self): return True
    
    def __len__(self): return len(self._subcounters)
    
    def __getitem__(self, index): return self._subcounters[index]
    
    
    def set_parent(self, parent): self._parent = parent
            
    def set_work(self, work, weights=None, subwork=False):
        #full reset
        self.reset(); self._subcounters = []; self._weights = []
        #check work size
        if isinstance(work, Sized): work = len(work)
        assert work > 0, \
        ('Counter: %x: work size should be grater than zero' % id(self))
        #set work size
        if subwork: #create subcounters
            self._work     = 0
            self._own_work = 0
            self._subwork  = 0
            self.add_subwork(work, weights)
        else: 
            self._work     = work
            self._own_work = work
            self._subwork  = 0
    #end def
    
    def set_subwork(self, work, weights=None):
        self.set_work(work, weights=weights, subwork=True)
        
    def add_subwork(self, work=1, weights=None):
        if not weights: weights = [1]*work
        assert len(weights) == work, \
        ('Counter: %x: a weight should be provided for every subwork.' % id(self))
        assert min(weights) > 0, 'Weights should be positive.'
        subwork = sum(weights)
        self._work    += subwork
        self._subwork += subwork
        self._weights.extend(weights)
        for _w in xrange(work):
            self._subcounters.append(WorkCounter(parent=self))
    #end def
    
    
    def elapsed(self): return time()-self._time0
    
    def last_checked(self): return time()-self._time_1
    
    def changed(self): return self._changed
    
    
    def count(self, num=1):
        assert self._count+num <= self._own_work, \
        ('Counter: %x: the work counted exceeds total work amount: %.8f+%.8f>%.8f' %
         (id(self), self._count, num, self._own_work))
        self._count += num
    #end def
    
    def done(self): self._count = self._work-self._subwork
    
    def reset(self): 
        for counter in self._subcounters: counter.reset()
        self._count  = 0
        self._time0  = time()
        self._time_1 = self._time0
    #end def
    
    
    def total(self): return self._work

    def progress(self):
        p = (sum(c.ratio()*w 
                 for w, c in zip(self._weights, self._subcounters)) 
             + self._count)
        dt  = time()-self._time0
        ddt = dt if not self._history else dt - self._history[-1][0]
        dp  = p if not self._history else p - self._history[-1][1]
        self._history.append((dt, p, ddt, dp))
        return p
    #end def
        

    def ratio(self):
        ratio = self.progress()/float(self._work)
        self._changed    = (ratio != self._last_ratio)
        self._last_ratio = ratio
        return ratio
    #end def
    
    def percent(self):
        self._time_1 = time() 
        return self.ratio()*100.0
    
    def changed_percent(self):
        percent = self.percent()
        if self._changed: return percent
        return None
    #end def
    
#    def ETA(self):
#        if not self._history or not self._history[-1][1]: return '-:--:--'
#        progress = self._history[-1][1]
#        if progress == self._work: return '0:00:00'
#        dt = round((self._work/progress-1)*self._history[-1][0])
#        return str(timedelta(seconds=dt))
#    #end def
    
    
#    def ETA(self):
#        if not self._history or not self._history[-1][1]: return '-:--:--'
#        progress = self._history[-1][1]
#        if progress == self._work: return '0:00:00'
#        v  = np.array([dp/float(ddt) for _dt, _p, ddt, dp in self._history if dp > 0])
#        if len(v) > 1:
#            with warnings.catch_warnings():
#                warnings.simplefilter("ignore")
#                v_mean, _var, _std = stats.bayes_mvs(v, 0.9)
#            rtl = round((self._work-progress)/v_mean[1][1])
#            rtu = round((self._work-progress)/v_mean[1][0])
#        else: 
#            rtl = round((self._work-progress)/v[0])
#            rtu = -1
#        if rtu > 0 and rtu-rtl > 1:
#            return 'from %s to %s' % (str(timedelta(seconds=rtl)),
#                                             str(timedelta(seconds=rtu)))
#        else: return 'approximately %s' % str(timedelta(seconds=rtl))
    #end def
    
    def ETA(self):
        if not self._history or not self._history[-1][1]: return '-:--:--'
        progress = self._history[-1][1]
        if progress == self._work: return '0:00:00'
        h  = np.array(self._history, dtype=float).transpose()
        v  = h[3]/h[2]
        m  = np.mean(v)
        if len(v) > 1:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                lm = linregress(h[0], h[1])
            rtl = round((self._work-progress)/lm[0])
            rtu = round((self._work-progress)/m)
        else: 
            rtl = round((self._work-progress)/m)
            rtu = -1
        if rtu > 0 and rtu-rtl > 1:
            return 'from %s to %s' % (str(timedelta(seconds=rtl)),
                                      str(timedelta(seconds=rtu)))
        else: return 'approximately %s' % str(timedelta(seconds=rtl))
    #end def
#end class


#proxy for WorkerCounter
from UMP import UManager
from multiprocessing.managers import BaseProxy

class AutoProxyMeta(type):
    '''Metaclass that replicates multiprocessing.managers.MakeProxyType
    functionality, but allows proxy classes that use it to be pickable'''
    def __new__(cls, name, bases, attrs):
        dic = {}
        for meth in attrs.get('_exposed_', ()):
            exec '''def %s(self, *args, **kwds):
            return self._callmethod(%r, args, kwds)''' % (meth, meth) in dic
        dic.update(attrs)
        return super(AutoProxyMeta, cls).__new__(cls, name, bases, dic)
#end metaclass
        
class WorkCounterProxy(BaseProxy):
    __metaclass__ = AutoProxyMeta
    _exposed_ = ('__str__', 'rstr', '__repr__',
                 '__nonzero__', '__len__', '__getitem__',
                 'set_work', 'set_subwork', 'add_subwork',
                 'elapsed',  'last_checked', 'ETA', 'count', 'done',
                 'total', 'percent', 'changed_percent',)
    _method_to_typeid_ = dict(__getitem__='WorkCounter')
    
    def __nonzero__(self): return True
    
    def __reduce__(self):
        _unpickle, (cls, token, serializer, kwds) = BaseProxy.__reduce__(self)
        kwds['manager'] = self._manager
        return _unpickle, (cls, token, serializer, kwds)
    #end def
#end class

class WorkCounterManager(UManager):
    def __reduce__(self):
        return (RebuildWorkCounterManager,
                (self._address, None, self._serializer))
WorkCounterManager.register('WorkCounter', WorkCounter, WorkCounterProxy)

def RebuildWorkCounterManager(address, authkey, serializer):
    mgr = WorkCounterManager(address, authkey, serializer)
    mgr.connect()
    return mgr
#end def



#tests
from matplotlib import pyplot as plt
from scipy.stats import linregress, mode, hmean, gmean
def plot_history(counter):
    h = np.array(counter._history, dtype=float).transpose()
    v = h[3]/h[2]
    m = np.mean(v)#; gm = gmean(v)
    lm = linregress(h[0], h[1])
    plt.plot(h[0], h[1], 'b-', 
             h[0], h[0]*lm[0]+lm[1], 'r-',
             h[0], h[0]*lm[0], 'r--',
             h[0], h[0]*m, 'g-',
#             h[0], h[0]*hm, 'g--',
#             h[0], h[0]*gm, 'g.-',
             ) 
    plt.show()
#end def

if __name__ == '__main__':
    import sys
    from time import sleep
    import cPickle as pickle
    
    cmgr = WorkCounterManager(); cmgr.start()
    c = WorkCounter(10)
    c.add_subwork(10, weights=(1.1,2.2,3,4,5,6,7,8,9,10))
    
    for i in xrange(10):
        print c
        c[i].count()
        c.count()
        sleep(1)
    print c
    
    plot_history(c)
    