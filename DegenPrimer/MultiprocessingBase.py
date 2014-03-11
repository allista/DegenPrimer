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
Created on Mar 6, 2013

@author: Allis Tauri <allista@gmail.com>
'''

import errno
import multiprocessing as mp
from UMP import UProcess
from Queue import Empty
from collections import Sequence
from threading import Thread
from AbortableBase import AbortableBase


cpu_count = mp.cpu_count()

#decorators
def worker(func):
    def worker(queue, abort_event, *args, **kwargs):
        result = func(abort_event, *args, **kwargs)
        if abort_event.is_set(): queue.cancel_join_thread()
        else: queue.put(result)
        queue.put(None)
        queue.close()
    #end def
    return worker
#end def
    
def worker_method(func):
    def worker_method(self, queue, abort_event, *args, **kwargs):
        result = func(self, abort_event, *args, **kwargs)
        if abort_event.is_set(): queue.cancel_join_thread()
        else: queue.put(result)
        queue.put(None)
        queue.close()
    #end def
    return worker_method
#end def
    
def data_mapper(func):
    def mapper(queue, abort_event, start, end, data, args):
        for i, item in enumerate(data[start:end]):
            if abort_event.is_set(): break
            if args: result = func(item, *args)
            else: result = func(item)
            queue.put((start+i, result))
        if abort_event.is_set(): queue.cancel_join_thread()
        queue.put(None)
        queue.close()
    return mapper
#end def
    
def data_mapper_method(func):
    def mapper_method(self, queue, abort_event, start, end, data, args):
        for i, item in enumerate(data[start:end]):
            if abort_event.is_set(): break
            if args: result = func(self, item, *args)
            else: result = func(self, item)
            queue.put((start+i, result))
        if abort_event.is_set(): queue.cancel_join_thread()
        queue.put(None)
        queue.close()
    return mapper_method
#end def
    
def results_assembler(func):
    def assembler(result, *args):
        func(result[0], result[1], *args)
    return assembler
#end def
    
def results_assembler_methd(func):
    def assembler_method(self, result, *args):
        func(self, result[0], result[1], *args)
    return assembler_method
#end def


def even_chunks(work, num_jobs):
    work_len = len(work)
    if work_len <= num_jobs: 
        #if work size is less then number of requested jobs, run less jobs
        chunks = [1 for _n in xrange(work_len)]
    else:
        #distribute work evenly between jobs
        ch  = work_len/num_jobs
        rem = work_len%num_jobs
        if rem: 
            chunks   = [ch+1 for _r in xrange(rem)]
            chunks  += [ch for _n in xrange(num_jobs-rem)]
        else: chunks = [ch for _n in xrange(num_jobs)]
    return chunks
#end def

@results_assembler
def ordered_results_assembler(index, result, output):
    output[index] = result
    
@results_assembler
def unordered_results_assembler(index, result, output):
    output.append(result)
    

class Work(Sequence, Thread, AbortableBase):
    
    def __init__(self, abort_event, jobs=None, timeout=1, daemonic=1, **kwargs):
        Thread.__init__(self)
        AbortableBase.__init__(self, abort_event)
        self.daemon     = True
        self._jobs      = jobs if jobs else []
        self._daemonic  = daemonic
        self._timeout   = timeout
        self._assembler = None
        self._args      = None
        self._counter   = kwargs.get('counter', None)
        self._launched  = False
    #end defs
    
    def __len__(self): return len(self._jobs)
    
    def __getitem__(self, index): return self._jobs[index]
    
    def add_job(self, job, queue): self._jobs.append((job, queue))


    def set_assembler(self, assembler, *args):
        self._assembler = assembler
        self._args      = args
    #end def
    
    def set_jobs(self, jobs): self._jobs = jobs
        
    def prepare_jobs(self, worker, work, num_jobs, *args):
        if self._counter is not None: self._counter.set_work(len(work))
        #determine number of jobs to launch
        if not num_jobs: num_jobs = cpu_count
        chunks = even_chunks(work, num_jobs)
        #prepare processes
        self._jobs = []; start = 0
        for chunk in chunks:
            end   = start+chunk
            queue = mp.Queue()
            job   = UProcess(target=worker, args=(queue, self._abort_event, 
                                                  start, end, work, args))
            job.daemon = self._daemonic
            self._jobs.append((job,queue))
            start = end
    #end def
    
    def launch(self):
        if self._launched: return
        for j, _q in self._jobs: j.start()
        self._launched = True
    #end def
    
    def get_result(self):
        assert hasattr(self._assembler, '__call__'), 'Assembler should be a callable'
        while self._jobs:
            finished_job = None
            for i,job in enumerate(self._jobs):
                try: 
                    out = job[1].get(True, self._timeout)
                    if out is not None:
                        self._assembler(out, *self._args)
                        if self._counter is not None: 
                            self._counter.count()
                    else: 
                        job[0].join()
                        finished_job = i
                    break
                except Empty:
                    if not job[0].is_alive():
                        job[0].join()
                        finished_job = i
                    continue
                except IOError, e:
                    if e.errno == errno.EINTR:
                        continue
                    else:
                        print 'Unhandled IOError:', e.message
                        raise
                except Exception, e:
                    print 'Unhandled Exception:'
                    print str(e)
                    raise
            if finished_job is not None:
                del self._jobs[finished_job]
    #end def
    get_nowait = Thread.start
    
    def run(self):
        try: 
            self.launch()
            self.get_result()
        #Ctrl-C
        #EOF means that target activity was terminated in the Manager process
        except (KeyboardInterrupt, EOFError): return
        #IO code=4 means the same
        except IOError, e:
            if e.errno == errno.EINTR: return
            raise
    #end def
    
    def wait(self):
        if self.is_alive(): self.join()
        return not self._abort_event.is_set()
#end class        
        

class MultiprocessingBase(AbortableBase):
    '''Base class for classes that use job parallelization through multiprocessing'''

    #aliases
    _Process = staticmethod(UProcess)
    _Queue   = staticmethod(mp.Queue)
    
    worker                  = staticmethod(worker)
    worker_method           = staticmethod(worker_method)
    data_mapper             = staticmethod(data_mapper)
    data_mapper_method      = staticmethod(data_mapper_method)
    results_assembler       = staticmethod(results_assembler)
    results_assembler_methd = staticmethod(results_assembler_methd)
    
    ordered_results_assembler   = staticmethod(ordered_results_assembler)
    unordered_results_assembler = staticmethod(unordered_results_assembler)
    
    #class body
    def __init__(self, abort_event, daemonic=True):
        AbortableBase.__init__(self, abort_event)
        self._daemonic = daemonic
    #end def


    def Work(self, jobs=None, timeout=1, **kwargs): 
        return Work(self._abort_event, jobs, timeout, self._daemonic, **kwargs)

    def start_work(self, *work):
        for w in work: w.launch()
        
    def wait(self, *work):
        for w in work: w.get_nowait()
        ret = True
        for w in work: ret &= w.wait()
        return ret
    #end def
        
        
    def parallelize(self, timeout, worker, work, *args, **kwargs):
        #prepare and start jobs
        w = self.Work(**kwargs)
        w.prepare_jobs(worker, work, None, *args)
        w.launch()
        #allocate result container, get result
        result = [None]*len(work)
        w.set_assembler(ordered_results_assembler, result)
        w.get_nowait()
        #if aborted, return None
        if not w.wait(): return None
        return result
    #end def

    def parallelize_work(self, timeout, func, data, *args, **kwargs):
        mapper = MultiprocessingBase.data_mapper(func)
        return self.parallelize(timeout, mapper, data, *args, **kwargs)
    #end def    
    
    def parallelize_functions(self, timeout, func_list, *args, **kwargs):
        @MultiprocessingBase.worker
        def worker(abort_event, start, end, _func_list, args):
            result = []
            for fi, func in enumerate(_func_list[start:end]):
                if abort_event.is_set(): break
                if args:
                    result.append((start+fi, func(*args)))
                else:
                    result.append((start+fi, func()))
            return result
        #end def
        return self.parallelize(timeout, worker, func_list, *args, **kwargs)
    #end def
#end class


#decorators and shortcuts
from tmpStorage import tmpDict
from UMP import UManager

def shelf_result(func, *args, **kwargs):
    '''Decorator that saves result of func(*args) into a temporary shelf 
    under "result" key and returns the path to it's database'''
    def s_func(*args, **kwargs):
        result = func(*args, **kwargs)
        if result is None: return None
        d = tmpDict(persistent=True)
        d['result'] = result 
        d.close()
        return d.filename
    return s_func
#end def 

def queue_result(queue, func, *args, **kwargs):
    '''Decorator that puts result of func(*args) into a queue'''
    def q_func(*args, **kwargs): queue.put(func(*args, **kwargs))
    return q_func
#end def


#shortcut functions
def parallelize_work(abort_event, daemonic, timeout, func, data, *args, **kwargs):
    '''Parallel map implementation by MultiprocessingBase'''
    return MultiprocessingBase(abort_event, daemonic).parallelize_work(timeout, func, data, *args, **kwargs)

def parallelize_functions(abort_event, daemonic, timeout, func_list, *args, **kwargs):
    '''Execute functions from func_list in parallel using MultiprocessingBase.'''
    return MultiprocessingBase(abort_event, daemonic).parallelize_functions(timeout, func_list, *args, **kwargs)

def FuncManager(funcs, names=None):
    '''Return a subclass of UManager with funcs registered by names (if given)
    or by (func.__name__).lstrip('_') for func in funcs''' 
    class _FuncManager(UManager): pass
    for i, func in enumerate(funcs):
        try: name = names[i]
        except: name = (func.__name__).lstrip('_')
        _FuncManager.register(name, func)
    return _FuncManager
#end def

Parallelizer = FuncManager((shelf_result(parallelize_work),
                            shelf_result(parallelize_functions)), 
                           ('parallelize_work',
                            'parallelize_functions'))



#tests
if __name__ == '__main__':
    from multiprocessing import Event
    abort_event = Event()
    
    class Test(MultiprocessingBase):
        def __init__(self, abort_event):
            MultiprocessingBase.__init__(self, abort_event)
            self._a = 5
            
        def _method(self, i):
            return self._a+abs(i)
        
        def method(self, data):
            return self.parallelize_work(1, self._method, data)
    
    t = Test(abort_event)
    print t.method([-1, -2, -3, -4, -5, -6, -7, -8])