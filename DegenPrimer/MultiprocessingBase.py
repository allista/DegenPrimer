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


class Work(Sequence):
    def __init__(self, assembler, args=None, counter=None):
        self._jobs      = None
        self._counter   = counter
        self._assembler = assembler
        self._args      = args
    #end defs
        
    def __len__(self): return len(self._jobs)
    
    def __getitem__(self, index): return self._jobs[index]
        
    def prepare_jobs(self, worker, work, num_jobs, *args):
        #determine number of jobs to launch
        if not num_jobs: num_jobs = self.cpu_count
        chunks = self.even_chunks(work, num_jobs)
        #prepare processes
        self._jobs  = []; start = 0
        for chunk in chunks:
            end   = start+chunk
            queue = self._Queue()
            job   = self._Process(target=worker, args=(queue, self._abort_event, 
                                                       start, end, 
                                                       work, args))
            job.daemon = self._daemonic
            self._jobs.append((job,queue))
            start = end
    #end def
    
    def start(self): 
        for j, _q in self._jobs: j.start()
#end class        
        

class MultiprocessingBase(AbortableBase):
    '''Base class for classes that use job parallelization through multiprocessing'''

    #aliases
    _Process   = staticmethod(UProcess)
    _Queue     = staticmethod(mp.Queue)
    _Lock      = staticmethod(mp.Lock)
    
    #number of cpu cores
    cpu_count  = mp.cpu_count()
    
    
    #class body
    def __init__(self, abort_event, daemonic=True):
        AbortableBase.__init__(self, abort_event)
        self._all_jobs   = []
        self._abort_lock = None
        self._daemonic   = daemonic
    #end def
    
    def __del__(self): self.abort()
    
    
    #decorators
    @staticmethod
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
    
    @staticmethod
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
    
    @staticmethod
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
    
    @staticmethod
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
    
    @staticmethod
    def results_assembler(func):
        def assembler(result, *args):
            func(result[0], result[1], *args)
        return assembler
    #end def
    
    @staticmethod
    def results_assembler_methd(func):
        def assembler_method(self, result, *args):
            func(self, result[0], result[1], *args)
        return assembler_method
    #end def
    
    
    def get_result(self, jobs, timeout, parser, *args, **kwargs):
        try: counter = kwargs['counter']
        except KeyError: counter = None
        jobs = list(jobs)
        while jobs:
            finished_job = None
            for i,job in enumerate(jobs):
                try: 
                    out = job[1].get(True, timeout)
                    if out is not None:
                        parser(out, *args)
                        if counter: counter.count()
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
                del jobs[finished_job]
    #end def
    
    
    def register_job(self, job):
        if self._abort_lock is None:
            self._abort_lock = self._Lock()
        with self._abort_lock:
            self._all_jobs.append(job)
    #end def
    
    
    def clean_jobs(self, jobs):
        if not self._all_jobs:
            raise ValueError('MultiprocessingBase.clean_jobs: no jobs has been regestered.')
        with self._abort_lock:
            for job in jobs: 
                if job in self._all_jobs:
                    self._all_jobs.remove(job)
        if not self._all_jobs: self._abort_lock = None
    #end def


    @staticmethod
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

    def prepare_jobs(self, worker, work, num_jobs, *args):
        #determine number of jobs to launch
        if not num_jobs: num_jobs = self.cpu_count
        chunks = self.even_chunks(work, num_jobs)
        #prepare processes
        jobs  = []; start = 0
        for chunk in chunks:
            end   = start+chunk
            queue = self._Queue()
            job   = self._Process(target=worker, args=(queue, self._abort_event, 
                                                       start, end, 
                                                       work, args))
            job.daemon = self._daemonic
            jobs.append((job,queue))
            self.register_job(jobs[-1])
            start = end
        return jobs
    #end def
    
    
    @staticmethod
    @results_assembler.__func__
    def ordered_results_assembler(index, result, output):
        output[index] = result
        
    @staticmethod
    @results_assembler.__func__
    def unordered_results_assembler(index, result, output):
        output.append(result)
    
    def start_jobs(self, *jobs):
        for job in jobs:
            for j, _q in job: j.start()
    #end def
    
    def parallelize(self, timeout, worker, work, *args, **kwargs):
        #setup counter, if it is given
        try: kwargs['counter'].set_work(len(work))
        except: pass
        #prepare and start jobs
        jobs = self.prepare_jobs(worker, work, None, *args)
        self.start_jobs(jobs)
        #allocate results container
        result = [None]*len(work)
        #assemble results
        self.get_result(jobs, timeout, self.ordered_results_assembler, result, **kwargs)
        #if aborted, return None
        if self._abort_event.is_set(): return None
        #else, cleanup and return result
        self.clean_jobs(jobs)
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
    
    
    def abort(self):
        if self._all_jobs:
            with self._abort_lock:
                while self._all_jobs:
                    finished_job = None
                    for i,job in enumerate(self._all_jobs):
                        try: 
                            while True: job[1].get_nowait()
                        except Empty: pass
                        if not job[0].is_alive():
                            job[0].join()
                            finished_job = i
                            break
                    if finished_job is not None:
                        del self._all_jobs[finished_job]
            self._abort_lock = None
    #end def
#end class


#aliases
even_chunks = MultiprocessingBase.even_chunks
cpu_count   = MultiprocessingBase.cpu_count

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