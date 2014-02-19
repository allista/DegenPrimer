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
import signal
import multiprocessing as mp
from multiprocessing.managers import BaseManager 
from Queue import Empty
from AbortableBase import AbortableBase


def ignore_interrupt():
    signal.signal(signal.SIGINT,  signal.SIG_IGN)
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
    signal.signal(signal.SIGQUIT, signal.SIG_IGN)
#end def


class UManager(BaseManager):
    '''Multiprocessing Manager that ignores interrupt, term and quit signals'''
    def start(self, initializer=ignore_interrupt, initargs=()):
        BaseManager.start(self, initializer=initializer, initargs=initargs)
#end class


class UProcess(mp.Process):
    '''Process that ignores interrupt, term and quit signals'''
    def run(self):
        ignore_interrupt()
        mp.Process.run(self)
#end class


class MultiprocessingBase(AbortableBase):
    '''Base class for classes that use job parallelization through multiprocessing'''

    #aliases
    _Process   = staticmethod(UProcess)
    _Queue     = staticmethod(mp.Queue)
    _Lock      = staticmethod(mp.Lock)
    
    #number of cpu cores
    cpu_count  = mp.cpu_count()
    
    #decorators
    @staticmethod
    def _worker(func):
        def worker(queue, abort_event, *args):
            result = None
            if args: result = func(abort_event, *args)
            else: result = func(abort_event)
            if abort_event.is_set(): queue.cancel_join_thread()
            else: queue.put(result)
            queue.put(None)
            queue.close()
        #end def
        return worker
    #end def
    
    @staticmethod
    def _worker_method(func):
        def worker_method(self, queue, abort_event, *args):
            result = None
            if args: result = func(self, abort_event, *args)
            else: result = func(self, abort_event)
            if abort_event.is_set(): queue.cancel_join_thread()
            else: queue.put(result)
            queue.put(None)
            queue.close()
        #end def
        return worker_method
    #end def
    
    @staticmethod
    def _data_mapper(func):
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
    def _data_mapper_method(func):
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
    def _results_assembler(func):
        def assembler(result, *args):
            func(result[0], result[1], *args)
        return assembler
    #end def
    
    @staticmethod
    def _results_assembler_methd(func):
        def assembler_method(self, result, *args):
            func(self, result[0], result[1], *args)
        return assembler_method
    #end def
    

    #class body
    def __init__(self, abort_event):
        AbortableBase.__init__(self, abort_event)
        self._all_jobs     = []
        self._abort_lock   = None
    #end def
    
    def __del__(self): self._abort()
    
    
    def _join_jobs(self, jobs, timeout, parser, *args):
        jobs = list(jobs)
        while jobs:
            finished_job = None
            for i,job in enumerate(jobs):
                try: 
                    out = job[1].get(True, timeout)
                    if out is not None:
                        parser(out, *args)
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
    
    
    def _register_job(self, job):
        if self._abort_lock is None:
            self._abort_lock = self._Lock()
        with self._abort_lock:
            self._all_jobs.append(job)
    #end def
    
    
    def _clean_jobs(self, jobs):
        if not self._all_jobs:
            raise ValueError('MultiprocessingBase._clean_jobs: no jobs has been regestered.')
        with self._abort_lock:
            for job in jobs: 
                if job in self._all_jobs:
                    self._all_jobs.remove(job)
        if not self._all_jobs: self._abort_lock = None
    #end def


    def _prepare_jobs(self, worker, work, num_jobs, *args):
        #process work in chunks
        work_len = len(work)
        jobs     = []
        start    = 0
        if not num_jobs: num_jobs = self.cpu_count*2
        chunk_size = max(work_len/num_jobs, 1)
        while start < work_len:
            end   = min(start+chunk_size, work_len)
            queue = self._Queue()
            job   = self._Process(target=worker, args=(queue, self._abort_event, 
                                                       start, end, 
                                                       work, args))
            job.daemon = 1
            jobs.append((job,queue))
            self._register_job(jobs[-1])
            start += chunk_size
        return jobs
    #end def
    
    
    @staticmethod
    @_results_assembler.__func__
    def _ordered_results_assembler(index, result, output):
        output[index] = result
        
    @staticmethod
    @_results_assembler.__func__
    def _unordered_results_assembler(index, result, output):
        output.append(result)
    
    def _start_jobs(self, jobs):
        for job, _queue in jobs: job.start()
    

    def _parallelize(self, timeout, worker, work, *args):
        #prepare and start jobs
        jobs = self._prepare_jobs(worker, work, None, *args)
        self._start_jobs(jobs)
        #allocate results container
        results = [None]*len(work)
        #assemble results
        self._join_jobs(jobs, timeout, self._ordered_results_assembler, results)
        #if aborted, return None
        if self._abort_event.is_set(): return None
        #else, cleanup and return results
        self._clean_jobs(jobs)
        return results
    #end def

    def _parallelize_work(self, timeout, func, data, *args):
        mapper = MultiprocessingBase._data_mapper(func)
        return self._parallelize(timeout, mapper, data, *args)
    #end def    
    
    def _parallelize_functions(self, timeout, func_list, *args):
        @MultiprocessingBase._worker
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
        return self._parallelize(timeout, worker, func_list, *args)
    #end def
    
    
    def _abort(self):
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