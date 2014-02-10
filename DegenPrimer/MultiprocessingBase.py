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
from Queue import Empty
from time import sleep


class MultiprocessingBase(object):
    '''Base class for classes that use job parallelization through multiprocessing'''

    #aliases
    _Process = mp.Process
    
    #number of cpu cores
    _cpu_count   = mp.cpu_count()
    
    
    #decorators
    @staticmethod
    def _worker(func):
        def worker(queue, abort_e, *args):
            result = None
            if args: result = func(abort_e, *args)
            else: result = func(abort_e)
            if abort_e.is_set():
                queue.cancel_join_thread()
            else: queue.put(result)
            queue.close()
            return
        #end def
        return worker
    #end def
    
    @staticmethod
    def _data_mapper(func):
        @MultiprocessingBase._worker
        def mapper(abort_e, start, end, data, args):
            result = []
            for i, item in enumerate(data[start:end]):
                if abort_e.is_set(): break
                if args:
                    result.append((start+i, func(item, *args)))
                else:
                    result.append((start+i, func(item)))
            return result
        return mapper
    #end def


    #class body
    def __init__(self):
        self._all_jobs     = []
        self._abort_lock   = mp.Lock()
        self._abort_events = []
    #end def
    
    def __del__(self): self.abort()
    
    
    @classmethod
    def _join_jobs(cls, abort_e, jobs, timeout, parser, *args):
        jobs = list(jobs)
        while jobs and not abort_e.is_set():
            finished_job = None
            for i,job in enumerate(jobs):
                try: 
                    out = job[1].get(True, timeout)
                    parser(out, *args)
                    job[0].join()
                    finished_job = i
                    break
                except Empty: 
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
    
    
    @classmethod
    def _new_queue(cls): return mp.Queue()
    
    
    def _new_event(self):
        abort_event = mp.Event()
        self._abort_events.append(abort_event)
        return abort_event
    #end def
    
    
    def _clean_jobs(self, jobs):
        with self._abort_lock:
            for job in jobs: 
                if job in self._all_jobs:
                    self._all_jobs.remove(job)
    #end def
    
    
    def _clean_abort_even(self, abort_event):
        with self._abort_lock:
            if abort_event in self._abort_events:
                self._abort_events.remove(abort_event)
    #end def


    def _prepare_jobs(self, abort_e, worker, work, num_jobs, *args):
        #process work in chunks
        work_len = len(work)
        jobs     = []
        start    = 0
        if not num_jobs: num_jobs = self._cpu_count*2
        chunk_size = int(work_len/(num_jobs)+1)
        #prepare jobs
        while start < work_len and not abort_e.is_set():
            end   = min(start+chunk_size, work_len)
            queue = self._new_queue()
            job   = self._Process(target=worker, args=(queue, abort_e, 
                                                       start, end, 
                                                       work, args))
            jobs.append((job,queue))
            self._all_jobs.append(jobs[-1])
            start += chunk_size
        #if aborted, return None
        if abort_e.is_set(): return None
        return jobs
    #end def
    
    
    @classmethod
    def _ordered_results_assembler(cls, out, results):
        for result in out: results[result[0]] = result[1]
        
    @classmethod
    def _unordered_results_assembler(cls, out, results):
        for result in out: results.append(result[1])
    
    @classmethod
    def _start_jobs(cls, jobs):
        for job, _queue in jobs: job.start()
    

    def _parallelize(self, abort_e, timeout, worker, work, *args):
        #prepare and start jobs
        jobs = self._prepare_jobs(abort_e, worker, work, None, *args)
        self._start_jobs(jobs)
        #allocate results container
        results = [None]*len(work)
        #assemble results
        self._join_jobs(abort_e, jobs, timeout, 
                        self._ordered_results_assembler, results)
        #if aborted, return None
        if abort_e.is_set(): return None
        #else, cleanup and return results
        self._clean_jobs(jobs)
        return results
    #end def

    def _parallelize_work(self, abort_e, timeout, func, data, *args):
        mapper = MultiprocessingBase._data_mapper(func)
        return self._parallelize(abort_e, timeout, mapper, data, *args)
    #end def    
    
    def _parallelize_functions(self, abort_e, timeout, func_list, *args):
        @MultiprocessingBase._worker
        def worker(abort_e, start, end, _func_list, args):
            result = []
            for fi, func in enumerate(_func_list[start:end]):
                if abort_e.is_set(): break
                if args:
                    result.append((start+fi, func(*args)))
                else:
                    result.append((start+fi, func()))
            return result
        #end def
        return self._parallelize(abort_e, timeout, worker, func_list, *args)
    #end def
    
    
    def abort(self):
        if self._abort_events:
            with self._abort_lock:
                for event in self._abort_events:
                    event.set()
            sleep(1)
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
                self._abort_events = []
    #end def
#end class        