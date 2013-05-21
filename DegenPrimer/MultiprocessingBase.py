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
from StringTools import print_exception


class MultiprocessingBase(object):
    '''Base class for classes that use job parallelization through multiprocessing'''
    
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


    #number of cpu cores
    _cpu_count   = mp.cpu_count()


    def __init__(self):
        self._all_jobs     = []
        self._abort_lock   = mp.Lock()
        self._abort_events = []
    #end def
    
    
    @classmethod
    def _join_jobs(cls, abort_e, jobs, timeout, parser, *args):
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
                    print_exception(e)
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


    def _parallelize(self, abort_e, timeout, worker, work, *args):
        #process work in chunks
        work_len   = len(work)
        chunk_size = int(work_len/(self._cpu_count*2)+1)
        results    = [None]*work_len
        jobs       = []
        start      = 0
        while start < work_len and not abort_e.is_set():
            end   = min(start+chunk_size, work_len)
            queue = mp.Queue()
            job   = mp.Process(target=worker, args=(queue, abort_e, start, 
                                                    work[start:end], args))
            job.start()
            jobs.append((job,queue))
            self._all_jobs.append(jobs[-1])
            start += chunk_size
        #join all jobs
        def parse_out(out, results):
            for result in out: 
                results[result[0]] = result[1]
        self._join_jobs(abort_e, list(jobs), timeout, parse_out, results)
        #if aborted, return None
        if abort_e.is_set(): return None
        #else, cleanup
        self._clean_jobs(jobs)
        return results
    #end def
    
    
    def _parallelize_functions(self, abort_e, timeout, func_list, *args):
        @MultiprocessingBase._worker
        def worker(abort_e, offset, _func_list, args):
            result = []
            for fi, func in enumerate(_func_list):
                if abort_e.is_set(): break
                if args:
                    result.append((offset+fi, func(*args)))
                else:
                    result.append((offset+fi, func()))
            return result
        #end def
        return self._parallelize(abort_e, timeout, worker, func_list, *args)
    #end def
    
    
    def _parallelize_work(self, abort_e, timeout, func, data, *args):
        @MultiprocessingBase._worker
        def worker(abort_e, offset, data, args):
            result = []
            for i, item in enumerate(data):
                if abort_e.is_set(): break
                if args:
                    result.append((offset+i, func(item, *args)))
                else:
                    result.append((offset+i, func(item)))
            return result
        #process matches in chunks
        return self._parallelize(abort_e, timeout, worker, data, *args)
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