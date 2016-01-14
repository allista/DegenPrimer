'''
Created on Mar 17, 2014

@author: allis
'''


import os
import sys
import errno
import traceback as tb
import multiprocessing.connection as mpc
from threading import Thread, Event
from time import sleep


class AbortEvent(object): pass

class SignalListener(Thread):
    def __init__(self, connection, event):
        Thread.__init__(self)
        self.daemon = 1
        self._con   = connection
        self._event = event
        self._abort = Event()
    #end def
    
    def run(self):
        try: 
            sig = self._con.recv()
            print sig
            if isinstance(sig, AbortEvent):
                self._event.set()
        except (KeyboardInterrupt, EOFError, IOError): return
        except Exception: tb.print_exc()
    #end def
    
    def stop(self): self._abort.set()
#end class


if __name__ == '__main__':
    a, b = mpc.Pipe()
    e = Event()
    sl = SignalListener(a, e)
    sl.start()
    sleep(1)
    b.send(AbortEvent())
    sleep(1)
    b.close()
#    sl.stop()
    sl.join(1)
    print 'Event is set: %s' % e.is_set()
    print 'Thread is alive: %s' % sl.is_alive()
