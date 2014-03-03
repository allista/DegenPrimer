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
Uninterruptable MultiProcessing

Created on Mar 2, 2014

@author: Allis Tauri <allista@gmail.com>
'''

import signal
import multiprocessing as mp
from multiprocessing.managers import BaseManager 

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