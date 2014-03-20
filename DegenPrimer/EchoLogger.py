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

import logging, sys, traceback
from StringTools import print_exception
from Output import OutIntercepter

class EchoLogger(OutIntercepter):
    '''Wrapper around logging module to capture stdout-err into a log file
    while still print it to std'''

    def __init__(self, name, level=logging.INFO):
        OutIntercepter.__init__(self)
        self._name      = name
        self._log       = name+'.log'
        self._level     = level
        self._logger    = logging.getLogger(name)
        self._handler   = logging.FileHandler(self._log, encoding='UTF-8')
        self._formatter = logging.Formatter('[%(asctime)s] %(message)s')
        #set log level
        self._handler.setLevel(self._level)
        self._logger.setLevel(self._level)
        #assemble pipeline
        self._handler.setFormatter(self._formatter)
        self._logger.addHandler(self._handler)
    #end def
        
    def __del__(self):
        self._handler.close()
        
        
    def __enter__(self):
        OutIntercepter.__enter__(self)
        self._logger.log(self._level, '=== START LOGGING ===')
        return self
    #end def
    
    def __exit__(self, _type, _value, _traceback):
        if _type is not None:
            print_exception(_value)
            self._logger.error('Exception occured:', 
                               exc_info=(_type, _value, _traceback))
            traceback.print_exception(_type, _value, _traceback, file=self._olderr)
        sys.stdout = self._oldout
        sys.stderr = self._olderr
        self._logger.log(self._level, '=== END LOGGING ===')
        return True
    #end def
    
        
    def write(self, text):
        log_text = ' '.join(text.split())
        if log_text: self._logger.log(self._level, log_text, exc_info=False)
        self._oldout.write(text)
    #end def
    
    def flush(self): self._oldout.flush()
#end class



#tests
if __name__ == '__main__':
    with EchoLogger('test'):
        print 'message 1'
        print 'message 2'
        print 'message 3'
        raise Warning('Tadam!')
        print 'message 4'
        