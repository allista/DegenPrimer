#!/usr/bin/python
# coding=utf-8
#
# Copyright (C) 2012 Allis Tauri <allista@gmail.com>
# 
# degen_primer is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# indicator_gddccontrol is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

'''
Created on Jun 30, 2012

@author: Allis Tauri <allista@gmail.com>
'''

import string
import random
from time import ctime

text_width = 80

def print_exception(e):
    print '\n%s: %s\n' % (type(e).__name__, e.message or 'no message.')
    
    
def random_text(length):
    return ''.join(random.choice(string.ascii_uppercase + string.digits) for i in range(length))
#end def


def hr(s, symbol='-'):
    left_hr  = (text_width - len(s))/2
    right_hr = (text_width - len(s) - left_hr)
    return symbol*left_hr + s + symbol*right_hr + '\n\n'
#end def

def time_hr(symbol='#'):
    return hr(' %s ' % ctime(), symbol)


def wrap_text(_text):
    wrapped_text = ''
    strings = _text.split('\n')
    for s in strings:
        while len(s) > text_width:
            new_line      = s[:text_width]
            last_space    = new_line.rfind(' ')
            wrapped_text += new_line[:last_space+1]+'\n'
            s = s[last_space+1:]
        wrapped_text += s+'\n'
    return wrapped_text[:-1]
#end def


#tests
if __name__ == '__main__':
    text_width = 40
    sample_text = 'The quick broun fox jumped over a very, very lazy dog forty two times in a row!'
    print wrap_text(sample_text)