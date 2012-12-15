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
from math import log

text_width = 80

def print_exception(e):
    print '\n%s: %s\n' % (type(e).__name__, e.message or 'no message.')
    
    
def random_text(length):
    return ''.join(random.choice(string.ascii_uppercase + string.digits) for i in range(length))
#end def


def hr(s, symbol='-'):
    return s.center(text_width, symbol) + '\n\n'
#end def

def time_hr(symbol='#'):
    return hr(' %s ' % ctime(), symbol)


def wrap_text(_text):
    wrapped_text = ''
    lines = _text.splitlines()
    for line in lines:
        while len(line) > text_width:
            new_line      = line[:text_width]
            last_space    = new_line.rfind(' ')
            wrapped_text += new_line[:last_space+1]+'\n'
            line = line[last_space+1:]
        wrapped_text += line+'\n'
    return wrapped_text
#end def


def print_dict(_dict, delimiter=':'):
    if not _dict: return ''
    max_key_len   = max(len(key) for key in _dict.keys())
    max_value_len = max(len(str(val)) for val in _dict.values())
    dict_str = ''
    for key in _dict:
        val = _dict[key]
        dict_str += '%s%s %s%s %s\n' % (key,
                                        max_key_len-len(key), 
                                        delimiter,
                                        max_value_len-len(str(val)),
                                        str(val))
    return dict_str
#end def


def print_table(table, delimiter=':'):
    if not table: return ''
    if len(set(len(row) for row in table)) > 1: 
        raise ValueError('StringTools.print_table: all rows in a table should be of equal size.')
    max_col_len = [max(len(str(table[r][c])) for r in range(len(table))) for c in range(len(table[0]))]
    table_str = ''
    for row in table:
        for c in range(len(row)):
            if c == 0: #first column left-aligned
                table_str += str(row[c]) #data
                if len(row) > 1:
                    table_str += ' '*(max_col_len[c]-len(row[c])) #spacer
            else: #others right-aligned
                table_str += ' '+delimiter
                table_str += ' '*(max_col_len[c]-len(row[c])+1) #spacer
                table_str += str(row[c]) #data
        table_str += '\n'
    return table_str
#end def


def format_quantity(quantity, unit='U'):
    '''Given quantity in units, return it's string representation using 
    prefixes m, u, n, p, etc.'''
    if not quantity: return '0.0  %s' % unit
    if quantity < 0: return 'N/A  %s' % unit
    mag = -1*log(quantity, 10)
    if 0  <  max < 1:  return '%.1f  %s' % (quantity,      unit)
    if 1  <= mag < 3:  return '%.1f m%s' % (quantity*1e3,  unit)
    if 3  <= mag < 6:  return '%.1f u%s' % (quantity*1e6,  unit)
    if 6  <= mag < 9:  return '%.1f n%s' % (quantity*1e9,  unit)
    if 9  <= mag < 12: return '%.1f p%s' % (quantity*1e12, unit)
    if 12 <= mag < 15: return '%.1f f%s' % (quantity*1e15, unit)
    if 15 <= mag:      return '%.1f a%s' % (quantity*1e18, unit)
#end def