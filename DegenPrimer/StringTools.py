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
# degen_primer is distributed in the hope that it will be useful, but
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

from math import log
import random
import string
from textwrap import TextWrapper
from time import ctime


text_width = 80
wrapper = TextWrapper(width=text_width, expand_tabs=False)

def print_exception(e):
    print('\n%s: error code %s; %s\n' % 
          (type(e).__name__,
           hasattr(e, 'errno') and str(e.errno) or 'N/A',  
           e.message or 'no message.'))
#end def
    
    
def random_text(length):
    return ''.join(random.choice(string.ascii_uppercase + string.digits) 
                   for _unused in xrange(length))
#end def


def hr(s='', symbol='-', width=text_width):
    return s.center(width, symbol) + '\n\n'
#end def

def time_hr(symbol='#'):
    return hr(' %s ' % ctime(), symbol)


def wrap_text(text, width=None):
    #replace intermediate whitespaces with single spaces
#    text = ' '.join(text.split())
    #wrap text
    if width:
        old_width = wrapper.width
        wrapper.width = width
        wrapped_text  = wrapper.fill(text)
        wrapper.width = old_width
    else: wrapped_text = wrapper.fill(text)
    return wrapped_text+'\n'
#end def


def line_by_line(texts, widths, divider='|', filler=' ', j_down=False, j_center=False):
    if len(texts) != len(widths):
        raise ValueError('StringTools.line_by_line: for each text a line width '
                         'should be provided')
    _texts = [wrap_text(t, widths[w]-1)[:-1].splitlines() for w,t in enumerate(texts)]
    _lines = max(len(t) for t in _texts)
    output_text = ''
    for l in xrange(_lines):
        for w,t in enumerate(_texts):
            if not j_down:
                if l < len(t):
                    if j_center:
                        output_text += t[l].center(widths[w]-1, filler) + divider
                    else: 
                        output_text += t[l].ljust(widths[w]-1, filler) + divider
                else: output_text += filler*(widths[w]-1) + divider
            else:
                l_start = _lines-len(t)
                if l >= l_start:
                    if j_center:
                        output_text += t[l-l_start].center(widths[w]-1, filler) + divider
                    else: 
                        output_text += t[l-l_start].ljust(widths[w]-1, filler) + divider
                else: output_text += filler*(widths[w]-1) + divider
        output_text += '\n'
    return output_text
#end if


def print_dict(_dict, header=None, delimiter=':'):
    if not _dict: return ''
    keys   = _dict.keys() + [header[0],] if header else _dict.keys()
    values = _dict.values() + [header[1],] if header else _dict.values()
    max_key_len   = max(len(str(key)) for key in keys)
    max_value_len = max(len(str(val)) for val in values)
    dict_items    = _dict.items()
    dict_items.sort(key=lambda x: x[0])
    if header: dict_items.insert(0, header)
    dict_str = ''
    for key, val in dict_items:
        dict_str += '%s%s %s%s %s\n' % (str(key),
                                        ' '*(max_key_len-len(str(key))), 
                                        delimiter,
                                        ' '*(max_value_len-len(str(val))),
                                        str(val))
    return dict_str
#end def


def print_table(table, delimiter=':'):
    if not table: return ''
    if len(set(len(row) for row in table)) > 1: 
        raise ValueError('StringTools.print_table: all rows in a table should be of equal size.')
    max_col_len = [max(len(str(table[r][c])) for r in xrange(len(table))) for c in xrange(len(table[0]))]
    table_str = ''
    for row in table:
        for c in xrange(len(row)):
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
    if quantity is None or quantity < 0: return 'N/A  %s' % unit
    if quantity == 0: return '0.0  %s' % unit
    mag = -1*log(quantity, 10)
    if       mag < 0:  return '%.1f  %s' % (quantity,      unit)
    if 0  <= mag < 3:  return '%.1f m%s' % (quantity*1e3,  unit)
    if 3  <= mag < 6:  return '%.1f u%s' % (quantity*1e6,  unit)
    if 6  <= mag < 9:  return '%.1f n%s' % (quantity*1e9,  unit)
    if 9  <= mag < 12: return '%.1f p%s' % (quantity*1e12, unit)
    if 12 <= mag < 15: return '%.1f f%s' % (quantity*1e15, unit)
    if 15 <= mag:      return '%.1f a%s' % (quantity*1e18, unit)
    #if everything fails for some unknown reason
    return 'N/A  %s' % unit
#end def