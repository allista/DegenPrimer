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
Created on Jul 19, 2012

@author: Allis Tauri <allista@gmail.com>
'''


from math import log, exp
import StringTools 

window_percent = 0.05
precision      = 1e6

def _format_electrophoresis(phoresis, window):
    text_width  = StringTools.text_width
    max_line    = max(l[1] for l in phoresis)
    max_mark    = max(len(str(l[2])) for l in phoresis)*2
    line_width  = text_width - max_mark - 7 #mark b :###   :
    #format phoresis
    phoresis_text = ''
    phoresis_text += ' '*(max_mark+5)+':'+'-'*line_width+':\n'
    phoresis_text += ' '*(max_mark+5)+':'+' '*line_width+':\n'
    for l in range(len(phoresis)):
        line = phoresis[l]
        next_mark   = str(phoresis[l+1][2]) if l < len(phoresis)-1 else str(line[2]+window) 
        mark_spacer = max_mark - len(str(line[2])) - len(next_mark)
        line_value  = (line_width*line[1])/max_line
        line_spacer = line_width - line_value 
        phoresis_text += '%s-%d%s bp :%s%s:\n' % (next_mark,
                                                  line[2], 
                                                  ' '*mark_spacer,
                                                  '#'*line_value,
                                                  ' '*line_spacer) 
    phoresis_text += ' '*(max_mark+5)+':'+' '*line_width+':\n'
    phoresis_text += ' '*(max_mark+5)+':'+'-'*line_width+':\n'
    return phoresis_text
#end def


def print_electrophoresis(products):
    global precision, window_percent
    max_len      = max(p['len'] for p in products)
    window       = int(max_len*window_percent)
    nearest_srip = max_len-window
    max_len_log  = int(log(max_len)*precision)
    min_len_log  = int(log(min(p['len'] for p in products))*precision)
    window_log   = int((log(max_len)-log(nearest_srip))*precision)
    #construct phoresis
    phoresis    = [[l,0,int(exp(l/precision))] for l in range(min_len_log, max_len_log+window_log, window_log)]
    for product in products:
        l = int(log(product['len'])*precision)
        p = (l - min_len_log)/window_log
        phoresis[p][1] += product['quantity']
    #format phoresis
    return _format_electrophoresis(phoresis, window)
#end def