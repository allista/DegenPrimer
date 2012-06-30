#!/usr/bin/python
# coding=utf-8

# Copyright (C) 2012 Allis Tauri <allista@gmail.com>
# 
# indicator_gddccontrol is free software: you can redistribute it and/or modify it
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

text_width = 80


def hr(string, symbol='-'):
    global text_width
    left_hr  = (text_width - len(string))/2
    right_hr = (text_width - len(string) - left_hr)
    return symbol*left_hr + string + symbol*right_hr + '\n\n'
#end def


def wrap_text(text):
    global text_width
    wrapped_text = ''
    strings = text.split('\n')
    for string in strings:
        while len(string) > text_width:
            new_line   = string[:text_width]
            last_space = new_line.rfind(' ')
            wrapped_text += new_line[:last_space+1]+'\n'
            string = string[last_space+1:]
        wrapped_text += string+'\n'
    return wrapped_text[:-1]
#end def

if __name__ == '__main__':
    text_width = 40
    text = 'The quick broun fox jumped over a very, very lazy dog forty two times in a row!'
    print wrap_text(text)