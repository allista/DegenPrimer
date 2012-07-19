#!/usr/bin/python
# coding=utf-8
#
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

import string
import random
from time import ctime

text_width = 80

def print_exception(e):
    print '\nException occurred: %s\n%s\n' % (str(type(e)), e.message)
    
    
def random_text(length):
    text = ''
    for x in range(length):
        text += random.choice(string.ascii_uppercase + string.digits)
    return text 
#end def


def hr(string, symbol='-'):
    global text_width
    left_hr  = (text_width - len(string))/2
    right_hr = (text_width - len(string) - left_hr)
    return symbol*left_hr + string + symbol*right_hr + '\n\n'
#end def

def time_hr(symbol='#'):
    return hr(' %s ' % ctime())


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


def format_histogram(histogram, col_names, hist_width=10):
    global text_width
    histogram_string = ''
    #maximum column name width
    max_name = text_width-hist_width-2
    #cut name column title if necessary 
    if len(col_names[0]) > max_name:
        names = col_names[0][:max_name]
    else: names = col_names[0]
    #spacers
    names_spacer = max_name - len(names)
    hists_spacer = hist_width - len(col_names[1])
    #column_titles 
    histogram_string += '-'*(names_spacer/2)+names + \
                        '-'*(names_spacer-names_spacer/2) + '|' + \
                        '-'*(hists_spacer/2)+col_names[1] + \
                        '-'*(hists_spacer-hists_spacer/2) + '|\n'
    #histogram lines
    for col in histogram:
        if len(col[0]) > max_name:
            name = col[0][:max_name]
        else: name = col[0]
        name_spacer = max_name - len(name)
        hist_value  = int(col[1]*hist_width)
        col_spacer  = hist_width - hist_value
        histogram_string += name + ' '*name_spacer
        histogram_string += ':' + '#'*hist_value + ' '*col_spacer + ':' + '\n'
    histogram_string += '\n'
    return histogram_string
#end def


#tests
if __name__ == '__main__':
    text_width = 40
    text = 'The quick broun fox jumped over a very, very lazy dog forty two times in a row!'
    print wrap_text(text)