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
Created on 2016-01-14

@author: Allis Tauri <allista@gmail.com>
'''

from DegenPrimer.StringTools import wrap_text

if __name__ == '__main__':
    txt = '''If true, TextWrapper attempts to detect sentence endings and ensure 
    that sentences are always separated by exactly two spaces. This is generally 
    desired for text in a monospaced font. However, the sentence detection 
    algorithm is imperfect: it assumes that a sentence ending consists of a 
    lowercase letter followed by one of '.', '!', or '?', possibly followed by 
    one of '"' or "'", followed by a space. One problem with this is algorithm 
    is that it is unable to detect the difference between “Dr.” in'''
    print wrap_text(txt)
    print '='*80
    