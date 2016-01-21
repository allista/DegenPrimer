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
Created on Jan 17, 2016

@author: allis

A quick hack to run tests before implementing proper UnitTesting
'''

import os

if __name__ == '__main__':
    dirname = 'DegenPrimer'
    sources = [f for f in os.listdir(dirname) if f.endswith('.py')]
    os.chdir(dirname)
    for s in sources:
        print 'Running %s\n' %s
        os.system('python %s' % s)
        print '='*80+'\n'
    print 'Done'