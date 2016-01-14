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
Created on Jan 14, 2016

@author: allis
'''

import os, datetime

header = r"""# coding=utf-8
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
Created on %s

@author: Allis Tauri <allista@gmail.com>
'''

"""

if __name__ == '__main__':
    dirname = 'DegenPrimer'
    sources = [f for f in os.listdir(os.path.join('..', dirname)) if f.endswith('.py')]
    try: os.mkdir(dirname)
    except: pass
    for s in sources:
        test = ''
        with open(os.path.join('..', dirname, s)) as inp:
            for l in inp:
                if l.find('__main__') >=0: 
                    print '%s has __main__ clause' % s
                if test: test += l
                elif l.startswith('#test'):
                    test += header % datetime.date.today()
        if not test: continue
        test_name = os.path.join(dirname, 'Test_'+s)
        if os.path.isfile(test_name): continue 
        with open(test_name, 'w') as out:
            print 'Saving %s' % test_name
            out.write(test)
    print 'Done'