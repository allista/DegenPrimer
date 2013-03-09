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
Created on Jun 19, 2012

@author: Allis Tauri <allista@gmail.com>
'''

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
import os
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

from DegenPrimer.UnifiedNN import UnifiedNN

from distutils.core import setup
setup(name='degen-primer',
      version='2.2-rc1',
      description='Degenerate primer characterization and '
                  'in silico PCR simulation software.',
      long_description=read('README'),
      license='GPL-3',
      author='Allis Tauri',
      author_email='allista@gmail.com',
      url='https://launchpad.net/degenprimer',
      classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX',
        'Programming Language :: Python'],
      packages=['DegenPrimer'],
      scripts=['degen_primer', 
               'align_as_amino'],
      data_files=[('share/degen_primer', [UnifiedNN.internal_NN_filename,
                                          UnifiedNN.dangling_ends_filename,
                                          UnifiedNN.loops_filename,
                                          UnifiedNN.tri_tetra_hairpin_loops_filename])]
      )
