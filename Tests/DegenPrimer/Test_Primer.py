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

from DegenPrimer.Primer import Primer
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

if __name__ == '__main__':
#    import cProfile
#    primer = Primer(SeqRecord(Seq('ATARTCTYCGAMGGCTATKCAGNCTGGGANGGNTACGNGGGTAAANAAACG'),id='primer1'), 0.9e-6)
    primer = Primer(SeqRecord(Seq('ATARTCTYCGAMGGCNATKCAGGNCTGRGGA'),id='primer1'), 0.9e-6)
    primer.generate_components()
    primer.calculate_Tms()
    print primer.str_sequences
    print primer.concentration
    print repr(primer)
    
#    from timeit import timeit
#    from tests.violin_plot import violin_plot
#    from matplotlib.pyplot import figure, show
#    
#    gen = [timeit('primer.generate_components()', 'from __main__ import primer', number=1) for _i in xrange(1)]
#    gen_mp = [timeit('primer.generate_components_mp()', 'from __main__ import primer', number=1) for _i in xrange(1)]
#    data = [gen, gen_mp]
#    print data
#    
#    fig=figure()
#    ax = fig.add_subplot(111)
#    violin_plot(ax,data,range(len(data)),bp=1)
#    show()
    
    
#    cProfile.run('''
#for _n in xrange(100): primer.generate_components_mp()
#for _n in xrange(100): primer.generate_components()
#''',
#    'gen_components.profile')
    print 'Done.'
    
