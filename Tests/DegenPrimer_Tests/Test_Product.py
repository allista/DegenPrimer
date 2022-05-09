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
Created on Jan 22, 2016

@author: allis
'''

from DegenPrimer.SecStructures import reverse_complement, Duplex
from DegenPrimer.Product import Product

def test():
    fwd_primer = 'ATARTCTYCGAMGGCTATCC'
    rev_primer = 'NAAGGYTTAGAKGCGGAAG'
    d1 = Duplex(fwd_primer, reverse_complement(fwd_primer))
    d2 = Duplex(rev_primer, reverse_complement(rev_primer))
    f = 'Product_test.json.gz'
    p = Product('test_template', 50, 150, (d1,), (d2,))
    f = p.save(f)
    p1 = Product.load(f)
    print p
    print '-'*80
    print p1
    assert str(p) == str(p1)

if __name__ == '__main__':
    test()