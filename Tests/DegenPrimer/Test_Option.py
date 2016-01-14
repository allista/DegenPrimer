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

from DegenPrimer.Option import Option

if __name__ == '__main__':
    opt = Option(name='optimization_parameter', 
                desc='PCR parameter to optimize. Optimization is performed '
                'in terms of maximization of an objective function which is: '
                'sum_over_target_products(product_concentration*product_fraction^3).',
                nargs='*',
                py_type=None,
                options=(Option(name='name',
                                   desc='Name of the parameter. '
                                   'Any of PCR conditions as well as '
                                   'primer and polymerase concentrations may be used. '
                                   'To optimize primer concentration, '
                                   'set this to the ID of the primer.',
                                   nargs=1,
                                   py_type=str,
                                   ),
                            Option(name='min',
                                   desc='Lower boundary of optimization interval.',
                                   nargs=1,
                                   py_type=str,
                                   ),
                            Option(name='ini',
                                   desc='Initial value of the parameter. Optimization '
                                   'starts from this value and tries to find the nearest '
                                   'optimum. Thus changing it may change the outcome of '
                                   'the optimization.',
                                   nargs=1,
                                   py_type=str,
                                   ),
                            Option(name='max',
                                   desc='Upper boundary of optimization interval',
                                   nargs=1,
                                   py_type=str,
                                   ),
                            ))
    print vars(opt)