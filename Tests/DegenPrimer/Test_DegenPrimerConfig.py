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

from DegenPrimer.DegenPrimerConfig import DegenPrimerConfig

def test():
    class TestConfig(DegenPrimerConfig):
        def _override_option(self, option):
            return DegenPrimerConfig._override_option(self, option)

    conf = TestConfig()
    conf.parse_configuration('../F-TGAM_0057-268_d1-R-TGAM_0055-624-d4.cfg')
    opts = conf.options
    print conf.find_option('primers/primer/sequence').full_name
    conf.from_options(opts)
    conf.job_id = 'test'
    conf.save_configuration()
