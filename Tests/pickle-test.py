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
Created on Mar 2, 2014

@author: Allis Tauri <allista@gmail.com>
'''

from BioUtils.Tools.tmpStorage import tmpDict, roDict
import shelve, os

if __name__ == '__main__':
    single = '../tests/DP-PCR-single'
    multi  = 'DP-PCR-multi'
    
    from pympler.asizeof import asizeof
    def mem(obj):
        return asizeof(obj)/1024.0/1024.0
    
    
#    with profile.timestamp('load_single'):
#    sd = tmpDict(single)
#    mixture = sd.itervalues().next()
#    sd.close()

#    print mem(mixture)        
#    print len(mixture.annealings), mem(mixture.annealings)
#    print len(mixture.templates[1]), mem(mixture.templates[1])
#    print len(mixture.templates[0]), mem(mixture.templates[0])
#    print len(mixture.products), mem(mixture.products)
#        
##    with profile.timestamp('save_single'):
##        sd = tmpDict(single)
##        sd['result'] = mixture
##        sd.close()
#    
#    with profile.timestamp('save_multiple'):
#        md = shelve.open(multi, protocol=-1)
#        md['reaction_id'] = mixture.reaction_id
#        md['fwd_templates'] = mixture.templates[1]
#        md['rev_templates'] = mixture.templates[0]
##        for t_id, template in enumerate(mixture.templates[1]):
##            md['tf%d'%t_id] = template
##        for t_id, template in enumerate(mixture.templates[0]):
##            md['tr%d'%t_id] = template
#        for p_id, product in mixture.products.iteritems():
#            md['p%d'%p_id] = product
#        for i, annealing in enumerate(mixture.annealings):
#            md['a%d'%i] = annealing
#        md.close()
    
    with profile.timestamp('open_multiple'):
        md = tmpDict(multi)
#    with profile.timestamp('read_keys'):
#        for k in md.iterkeys(): pass
    with profile.timestamp('load_templates'):
        templates = [[],[]]
        templates[1] = md['fwd_templates']
        templates[0] = md['rev_templates']
    with profile.timestamp('load_multiple'):
        products    = dict()
        reaction_id = md['reaction_id']
        for k in md.iterkeys():
            if k.startswith('p'):
                products[int(k.strip('p'))] = md[k]
#            elif k.startswith('tf'):
#                templates[1].append(v)
#            elif k.startswith('tr'):
#                templates[0].append(v)
            elif k.startswith('a'):
                #do something with annealing
                pass
    md.close()

    print 'Done'