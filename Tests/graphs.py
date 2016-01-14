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
Created on Feb 12, 2014

@author: Allis Tauri <allista@gmail.com>
'''

import numpy as np
from scipy.sparse.csgraph import connected_components, csgraph_from_dense

if __name__ == '__main__':
    
    def _matrix_from_edges(edges, vertices=None):
        #extract vertices from edges if not provided
        if not vertices:
            vertices = set()
            for edge in edges:
                vertices.add(edge[0])
                vertices.add(edge[1])
            vertices = sorted(vertices)
        #construct dense graph
        N = max(vertices)+1
        graph = [[0]*N for _n in xrange(N)]
        for edge in edges: 
            graph[edge[0]][edge[1]] = 1
        return np.array(graph)
    #end def
    
    def _subgroups(edges):
        graph = _matrix_from_edges(edges)
        n_comp, comp = connected_components(csgraph_from_dense(graph), directed=False)
        if n_comp < 2: return [edges]
        groups = [[] for _n in xrange(n_comp)]
        for edge in edges: groups[comp[edge[0]]].append(edge)
        groups = [group for group in groups if group]
        return groups
    #end def
    
    edges = [[1,2], [3,4], [5,6], [1,4], [3,4], [5,2], [7,8], [8,9], [7,10]]
    
    graph = _matrix_from_edges(edges)
    print graph
    
    components = connected_components(csgraph_from_dense(graph), directed=False)
    print components
    
    groups = _subgroups(edges)
    print np.array(groups)
    
    
#    from timeit import timeit
#    from tests.violin_plot import violin_plot
#    from matplotlib.pyplot import figure, show
#    
#    lfor = [timeit('''
#A = 0
#for edge in edges:
#    A += edge[0]*edge[1]
#    ''', 'from __main__ import edges') for _i in xrange(100)]
#    
#    lsum = [timeit('sum(edge[0]*edge[1] for edge in edges)', 'from __main__ import edges') for _i in xrange(100)]
#    
#    fig=figure()
#    ax = fig.add_subplot(111)
#    data = [lfor, lsum]
#    violin_plot(ax,data,range(len(data)),bp=1)
#    show()
    
    print 'Done'
