'''
Created on Nov 23, 2012

@author: allis
'''

import sys, os
import argparse
from time import sleep, time
from datetime import timedelta
import errno
import cProfile
from DegenPrimer.StringTools import print_exception
    
    
class PrimerAction(argparse.Action):
    def __init__(self, **kwargs):
        argparse.Action.__init__(self, **kwargs)
        
    def __call__(self, parser, namespace, values, option_string=None):
        primer = dict()
        primer['sequence'] = values[0]
        primer['id'] = values[1]
        primer['concentration'] = values[2]
        primers = getattr(namespace, self.dest)
        if primers is None:
            setattr(namespace, self.dest, [primer])
        else: primers.append(primer)
#end class
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('foo', default='bar')
    parser.add_argument('--primer', metavar='*', dest='primers', action=PrimerAction(), nargs=3)
    args = parser.parse_args('erty --primer ATGC id 5 --primer GGCA id1 9'.split())
    print args
    
    args = parser.parse_args('qwer'.split())
    print args
    
#    import matplotlib.pyplot as plt
#    import numpy as np
#
#    # evenly sampled time at 200ms intervals
#    t = np.arange(0., 5., 0.2)
#    # red dashes, blue squares and green triangles
#    plt.plot(t, t, 'r--', t, t**2, 'bs', t, t**3, 'g^')
#
#    plt.ylabel('some numbers')
#    
#    plt.show()

    
    
    #cProfile.run('pass', 'test.profile')
    #print 'Done'
                