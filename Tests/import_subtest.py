'''
Created on Mar 13, 2014

@author: allis
'''
from DegenPrimer.TD_Functions import PCR_P
import DegenPrimer.TD_Functions as tdf

def f():
    print PCR_P
    print tdf.PCR_P
    
from multiprocessing.managers import BaseManager

class FManager(BaseManager): pass
FManager.register('f', f)

fm = FManager(); fm.start()


    