'''
Created on Mar 13, 2014

@author: allis
'''

from DegenPrimer.TD_Functions import PCR_P
import DegenPrimer.TD_Functions as tdf
from import_subtest import f, fm

if __name__ == '__main__':
    print PCR_P
    print tdf.PCR_P
    
    PCR_P.PCR_T = 53
    
    print PCR_P
    print tdf.PCR_P
    
    f()
    fm.f()