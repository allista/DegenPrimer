'''
Created on Jan 25, 2016

@author: allis
'''

import os
import numpy as np
import pandas as p
from scipy.stats import linregress
from matplotlib import pyplot as plt


if __name__ == '__main__':
    datadir = 'data'
    filenames = [f for f in os.listdir(datadir) if f.endswith('.sol')]
    dfs = [p.read_csv(os.path.join(datadir, f)) for f in filenames]
    plt.figure()
    s = []
    for df in dfs:
        df = df.sort_values('r')
        df['x'] = df.K
        m = np.max(df.x)
        if m > 1:
            df['x'] /= np.max(df.x)
        lr = linregress(df.x, df.r)
        plt.plot(df.r, df.x, 'o-')
        print 'sum(Err^2)=%f' % np.sum((df.x-df.r)**2) 
        print lr
        print
        
#    s.sort(key=lambda x: x[1])
#    lr = linregress([p[0] for p in s], [p[1] for p in s])
#    print lr
#    plt.plot([p[0] for p in s], [p[1] for p in s], 'o')
    plt.xlabel('r')
    plt.ylabel('normed X')
    plt.show()