# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 15:38:14 2017

@author: erwan
"""

from __future__ import print_function
from __future__ import absolute_import
import sys
from time import time
from six.moves import range

class ProgressBar():
    ''' A console progress-bar
    
    Example
    -------
    
    add a progress bar in a loop::
    
        pb = ProgressBar(N)
        for i, ... in enumerate(...):
            (...)
            pb.update(i)
        pb.done()
    
    See test in progress_bar.py
    
    Notes
    -----
    
    Todo:
    
    - One day extend for multiprocss with several progress values? 
      https://stackoverflow.com/questions/7392779/is-it-possible-to-print-a-string-at-a-certain-screen-position-inside-idle
    
    '''
    
    def __init__(self, N, active=True):
        ''' 
        write to progress bar completion status i/N. 
        
        
        Parameters    
        ----------
        
        N: int
            total number of iterations
        
        '''
        self.t0 = time()
        self.N = N
        self.active = active
        
    def set_active(self, active=True):
        ''' Option to activate/deactivate the ProgressBar. Used not to make it 
        appear on small processes (based on a condition) without changing most
        of the code'''

        self.active = active

    def update(self, i):
        ''' 
        write to progress bar completion status i/N. If t0 is not None, also 
        write the time spent
        
        
        Parameters    
        ----------
        
        i: int
            current iteration
        '''
        
        if not self.active: return
        
        N = self.N
        t0 = self.t0
        
        if t0 is None:
            msg = '{0:.1f}%'.format(i/N*100)
        else:
            msg = '({0:.0f}s)\t{1:.1f}%'.format(time()-t0, i/N*100)
            
        if sys.stdout is not None:
            sys.stdout.write('\r'+msg)
            sys.stdout.flush()
            
    def done(self):
        
        if not self.active: return
        
        self.update(self.N)
        # make new line
        if sys.stdout is not None:
            sys.stdout.write('\n')
            sys.stdout.flush()
            
        
# %% Test     
def test_progress_bar(*args, **kwargs):
    ''' Minimal example of a progress bar '''
    
    from radis.misc.progress_bar import ProgressBar
    from time import sleep
    from numpy.random import rand
    
    print('Testing progress bar')
    
    a = 0
    r = list(range(1000))
    N = len(r)
    pb = ProgressBar(N, t0=time())
    for i in r:
        if i % 10: 
            pb.update(i)
        a += i
        sleep(rand()*3e-3)
    pb.done()
    
    return True  # nothing implemented

if __name__ == '__main__':
    test_progress_bar()