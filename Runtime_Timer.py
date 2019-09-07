# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 18:00:28 2018

@author: ark47
"""
import timeit
class Timer:
    def __init__(self):
        self.start = 0
       
    def Start():
        Timer.start = timeit.default_timer()
        
    def End():
        stop = timeit.default_timer()
        seconds = stop - Timer.start
        m, s = divmod(seconds, 60)
        h, m = divmod(m, 60)
        print ("Runtime %d:%02d:%02d" % (h, m, s))

