# -*- coding: utf-8 -*-
"""
random utility functions

@author: Thomas Elston
"""
import numpy as np

def makehistbinwidths(binwidth,data):
    bins=np.arange(min(data), max(data) + binwidth, binwidth)
    return bins

