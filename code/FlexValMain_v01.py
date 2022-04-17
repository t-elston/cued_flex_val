# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#%% imports

from scipy.ndimage.filters import uniform_filter1d
import CuedFlexValBehavior as bhv


#%% analysis
# make sure data are in the path
#datadir = ("C:/Users/Thomas Elston/Documents/PYTHON/CuedFlexVal/bhvdata/singleday/")
datadir = ("C:/Users/Thomas Elston/Documents/PYTHON/CuedFlexVal/bhvdata/4cues8stim/")





alldata = bhv.getMLH5(datadir,run_debug=False)

modelresults,choicegrid = bhv.assesschoice(alldata,run_debug=True)

bhv.choice_or_rt_by_val(alldata,debug_ = False)