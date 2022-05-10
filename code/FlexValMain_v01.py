# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#%% imports
import CuedFlexValBehavior as bhv


#%% analysis
# make sure data are in the path
datadir = ("C:/Users/Thomas Elston/Documents/PYTHON/CuedFlexVal/bhvdata/Don_6633/main data/")
datadir = ("C:/Users/Thomas Elston/Documents/PYTHON/CuedFlexVal/bhvdata/4cues8stim/")

alldata = bhv.getMLH5(datadir,run_debug=False)

modelresults,choicegrid = bhv.assesschoice(alldata,run_debug=False)

bhv.choice_or_rt_by_val(alldata,debug_ = False)

bhv.plot_rt_by_value(alldata, debug_= False)

