# -*- coding: utf-8 -*-
"""
FlexValMain

This is the top of the stack for assessing behavior in a context-dependent 
valuation task
"""
#%% imports
import cued_flex_val_behavior as bhv


#%% analysis
# make sure data are in the path
datadir = ("C:/Users/Thomas Elston/Documents/PYTHON/CuedFlexVal/bhvdata/4cues8stim/")
datadir2 = ("C:/Users/Thomas Elston/Documents/PYTHON/CuedFlexVal/bhvdata/Don_6633/main data/")
#datadir = ("C:/Users/Thomas Elston/Documents/PYTHON/CuedFlexVal/bhvdata/Don_6633/single day/")

alldata = bhv.getMLH5(datadir)
overrep_data = bhv.getMLH5(datadir2)


modelresults,choicegrid = bhv.assesschoice(alldata)

bhv.choice_or_rt_by_val_v02(alldata)

bhv.choice_or_rt_by_val_v01(alldata)

bhv.plot_rt_by_value(alldata)

bhv.check_switch_cost(alldata)

bhv.compare_choice_rt_sacc1_x_conditions(alldata,overrep_data)


# %%
