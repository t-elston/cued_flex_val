# -*- coding: utf-8 -*-
"""
CFV_preprocess_main.py
@author: Thomas Elston

Top of the stack for assessing behavior and neurophysiology
in a context-dependent valuation task.
"""
#%% imports
import cued_flex_val_behavior as bhv
import read_and_wrangle_pl2_spikedata as get_spk
import importlib

#%% update modules
importlib.reload(bhv)
importlib.reload(get_spk)


#%% define directories and extract behavior
raw_bhv_dir = 'C:/Users/Thomas Elston/Documents/PYTHON/CuedFlexVal/recordings/Don/raw_bhv/'
spk_dir ='C:/Users/Thomas Elston/Documents/PYTHON/CuedFlexVal/recordings/Don/sorted/'
save_bhv_dir = 'C:/Users/Thomas Elston/Documents/PYTHON/CuedFlexVal/recordings/Don/preprocessed_bhv/'
save_spk_dir = 'C:/Users/Thomas Elston/Documents/PYTHON/CuedFlexVal/recordings/Don/preprocessed_spk/'

#%% extract behavior
all_bhv = bhv.getMLH5(raw_bhv_dir, save_dir = save_bhv_dir)

#%% extract spikes
get_spk.make_spike_table(spk_dir, 2, 40,38, 100, 25, offsets = [1000,1000], save_dir = save_spk_dir)


#%% summarize behavior
modelresults,choicegrid = bhv.assesschoice(all_bhv)




# %%
