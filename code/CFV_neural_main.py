'''
CFV_neural_main.py

Top of stack for neural analysis
'''
#%% import modules
import importlib
import CFV_neural_functions as neur

#%% reload modules after updating them
importlib.reload(neur)

#%% where are the data?
datadir = 'C:/Users/Thomas Elston/Documents/PYTHON/CuedFlexVal/recdata/'

# %% find selective units
cue_ttests, choice_glms, all_bhv = neur.unit_selectivity_analysis(datadir)


# %% plot cue selective units
neur.plot_cue_selective_units(cue_ttests)

# %%
neur.plot_cue_pop_summary(cue_ttests)

# %% let's make a psuedo population decoder

all_train, all_test, train_state_labels, test_state_labels, ts, all_ofc = neur.make_cue_psuedo_pop(datadir, 500)
neur.decode_cue_psuedo_pop(all_train, all_test, train_state_labels, test_state_labels, ts, all_ofc)

# %%
neur.plot_choice_selective_units(choice_glms)
# %%

#all_train, all_test, train_state_labels, test_state_labels, ts, all_ofc = neur.make_choice_psuedo_pop(datadir, 500)
all_train, all_test, train_val_labels, test_val_labels, ts, all_ofc = neur.make_choice_value_psuedo_pop(datadir, 500)
neur.decode_choice_val_psuedo_pop(all_train, all_test, train_val_labels, test_val_labels, ts, all_ofc)
# %%

