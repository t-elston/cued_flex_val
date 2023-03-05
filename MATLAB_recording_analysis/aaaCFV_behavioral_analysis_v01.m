% CFV_behavioral_analysis_v01

bhvdir = 'E:\CuedFlexVal\bhv2_files\';

% get data from the individual files
[all_bhv] = CFV_extract_all_bhv_v01(bhvdir);


% now plot responses as a function of value difference
plotCFVchoices_v01(all_bhv);
