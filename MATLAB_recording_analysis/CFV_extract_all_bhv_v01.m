function [all_bhv] = CFV_extract_all_bhv_v01(bhvdir)

all_bhv = table;

bhv_files = dir([bhvdir '*.bhv2']);

n_files = numel(bhv_files);

for f = 1:n_files 
    
    [bhv, ~] = extractCFVbhv_vRecorded(bhvdir, bhv_files(f).name);
    
    all_bhv = [all_bhv; bhv];
    
end % of looping over files


all_bhv(all_bhv.UseTrial==0,:) = [];

end% of function