function [down_FRs, down_zFRs, down_ts, mean_FR] = downsample_spikes(raw_firing_rates, win_size, step_size, raw_ts, total_n_units) 


% make array of window start/ends for downsampling
win_centers = win_size/2:step_size:numel(raw_firing_rates(1,:,1));
win_starts = win_centers-round(win_size/2);
win_starts(win_starts==0)=1;
win_ends   = win_centers+round(win_size/2);
win_ends(win_ends > numel(raw_firing_rates(1,:,1))) = numel(raw_firing_rates(1,:,1));

down_FRs = [];
down_ts = [];
down_zFRs=[];
mean_FR = [];

for u = 1:total_n_units
    
    for w = 1:numel(win_centers)
        down_FRs(:,w,u) =  nanmean(raw_firing_rates(:,win_starts(w) : win_ends(w),u),2);
        down_ts(w) = nanmean(raw_ts(win_starts(w) : win_ends(w)));
    end % of cylcing over win_centers
    
    % zscore the firing rates
    u_mean = nanmean(down_FRs(:,:,u),'all');
    u_std  = nanstd(down_FRs(:,:,u),[],'all');
    
    down_zFRs(:,:,u) = (down_FRs(:,:,u) - u_mean) / u_std;  
    mean_FR(u) = u_mean;
     
end % of downsampling firing rates

end % of function