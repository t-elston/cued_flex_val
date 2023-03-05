function [coh_z_map, raw_coh_map] = CFV_make_coh_map_v01(trial_data,iti_data)

trial_data = circshift(trial_data,3,1);
iti_data = circshift(iti_data,3,1);


iti_means = squeeze(mean(iti_data,2));
trial_means = squeeze(mean(trial_data,2));


[n_freq,n_times,n_trials] = size(trial_data);
coh_z_map = NaN(n_freq,n_times);

raw_coh_map = mean(trial_data,3);

for f = 1:n_freq
    
    for t = 1:n_times
        
              
       [pval,~,stats]=signrank(squeeze(trial_data(f,t,:)),iti_means(f,:)');
       %[~,pval,~,stats]=ttest(squeeze(trial_data(f,t,:)),iti_means(f,:)');


        
       if pval < .001
        coh_z_map(f,t)   =stats.zval;
       else
           coh_z_map(f,t)=0;
       end
           
    end % of cycling through time steps
      
end % of cycling through frequencies


end