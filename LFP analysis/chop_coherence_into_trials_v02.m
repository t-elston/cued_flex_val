function [trial_data, coh_ts, z_coh] = chop_coherence_into_trials_v02(coh, ts, pl2_event_times, freqs)


% what's the new sampling frequency of the coherence?
coh_fs = round(mean(diff(ts)));

neg_offset = 2000;
pos_offset = 2000;

% pull out coherence data aligned to pics on
choice_times = pl2_event_times.pics_on;
n_trials = numel(choice_times);

coh_ts = neg_offset*-1 : coh_fs : pos_offset;


% cycle through the trials
for t = 1:n_trials
          
        % find the closest trial start/stop times in the spectrogram timestamps
        [~,ts_start] = min(abs(ts - (choice_times(t)-neg_offset)));
        [~,ts_end] = min(abs(ts - (choice_times(t)+pos_offset))); 
                
        trial_data(:,:,t) = coh(:,ts_start: ts_end);
        
end % of cycling through trials

mean_coh = mean(trial_data,3);
z_coh = zscore(mean_coh, [], 'all');

% % plot for sanity checking
% figure;
% [t,f] = meshgrid(coh_ts,freqs);
% coh = surf(t,f,z_coh);
% view(0,90);
% ylabel('Freq (Hz)');
% xlabel('Time from Choice');
% title('Mean HPC-OFC Coherence')
% shading interp
% cb = colorbar;
% % cb.FontSize = 12;
% % cb.Position = [.25 .66 .01 .2];
% caxis([0, 5])
% axis tight
% set(gca,'FontSize',12,'LineWidth',1);
% xlim([-1000,1000])
% colormap(hot)



end % of function