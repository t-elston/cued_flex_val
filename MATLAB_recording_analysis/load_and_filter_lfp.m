function [lfp_signal] = load_and_filter_lfp(fname, channel, notch_filts)

[fs,~,~,~,lfp_signal] = plx_ad(fname,channel);

% for f = 1:length(notch_filts)
%     lfp_signal = fftfilt(notch_filts{f}.Coefficients,lfp_signal);
% end



end