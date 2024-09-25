%% FOOOF Matlab Wrapper Example - Plot a FOOOF model
%
% This example computes an example power spectrum model fit, 
% and then plots the result.
%

%% Run Example

% Load data
load('data/ch_dat_one.mat');

% Calculate a power spectrum with Welch's method
[psd, freqs] = pwelch(ChLFP.Values, 1000, [], [], 1000);

% Transpose, to make inputs row vectors
freqs = freqs';
psd = psd';

% FOOOF settings
settings = struct();  % Use defaults
f_range = [1, 50];

% Run FOOOF, also returning the model
fooof_results = fooof(freqs, psd, f_range, settings, true);

% Plot the resulting model
fooof_plot(fooof_results)