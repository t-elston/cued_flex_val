function filtered_signal = bandpass_filter(input_signal, fs, f_low, f_high, filter_order)
% BANDPASS_FILTER Bandpass filter using Butterworth filter
%
%   filtered_signal = bandpass_filter(input_signal, fs, f_low, f_high, filter_order)
%
%   input_signal: Input signal to be filtered
%   fs: Sampling frequency of the input signal (Hz)
%   f_low: Lower cutoff frequency of the bandpass filter (Hz)
%   f_high: Upper cutoff frequency of the bandpass filter (Hz)
%   filter_order: Order of the Butterworth filter (integer)

% Normalize the cutoff frequencies
Wn = [f_low, f_high] / (fs / 2);

% Design the Butterworth bandpass filter
[b, a] = butter(filter_order, Wn, 'bandpass');

% Apply zero-phase filtering using filtfilt
filtered_signal = filtfilt(b, a, input_signal);
end