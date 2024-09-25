% Power = SpikePower(Spikes, nFFT)
%
% A quick script to calculate the power spectra in a bunch of vectors
% which may or may not be individual spike waveforms
%
% Data is nSamples by nSpikes
% Power is nSamples/2 by nSpikes because you only have nSamples/2
% frequency bands for real waveforms
%
% At the moment this program works by FFT with a Hamming window.
%
% nFFT is an optional argument for zero padding.

function Power = SpikePower(Spikes, nFFT)

nSamples = size(Spikes, 1);
nSpikes = size(Spikes, 2);

% subtract mean
zSpikes = Spikes - repmat(mean(Spikes,1), nSamples, 1);
%wSpikes = zSpikes .*repmat(hamming(nSamples),1,nSpikes);
wSpikes = Spikes .*repmat(hamming(nSamples),1,nSpikes);

% do we need to zero pad?
if (nargin>=2 & nFFT > nSamples)
	wSpikes = [ wSpikes ; zeros(nFFT-nSamples, nSpikes) ];
else
	nFFT = nSamples;
end;
Power = abs(fft(wSpikes)).^2;
Power = Power(1:nFFT/2+1,:);

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu