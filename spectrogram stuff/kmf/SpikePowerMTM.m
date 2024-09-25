% PSpikes = SpikePowerMTM(Spikes, NW, NFFT)
%
% A second quick script to calculate the power spectra in a bunch of vectors
% which may or may not be individual spike waveforms.  This time it uses the
% multitaper method
%
% Data is nSamples by nSpikes
% Power is nSamples/2 by nSpikes because you only have nSamples/2
% frequency bands for real waveforms
%
% At the moment this program works by FFT with a Hamming window.

function PSpikes = SpikePowerMTM(Spikes, NW, NFFT)

nSamples = size(Spikes, 1);
nSpikes = size(Spikes, 2);

PSpikes = zeros(NFFT/2+1, nSpikes);

for i=1:nSpikes
	PSpikes(:,i) = pmtm(Spikes(:,i), NW, NFFT);
end;



% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu