% [y f] = Harmonix(x, nFFT, fS, WinLength, nOverlap, nHarms)
%
% Code to detect harmonics in an input time series.
% x - input 1d time series.  nFFT and fS as in psd
%
% output y is a 2D array with the first dimension being fundamental frequency
% (as indexed by f) and the second dimension being harmonic number.


function [y, f] = Harmonix(varargin);


% default argument BS
[x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers,nChannels,nSamples,...
	nFFTChunks,winstep,select,nFreqBins,f,t] = mtparam(varargin);

if (nargin<6) 
	nHarms = 2;
else 
	nHarms = NW;	% renamed so we can reuse the mtparam function - that's all
end
	
Fundamentals = 1:floor(nFFT/(2*nHarms));

HarmCo = zeros(length(Fundamentals),nHarms);
Power = zeros(nFFT,1);

Taper = hamming(WinLength);

for j=1:nFFTChunks
	Segment = x((j-1)*winstep+[1:WinLength]) .* Taper;
	
	Periodogram = fft(Segment,nFFT);
	Periodogram(find(Periodogram==0)) = eps;
	Magnitude = abs(Periodogram);
	Phase = Periodogram ./ Magnitude;

	for HarmNo = 1:nHarms;	
		hRange = 1+HarmNo*(Fundamentals-1);
	
		HarmCo(:,HarmNo) = HarmCo(:,HarmNo) + Phase(Fundamentals).^HarmNo .* conj(Phase(1+HarmNo*(Fundamentals-1))) ...
					.* Magnitude(Fundamentals) .* Magnitude(1+HarmNo*(Fundamentals-1));
		
	end
	
	Power = Power + abs(Periodogram).^2;

end

Power = Power / nFFTChunks;

for HarmNo = 1:nHarms

	hRange = 1+HarmNo*(Fundamentals-1);
	y(:,HarmNo) = HarmCo(:,HarmNo) ./ (nFFTChunks*sqrt(Power(Fundamentals) .* Power(1+HarmNo*(Fundamentals-1))));
	
end
	

if nargout==0
	plot(f(1:length(Fundamentals)),abs(y(:,2:end)));
	grid on;
	legends = (['1st'; '2nd'; '3rd'; '4th'; '5th'; '6th'; '7th'; '8th'; '9th'])     ;
	if (nHarms>2) legend(legends(2:nHarms,:)); end;
	clear y;
end

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu