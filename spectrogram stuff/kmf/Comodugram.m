% [Co nTimeBins] = Comodugram(x, nFFT, SampleRate, FreqRange, NW, Detrend)
%
% Takes an input sequence x and does a multi-pane
% plot showing correlated changes in power in frequency
% bands 
%
% nFFT and SampleRate are just like for specgram() function
%
% FreqRange = [fLow fHigh] allows you to view only a certain 
% frequency range.  Specify it in Hz.
%
% NB x is of the form x(Time, Channel)   !!!IT CHANGED!!!
%
% NW is an argument for the multitaper method, as is Detrend
%
% optional output argument nTimeBins gives number of points in the regression
% - so you can do significance testing.
 
function [Co, nTimeBins] = Comodugram(x, nFFT, SampleRate, FreqRange, NW, Detrend)

if (nargin<2 | isempty(nFFT)) nFFT = 256; end;
if (nargin<3 | isempty(SampleRate)) SampleRate = 2; end;
if (nargin<4 | isempty(FreqRange)) FreqRange = [0 SampleRate/2]; end;
if (nargin<5) NW = 3; end;
if (nargin<6) Detrend = []; end;

Clip = 0;
% if Clip is 1, correlations below 0 will be replaced by 0.

nChannels = size(x,2);
nSamples = size(x,1);

% compute number of time bins that will be produced by specgram(), and allocate array
%nTimeBins = fix((nSamples - nFFT/2) / (nFFT/2));
nTimeBins = round((nSamples-nFFT) / (nFFT/2));
spex = zeros(1 + nFFT/2, nTimeBins, nChannels);

% calculate spectrograms
for i=1:nChannels
%	[spex(:,:,i), f, t] = mtcsg(x(:,i), nFFT, SampleRate);
	[spex(:,:,i), f, t] = mtchgnorm(x(:,i), nFFT, SampleRate,[],[],NW,Detrend);
end;
spex = abs(spex);

% find frequency bins to consider
FreqBins = find(f >= FreqRange(1) & f <= FreqRange(2));
nFreqBins = length(FreqBins);

% calculate correlation coefficients
DataMat = reshape(permute(spex(FreqBins, :, :), [2 1 3]), ...
				[nTimeBins, nFreqBins*nChannels]);
				
CorrMat = corrcoef(DataMat);
if (Clip) CorrMat = clip(CorrMat, 0, 1); end;

% produce output array and plot(if required)
C = zeros(nFreqBins,nFreqBins,nChannels,nChannels);

for i=1:nChannels
	for j=1:nChannels
		
		C(:,:,i,j) = CorrMat((i-1)*nFreqBins + (1:nFreqBins), (j-1)*nFreqBins + (1:nFreqBins));
		
		if (nargout<1)
			subplot(nChannels, nChannels, j + (i-1) * nChannels);
			imagesc(f(FreqBins), f(FreqBins), C(:,:,i,j));
		end;
		
		drawnow;
	end;
end;

if nargout >=1
	Co = C;
end

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu