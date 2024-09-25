% Co = Comodugram2(x, nFFT, SampleRate, FreqRange, NW, Detrend)
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

 
function Co = Comodugram2(x, nFFT, SampleRate, FreqRange, NW, Detrend)

if (nargin<2 | isempty(nFFT)) nFFT = 256; end;
if (nargin<3 | isempty(SampleRate)) SampleRate = 2; end;
if (nargin<4 | isempty(FreqRange)) FreqRange = [0 SampleRate/2]; end;
if (nargin<5) NW = 3; end;
if (nargin<6) Detrend = []; end;


Clip = 0;
% if Clip is 1, correlations below 0 will be replaced by 0.


% check to see if the first 2 arguments are (x, nFFT) or (psg, f)
if all(size(nFFT) == [1 1])
	nChannels = size(x,2);
	nSamples = size(x,1);
	% compute number of time bins that will be produced by specgram(), and allocate array
	%nTimeBins = fix((nSamples - nFFT/2) / (nFFT/2));
	nTimeBins = round((nSamples-nFFT) / (nFFT/2));
%	Spex = zeros(1 + nFFT/2, nTimeBins, nChannels);
	% calculate spectrograms
	for i=1:nChannels
		[Spex(:,:,i), f, t] = mtcsg(x(:,i), nFFT, SampleRate);
	end;
else
	Spex = x;
	f = nFFT;
	nChannels = size(Spex,3);
	nTimeBins = size(Spex,2);
end

Spex = abs(Spex);
Means = mean(Spex,2);
%Means = std(Spex,0,2);

% normalize them
SpexNorm = Spex ./ repmat(Means, [1 nTimeBins 1]);

% find frequency bins to consider
FreqBins = find(f >= FreqRange(1) & f <= FreqRange(2));
nFreqBins = length(FreqBins);

% calculate covariances
DataMat = reshape(permute(SpexNorm(FreqBins, :, :), [2 1 3]), ...
				[nTimeBins, nFreqBins*nChannels]);
				
CorrMat = cov(DataMat);

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