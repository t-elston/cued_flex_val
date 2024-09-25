% Comoherogram(x, nFFT, SampleRate, FreqRange, Clip)
%
% This function is like Comodugram but it also looks
% at coherences.

function Comoherogram(x, nFFT, SampleRate, FreqRange)

if (nargin<2 | isempty(nFFT)) nFFT = 256; end;
if (nargin<3 | isempty(SampleRate)) SampleRate = 2; end;
if (nargin<4 | isempty(FreqRange)) FreqRange = [0 SampleRate/2]; end;
if (nargin<5) Clip = 0; end;

nChannels = size(x,2);
nSamples = size(x,1);

% calculate spectra and coherences
[spex, f, t] = mtcoherenorm(x, nFFT, SampleRate);
nTimeBins = length(t);


% Make array of indices to pick out power spectra then coherences
DiagIndices = find(eye(nChannels)) ; 
UpperIndices = find(triu(ones(nChannels))-eye(nChannels));
Indices = [DiagIndices ; UpperIndices];
nSpectra = length(Indices);

% find frequency bins to consider
FreqBins = find(f >= FreqRange(1) & f <= FreqRange(2));
nFreqBins = length(FreqBins);

% calculate correlation coefficients
DataMat = reshape(permute(spex(FreqBins, :, Indices), [2 1 3]), ...
				[nTimeBins, nFreqBins*nSpectra]);
				
CorrMat = corrcoef(DataMat);
if (Clip) CorrMat = clip(CorrMat, 0, 1); end;

% plot it

for i=1:nSpectra
	for j=1:nSpectra
		subplot(nSpectra, nSpectra, j + (i-1) * nSpectra);
		
		imagesc(f(FreqBins), f(FreqBins), ...
				CorrMat((i-1)*nFreqBins + (1:nFreqBins), (j-1)*nFreqBins + (1:nFreqBins)));
		
		drawnow;
	end;
end;


% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu