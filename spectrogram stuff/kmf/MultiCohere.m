% MultiCohere(x, nFFT, SampleRate, FreqRange)
%
% Takes an input sequence x and does a multi-pane
% plot showing power spectra on the diagonal
% and coherences off the diagonal
%
% nFFT and SampleRate are just like for cohere() function
%
% FreqRange = [fLow fHigh] allows you to view only a certain 
% frequency range.  Specify it in Hz.
%
% NB x is of the form x(Time, Channel)

function MultiCohere(x, nFFT, window, overlap, SampleRate, FreqRange)

if (nargin<2 | isempty(nFFT)) nFFT = 256; end;
if (nargin<3 | isempty(SampleRate)) SampleRate = 2; end;
if (nargin<4 | isempty(FreqRange)) FreqRange = [0 SampleRate/2]; end;


nChannels = size(x,2);

for i=1:nChannels
	
	% plot psd	
	subplot(nChannels, nChannels, i + (i-1)*nChannels);
	
	
	% coherences
	for j=1:nChannels
		subplot(nChannels, nChannels, i + (j-1)*nChannels);
		if (i==j)		
			pwelch(x(:,i), window, overlap, nFFT, SampleRate);
			set(gca, 'xlim', FreqRange);
			ylabel('Power');
			drawnow;
		else
			mscohere(x(:,i), x(:,j), window, overlap, nFFT, SampleRate);
			set(gca, 'xlim', FreqRange);
			ylabel('Coherence');
			drawnow;
		end;
	end
end;

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu