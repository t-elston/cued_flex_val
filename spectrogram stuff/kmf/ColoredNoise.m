% x = ColoredNoise(n, Csd, f)
%
% simulates n samples of colored gaussian noise with csd given by (Csd, f).
% Assumes sample rate is 2*max(f)
%
% NB psd is POWER cross-spectrum, not amplitude spectrum.
%
% NB THIS DOES NOT WORK FOR VECTOR TIME SERIES YET.
% But it's open-source software!  If you make it 
% do that, email me, harris@axon.rutgers.edu

function x = ColoredNoise(n, Csd, f)

nDims = size(Csd, 2);
if size(Csd,3) ~= nDims
	error('Csd must have dims 2 and 3 the same size');
end

FiltOrd = length(f)*2 - 1;
Nyquist = max(f);

% asd is square root of Csd (amplitude spectrum)
asd = zeros(size(Csd));
for fb = 1:length(f)
	asd(fb,:,:) = squeeze(Csd(fb,:,:))^ .5;
end

if nDims==1
	b = fir2(FiltOrd, f/Nyquist, sqrt(Csd));
else
for ch1=1:nDims
	for ch2=1:nDims
		b(:,ch1,ch2) = fir2(FiltOrd, f/Nyquist, asd(:,ch1,ch2));
	end
end
end
Gaussian = randn(n+FiltOrd,nDims);

% have to do a loop to avoid MATLAB filtering only works 1d at a time

Filtered = zeros(size(Gaussian));

for col=1:nDims
	Filtered = Filtered + fftfilt(b(:,:,col),Gaussian);
end

x = Filtered(FiltOrd+1:end,:);

	
[ct f] = mtcsd(x, FiltOrd, 2*Nyquist, FiltOrd, 0, 1);
%[ct f] = psd(x, FiltOrd, 2*Nyquist, FiltOrd, 0, 1);
%for ch1=1:nDims
%	for ch2=1:nDims
%		subplot(nDims,nDims,(ch1-1)*nDims+ch2)
%		plot(f,20*log10(abs([Csd(:,ch1,ch2) ct(:,ch1,ch2)])));
%	end
%end

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu