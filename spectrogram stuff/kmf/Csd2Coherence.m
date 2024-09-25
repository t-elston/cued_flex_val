% Coherence = Csd2Coherence(Csd)
%
% A simple function that takes an array of cross-spectra
% and converts it to (complex) coherences.
%
% The input array may be a frequency or a time-frequency 
% transform, and is indexed like Csd(...,Ch1,Ch2)
% Output array has the same indexing.

function Coherence = Csd2Coherence(Csd)

nDims = ndims(Csd);
OriginalSize = size(Csd);
nCh = size(Csd,nDims-1);

if (size(Csd,nDims) ~= nCh) 
	error('Last 2 dimensions of input array must be the same size');
end

if (nDims==2)
	NewSize = [1 OriginalSize];
else
	NewSize = [ prod(OriginalSize(1:nDims-2)) , nCh, nCh];
end

ReshapedCsd = reshape(Csd, NewSize);

for Ch1 = 1:nCh, for Ch2 = 1:nCh
	ReshapedCoherence(:,Ch1,Ch2) = ReshapedCsd(:,Ch1,Ch2) ...
			./ sqrt(ReshapedCsd(:,Ch1,Ch1) .* ReshapedCsd(:,Ch2,Ch2));
end,end

Coherence = reshape(ReshapedCoherence, OriginalSize);

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu