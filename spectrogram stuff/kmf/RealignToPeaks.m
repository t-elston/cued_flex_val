% RealignToPeaks(X, BeforeSamps, AfterSamps)
%
% Realigns the waveforms in X so that their peaks are all at the same point
% Time goes down the columns of X.  all other indices are treated as referencing different waves.
%
% BeforeSamps and AfterSamps give the number of samples before and after the peak to
% include in the output

function Out = RealignToPeaks(X, BeforeSamps, AfterSamps)

% Get original size
OriginalSize = size(X);
nWaves = prod(OriginalSize(2:end));
OriginalWaveLength = OriginalSize(1);

% reshape input into a 2d array
Waves = reshape(X, OriginalSize(1), nWaves);

% find peak positions
[PeakVal PeakPos] = max(Waves);

% Set up array of realigned waves
RealignedWaves = zeros(1+BeforeSamps+AfterSamps, nWaves);

% Fill it up
for i=1:nWaves

	% check whether entire requested area is actually there
	if (PeakPos(i)<=BeforeSamps)
	
		% peak position is too soon...
		warning(sprintf('wave %d was padded at the start with %d zeros', i, 1+BeforeSamps-PeakPos(i)));
		RealignedWaves(2+BeforeSamps-PeakPos(i):end, i) =Waves(1:PeakPos(i)+AfterSamps, i);
	
	elseif (PeakPos(i)>OriginalWaveLength-AfterSamps)
	
		% peak position is too late ....
		warning(sprintf('wave %d was padded at the end with %d zeros', i, PeakPos(i)-OriginalWaveLength+AfterSamps));
		RealignedWaves(1:1+OriginalWaveLength-PeakPos(i)+BeforeSamps, i) =Waves(PeakPos(i)-BeforeSamps:OriginalWaveLength, i);
		
	else
		% we have enough samples on both sides
		RealignedWaves(:,i) = Waves(PeakPos(i)-BeforeSamps:PeakPos(i)+AfterSamps, i);
	end;
end;

% Reshape output

Out = reshape(RealignedWaves, [1+BeforeSamps+AfterSamps, OriginalSize(2:end)]);

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu