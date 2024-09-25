%   [CrossingTimes, Waves, AfterSpikes, BeforeSpikes, NoSpikes] = ThresholdCrossings(In, Threshold, SpikeThreshold, UpDown);
%
% Finds a matrix of spikes which correspond to times that the waveform In crossed over the
% Threshold from below
%
%  CrossingTimes is a vector of times that the waves crossed this potential
%  AfterSpike is a list of those waveforms which later went on to produce a spike
%  as defined by exceeding SpikeThreshold.  BeforeSpikes and NoSpikes similarly.
%
%  UpDown is an optional argument that you set to 1 if you want to count up and
%  down crossings, not just upward ones.

function [CrossingTimes, Waves, AfterSpikes, BeforeSpikes, NoSpikes] = ThresholdCrossings(In, Threshold, SpikeThreshold, UpDown);

if (nargin<4) UpDown = 0; end;

% Find crossing times
fprintf('Finding crossing times...\n');

nCrossings = 0;

if (UpDown)
	for i=2:max(size(In))
		if ((In(i-1)<Threshold & In(i) >=Threshold) | (In(i-1)>Threshold & In(i) <=Threshold))
			nCrossings = nCrossings+1;
			CrossingTimes(nCrossings) =i;
		end;
	end;
else
	for i=2:max(size(In))
		if (In(i-1)<Threshold & In(i) >=Threshold)
			nCrossings = nCrossings+1;
			CrossingTimes(nCrossings) =i;
		end;
	end;
end;

% Check there are some ....
if (nCrossings == 0) error('No Crossings found!'); end;

% make sure CrossingTimes is a column
CrossingTimes = CrossingTimes(:);

% Now we have a vector of data and list of times to extract from it
% Use some dirty matlab tricks to get it individual spikes out of the data
fprintf('Extracting Spikes...\n');

BeforeSamps = 200;
AfterSamps = 200;
SampleRange = -BeforeSamps:AfterSamps;
nSamples = size(SampleRange, 2);

% Matrix of indices

Indices = CrossingTimes(:, ones(1, nSamples)) + SampleRange(ones(1,nCrossings),:);

Waves = In(Indices);

% Calculate those that lead to a spike and those that came from a spike

AfterSpikes =find(max(Waves(:,1+BeforeSamps:end),[], 2) > SpikeThreshold);
BeforeSpikes =find(max(Waves(:,1:BeforeSamps), [], 2) > SpikeThreshold);
NoSpikes =setdiff(1:nCrossings,union(BeforeSpikes,AfterSpikes));
% plot stuff

figure(101);
if ~isempty(BeforeSpikes) plot(Waves(BeforeSpikes,:)', 'g'); end;
hold on
if ~isempty(AfterSpikes) plot(Waves(AfterSpikes,:)', 'b'); end;
plot(Waves(NoSpikes, :)', 'r')
hold off


% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu