% [m Bins sd stderr] = MeanSmooth(x,y,Bins)
% splits x into bins, and calculates the mean of y in each bin.
%
% Bins may be a single number - in which case the range of x is divided into
% that many equally sized bins, or it may be a vector in which case
% the nth bin is all those points greater than or equal to x(n) but less than
% x(n+1), and the last one is all those greater than x(end)
%
% third and fourth optional output arguments are the sd and standard error for
% the relevant bin.
%
% see also MedianSmooth

function [m, Bins, sd, stderr] = MeanSmooth(x,y,Bins)

if (length(Bins) == 1)
	nBins = Bins;
	Bins = min(x) + (0:nBins-1)*range(x)/nBins;
else	
	nBins = length(Bins);
end

m = zeros(nBins,1);
sd = zeros(nBins,1);
stderr = zeros(nBins,1);

for n=1:nBins-1
	PointsInBin = find(x>=Bins(n) & x<Bins(n+1));
	if ~isempty(PointsInBin)
		m(n) = mean(y(PointsInBin));
		sd(n) = std(y(PointsInBin));
		stderr(n) = sd(n) / sqrt(length(PointsInBin));
	end
end

% do last bin
PointsInBin = find(x>=Bins(nBins));
if ~isempty(PointsInBin)
	m(nBins) = mean(y(PointsInBin));
	sd(nBins) = std(y(PointsInBin));
	stderr(nBins) = sd(nBins) / sqrt(length(PointsInBin));
end

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu