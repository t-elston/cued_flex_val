% [m Bins iq Resid] = MedianSmooth(x,y,Bins)
% splits x into bins, and calculates the median of y in each bin.
%
% Bins may be a single number - in which case the range of x is divided into
% that many equally sized bins, or it may be a vector in which case
% the nth bin is all those points greater than or equal to x(n) but less than
% x(n+1), and the last one is all those greater than x(end)
%
% third optional output argument is the iqr of the relevant bin
%
% fourth optional output is the residual
%
% see also MeanSmooth

function [m, Bins, iq, Resid] = MedianSmooth(x,y,Bins)

if (length(Bins) == 1)
	nBins = Bins;
	Bins = min(x) + (0:nBins-1)*range(x)/nBins;
else	
	nBins = length(Bins);
end

m = zeros(nBins,1);
iq = zeros(nBins,1);
PointsInBin = cell(nBins, 1);

for n=1:nBins-1
	PointsInBin{n} = find(x>=Bins(n) & x<Bins(n+1));
	if ~isempty(PointsInBin{n})
		m(n) = median(y(PointsInBin{n}));
		iq(n) = iqr(y(PointsInBin{n}));
	end
end

% do last bin
PointsInBin{nBins} = find(x>=Bins(nBins));
if ~isempty(PointsInBin{nBins})
	m(nBins) = median(y(PointsInBin{nBins}));
	iq(nBins) = iqr(y(PointsInBin{nBins}));
end

% Calculate residuals
Resid = zeros(size(x));
for n=1:nBins
	Resid(PointsInBin{n}) = y(PointsInBin{n}) - m(n);
end

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu