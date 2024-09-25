% [m xBins yBins iq] = MedianSmooth2(x,y,z,xBins, yBins)
% splits the xy plane into bins, and calculates the mean of z in each bin
%
% xBins may be a single number - in which case the range of x is divided into
% that many equally sized bins, or it may be a vector in which case
% the nth bin is all those points greater than or equal to x(n) but less than
% x(n+1), and the last one is all those greater than x(end)
%
% fourth optional output arguments is the iqr for
% the relevant bin.
%
% see also MeanSmooth2, MedianSmooth

function [m, xBins, yBins, iq] = MedianSmooth(x,y,z,xBins,yBins)

% replace xBins and yBins if you need to
if (length(xBins) == 1)
	nxBins = xBins;
	xBins = min(x) + (0:nxBins-1)*range(x)/nxBins;
else	
	nxBins = length(xBins);
end

if (length(yBins) == 1)
	nyBins = yBins;
	yBins = min(y) + (0:nyBins-1)*range(y)/nyBins;
else	
	nyBins = length(yBins);
end

% z may be multidimensional
nDim = size(z,2);

m = zeros(nxBins,nyBins,nDim);
iq = zeros(nxBins,nyBins,nDim);


for nx=1:nxBins, 
	for ny=1:nyBins
		% find points in bin
		if (nx<nxBins) 	
			xIn = (x>=xBins(nx) & x<xBins(nx+1));
		else
			xIn = (x>=xBins(nx));
		end
	
		if (ny<nyBins) 	
			yIn = (y>=yBins(ny) & y<yBins(ny+1));
		else
			yIn = (y>=yBins(ny));
		end
	
		PointsInBin = find(xIn & yIn);
		if ~isempty(PointsInBin)
			m(nx,ny,:) = median(z(PointsInBin,:));
			iq(nx,ny,:) = iqr(z(PointsInBin,:));
		end
	end
end



% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu