% ImageMatrix(C, [gamma]) or ImageMatrix(X,Y,C, [gamma])
%
% Takes a 4D array C and displays a matrix of imagesc plots
%
% input arguments are as in imagesc - the only difference is that C is 4D
% the first 2 dimensions of C are X and Y the second 2 are X and Y of the subplots
% (i think).
%
% if C is complex it will call ComplexImage - this is where [gamma] comes in.

function ImageMatrix(varargin)

% sort out input arguments (BORING!) 
if nargin == 1
	AxesSpecified = 0;
	C = varargin{1};
	gamma = 1;
elseif nargin == 2
	AxesSpecified = 0;
	C = varargin{1};
	gamma = varargin{2};
elseif nargin == 3
	AxesSpecified = 1;
	X = varargin{1};
	Y = varargin{2};
	C = varargin{3};
	gamma = 1;
elseif nargin == 4
	AxesSpecified = 1;
	X = varargin{1};
	Y = varargin{2};
	C = varargin{3};
	gamma = varargin{4};
else 
	error('number of arguments needs to be 1 to 4');
end

nPlotsX = size(C,3);
nPlotsY = size(C,4);

% now make the plot matrix
for i=1:nPlotsX
	for j=1:nPlotsY
	
		% select correct subplot
		subplot(nPlotsY, nPlotsX, i + (j-1)*nPlotsX);
		% make individual plot
		ThisPlot = C(:,:,i,j);
		if AxesSpecified
			if (isreal(ThisPlot))
				imagesc(X,Y,ThisPlot);
			else
				ComplexImage(X,Y,ThisPlot,gamma);
			end
		else
			if (isreal(ThisPlot))
				imagesc(ThisPlot);
			else
				ComplexImage(ThisPlot,gamma);
            end
        end
        axis xy
        ylim([3.5 21.5]);
        colorbar
	end
end


% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu