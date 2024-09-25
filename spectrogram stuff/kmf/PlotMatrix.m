% PlotMatrix([x], Y, ...)
%
% Takes an array Y and displays a matrix of plots
%
% If Y is 3D then the first index is the X axis, the
% second and third are the subplot number
%
% If Y is 4D then it will plot multiple plots on each graph.
% The first index is the X axis, second is line number, third
% and fourth are subplot number
%
% final arguments are passed to plot.

function PlotMatrix(varargin)

% sort out input arguments (BORING!)
FirstStringArg = 1+nargin;
for i=1:nargin
	if(isstr(varargin{i}))
		FirstStringArg = i;
		break;
	end
end

if FirstStringArg == 2;
	AxesSpecified = 0;
	Y = varargin{1};
else
	AxesSpecified = 1;
	X = varargin{1};
	Y = varargin{2};
end

if (ndims(Y) == 3)
	nPlotsX = size(Y,2);
	nPlotsY = size(Y,3);
elseif (ndims(Y) == 4)
	nPlotsX = size(Y,3);
	nPlotsY = size(Y,4);
else
	error('Input Y must have 3 or 4 dimensions');
end

% now make the plot matrix
for i=1:nPlotsX
	for j=1:nPlotsY
	
		% select correct subplot
		subplot(nPlotsY, nPlotsX, i + (j-1)*nPlotsX);
		
		if (ndims(Y) == 3) 
			ThisGraph = Y(:,i,j);
		elseif (ndims(Y) ==4)
			ThisGraph = Y(:,:,i,j);
		else
			error('Input Y must have 3 or 4 dimensions');
		end
		
		% make individual plot
		if AxesSpecified
			plot(X,ThisGraph,varargin{FirstStringArg:end});
		else
			plot(ThisGraph,varargin{FirstStringArg:end});
		end
		axis tight
	end
end


% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu