% ComplexImage(C [,gamma]) or ComplexImage(X,Y,C [,gamma])
%
% does an hsv colormap of a complex array
% gamma is a power for the brightness - like on your monitor
function ComplexImage(varargin)

switch nargin
	case 1 
		AxesSpecified = 0;
		C = varargin{1};
		gamma = 1;
	case 2
		AxesSpecified = 0;
		C = varargin{1};
		gamma = varargin{2};		
	case 3
		AxesSpecified = 1;
		X = varargin{1};
		Y = varargin{2};
		C = varargin{3};
		gamma = 1;
	case 4
		AxesSpecified = 1;
		X = varargin{1};
		Y = varargin{2};
		C = varargin{3};
		gamma = varargin{4};
	otherwise	
		error('number of arguments needs to be 1 to 4');
end

Amp = abs(C);
AmpMin = min(Amp(:));
AmpMax = max(Amp(:));
Phase = angle(C);

Hsv = zeros([size(C) 1]);



Hsv(:,:,1) = (Phase+pi)/(2*pi);
Hsv(:,:,2) = 1;

if (AmpMax == AmpMin)
	Hsv(:,:,3) = 1;
else
	Hsv(:,:,3) = ((Amp-AmpMin)/(AmpMax-AmpMin)).^gamma;
end

if (AxesSpecified)
	image(X, Y, hsv2rgb(Hsv));
else
	image(hsv2rgb(Hsv));
end

tit = sprintf('Amplitude range %f to %f', AmpMin, AmpMax);
title(tit)

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu