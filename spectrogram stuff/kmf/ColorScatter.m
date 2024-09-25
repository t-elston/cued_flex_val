% ColorScatter(Data, Colors, Size, ColorOrder)
%
% draws a 2d scatter plot of the data with colors allocated according to the vector Colors
% in size Size (default value 5)
%
% ColorOrder is an optional vector specifying what the colors are
% which overrides the default
%
% optional output h returns a vector of handles for the colors

function h = ColorScatter(Data, Colors, Size, ColorOrder)

if (nargin<3) Size = 5; end;

MinColor = min(Colors);
MaxColor = max(Colors);

if (nargin<4) ColorOrder = get(gca, 'ColorOrder'); end;

if (size(ColorOrder,1)==1) ColorOrder=ColorOrder(:); end;

ih = ishold;
h = [];
for c=MinColor:MaxColor
	Pts = find(Colors==c);
	%Col =ColorOrder(1+c-MinColor,:);
	%Col = 'r';
	
	%if(size(ColorOrder,2) == 3)
	Col = ColorOrder(1+c-MinColor,:);
	%else Col = ColorOrder(1+c-MinColor);
	%end;
	
	hh = plot(Data(Pts,1), Data(Pts,2), 'o', 'MarkerSize', Size, ...
		'MarkerEdgeColor', Col, 'MarkerFaceColor', Col);
	h = [h ; hh];
	hold on
end;

%legend(num2str((MinColor:MaxColor)'));
	
if (ih==0) hold off; end;


% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu