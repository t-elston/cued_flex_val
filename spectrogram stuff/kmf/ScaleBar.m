% h = ScaleBar(Left,Bottom,Width,Height,xlabel,ylabel)
%
% Draws a scale bar on the current plot
% bottom left corner is (x0,y0)
% size is (Width, Height)
% labels given by xlabel, ylabel

function h = ScaleBar(Left,Bottom,Width,Height,xlabel,ylabel)

HoldState = ishold;

hold on

% Plot x line
hx = plot([Left, Left+Width], [Bottom, Bottom], 'k-');

% x label
hxl = text(Left + Width/2, Bottom, xlabel, ...
		'HorizontalAlignment', 'center', ...
		'VerticalAlignment', 'top');
		
% Plot y line
hy = plot([Left, Left], [Bottom, Bottom+Height], 'k-');

% y label
hyl = text(Left, Bottom + Height/2, ylabel, ...
		'HorizontalAlignment', 'right', ...
		'VerticalAlignment', 'middle');

if (HoldState==0)
	hold off
end

h = [hx hy hxl hyl];

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu