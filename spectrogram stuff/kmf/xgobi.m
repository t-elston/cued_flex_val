% [colors glyphs] = xgobi(data, colors, glyphs, labels)
%
% initiates an xgobi process from the array data
% colors is an optional array of colors.
% The entries of colors should be numbers in the
% range 1 to 10, which are then translated to the
% color names used by xgobi.
%
% glyphs is an optional array of glyph numbers -
% to find out what number means what glyph look at
% the 'Glyph' menu on the 'Brush' screen of xgobi.
% "Points" is 31.
%
% Any Infs or NaNs will cause removal of that data line
%
% labels gives labels for the points.  This must be
% an array of numbers due to matlab fprintf sucks.
%
% optional output arguments give final brushing - also not yet
% implemented.

function xgobi(data, colors, glyphs, labels)

% generate file name
FileBase = tempname;

% Look for missing data
Infinities = ~isfinite(data);
if (any(Infinities(:)))
	% find rows with infinities 
	InfinityRows = find(any(Infinities,2));
	% delete those rows
	data(InfinityRows, :) = [];	
	if exist('colors'), colors(InfinityRows) = []; end;
	if exist('glyphs'), glyphs(InifinityRows) = []; end;
	disp([num2str(length(InfinityRows)), ' rows containing infinities were deleted.']);
end;

nPoints = size(data,1);



if nargin>=2 & ~isempty(colors);
	% make color file
	
	if (length(colors)==1)
		colors = repmat(colors, nPoints, 1);
	elseif (length(colors) ~= nPoints)
		error('Number of colors should be same as number of points');
	end;
	
	if 1 % new colors
	ColorNames = [
	'Default           ';... % color 1
	'DeepPink          ';... % color 2
	'OrangeRed1        ';... % color 3
	'DarkOrange        ';... % color 4
	'Gold              ';... % color 5
	'Yellow            ';... % color 6
	'DeepSkyBlue1      ';... % color 7
	'SlateBlue1        ';... % color 8
	'YellowGreen       ';... % color 9
	'MediumSpringGreen ';... % color 10
	'MediumOrchid      ';... % color 11
	];
	else % old colors
	ColorNames = [
	'Default     ';...
	'Red         ';...
	'Green       ';...
	'SkyBlue     ';...
	'Orange      ';...
	'Yellow      ';...
	'YellowGreen ';...
	'HotPink     ';...
	'Orchid      ';...
	'Peru        ';...
	'SlateBlue   '...
	];
	end
	
	colors = colors(:);
	ColorIndex = 1+colors-min(colors);

	% print out color map
	for i=min(colors):max(colors)
		fprintf('Color %d -> %s\n', i, ColorNames(1+i-min(colors),:));
	end
	
	strmat = ColorNames(ColorIndex,:);
	
	fp = fopen([FileBase, '.colors'], 'w');
		
	fprintf(fp,'%s\n',strmat');
	
	fclose(fp);
	
if 0
	%make .resources file
	fp = fopen([FileBase, '.resources'], 'w');
	for i=2:size(ColorNames,1)
		fprintf(fp, '*brushColor%d:%s\n', i, ColorNames(i,:));
	end;
	fclose(fp);
end
	
end;

if nargin >=3 & ~isempty(glyphs)
	% make glyphs file
	if length(glyphs) == 1
		glyphs = repmat(glyphs, nPoints, 1);
	end;
	
	if length(glyphs) ~= nPoints
		error('number of glyphs should be same as number of points or 1');
	end;
	fp = fopen([FileBase, '.glyphs'], 'w');
	fprintf(fp, '%d\n', glyphs);
	fclose(fp);
end;

if (nargin >=4 & ~isempty(labels))
	fp = fopen([FileBase, '.row'], 'w');
	fprintf(fp, '%d\n', labels);
	fclose(fp);
end

msave([FileBase, '.dat'], data);

% now run xgobi and clean up any old files...
eval(['! ~/bin/xgobi ', FileBase, '&'])
%pause(0.1);
%eval(['! rm ', FileBase, '.*']);


% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu