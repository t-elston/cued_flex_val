% xgvis(data, dist, edges, colors, glyphs, labels)
%
% initiates an xgvis process.
% data gives the starting positions of the points (N by P matrix)
% dist gives distances between points (N by N matrix)
% colors is an array of colors (N dimensional).
% The entries of colors should be numbers in the
% range 1 to 10, which are then translated to the
% color names used by xgvis.
%
% glyphs is an optional array of glyph numbers -
% to find out what number means what glyph look at
% the 'Glyph' menu on the 'Brush' screen of xgvis.
%
% Just about all of these arguments are optional.  If
% you don't want to specify one, give it as [].
%
% labels gives labels for the points.  This must be
% an array of numbers due to matlab fprintf sucks.
%
function xgvis(data, dist, edges, colors, glyphs, labels)

% generate file name
FileBase = tempname;

if(~isempty(data)) 
	nPoints = size(data,1);
	msave([FileBase, '.dat'], data);
end

if (nargin>=2 & ~isempty(dist))
	if (~isempty(data) & nPoints~=size(dist) | size(dist,1)~=size(dist,2))
		error('dist should be a NxN matrix for n data points');
	else
		nPoints = size(dist,1);
	end

%	msave([FileBase, '.dist'], dist);
	fp=fopen([FileBase, '.dist'], 'w');
	[m n] = size(dist);
	for row=1:m 
		for col=1:n
			if isfinite(dist(row, col))
				fprintf(fp, '%d ', dist(row, col));
			else
				fprintf(fp, '. ');
			end
		end
		fprintf(fp, '\n');
	end
	fclose(fp);
end

if (nargin>=3 & ~isempty(edges))
	[x y] = ind2sub(size(edges),find(edges));
	msave([FileBase, '.edges'], [x y]);
end

if (nargin>=4 & ~isempty(colors))
	% make color file
	
	if (length(colors) ~= nPoints)
		error('Number of colors should be same as number of points');
	end;
	
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

	colors = colors(:);
	ColorIndex = 1+colors-min(colors);
	%ColorIndex = colors(:);
	
	strmat = ColorNames(ColorIndex,:);
	
	fp = fopen([FileBase, '.colors'], 'w');
		
	fprintf(fp,'%s\n',strmat');
	
	fclose(fp);
end;

if (nargin >=5 & ~isempty(glyphs))

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

if (nargin >=5 & ~isempty(labels))
	fp = fopen([FileBase, '.row'], 'w');
	fprintf(fp, '%d\n', labels);
	fclose(fp);
end


% now run xgvis and clean up any old files...
eval(['! source ~/.cshrc ; ~/bin/xgvis ', FileBase, '&'])
%pause(0.1);
%eval(['! rm ', FileBase, '.*']);


% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu