function StringArray = LoadStringArray(FileName, delimiters)
% StringArray = LoadStringArray(FileName)
%
% This function takes a text file and returns a
% 2d cellular array of strings.  Each line is
% split into tokens separated by one of the delimeters
%
% NB lines beginning with % are ignored

if (nargin == 1)
    delimiters = [9:13 32]; % White space characters
end

StringList = cell(1);
fp = fopen(FileName);
i=1;
while ~feof(fp)
	Line = fgetl(fp);
	if (Line(1) ~= '%')
		j=1;
		while ~isempty(Line)
			[tok Line] = strtok(Line, delimiters);
			StringArray{i,j} = tok;
			j = j+1;
		end
	i = i+1;
	end;
end;
fclose(fp);

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu