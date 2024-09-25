% [Data Headers] = LoadHeadedArray(FileName)
%
% Loads up a white space-separated numeric array which has a 
% single line on top giving column headings.

function [Data, Headers] = LoadHeadedArray(FileName)

fp = fopen(FileName);

HeadLine = fgetl(fp);

Headers = {};
cnt = 1;
while(HeadLine)
	[Token HeadLine] = strtok(HeadLine);
	Headers{cnt} = Token;
	cnt = cnt+1;
end;

nColumns = length(Headers);

Data = fscanf(fp, '%f', [nColumns inf])';

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu