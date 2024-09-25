% [nLines, nWords, nChars] = WordCount(FileName)
%
% Counts number of lines words and characters in 
% specified file.  Works by running unix wc command
%
% NB if you have only one output variable it is number
% of LINES, not words.

function [nLines, nWords, nChars] = WordCount(FileName)

Command = ['wc ', FileName];

[status, result] = unix(Command);

if (status~=0) 
	error(['WordCount failed!  Unix error: ', result]);
end;

Parsed =  sscanf(result, '%d');

nLines = Parsed(1);
nWords = Parsed(2);
nChars = Parsed(3);

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu