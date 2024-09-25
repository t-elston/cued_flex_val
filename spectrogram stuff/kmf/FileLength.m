% FileLength(FileName)
%
% Returns the length of the file in bytes
% If the file is not there, returns -1.

function n = FileLength(FileName);

fp = fopen(FileName, 'r');

if (fp == -1)
	n = -1;
	return;
end;

fseek(fp, 0, 'eof');
n = ftell(fp);

fclose(fp);

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu