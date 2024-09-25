% FileExists(FileName)
%
% returns 1 if the file exists, 0 otherwise

function exists = FileExists(FileName)

fp = fopen(FileName, 'r');

if (fp==-1)
	exists = 0;
else
	fclose(fp);
	exists = 1;
end;

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu