function bsave(fname, matrix, arg2)
% saves a binary file
% bsave(fname, matrix, arg2)
% uses command fwrite(fp, matrix, arg2);

fp = fopen(fname, 'w');

fwrite(fp, matrix, arg2);

fclose(fp);

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu