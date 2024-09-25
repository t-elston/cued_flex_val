% msave - save an integer matrix in an ascii file
%
% msave(filename,matrix)
%

function msave(filename,matrix)

[m,n] = size(matrix);
formatstring = '%d';
for ii=2:n,
  formatstring = [formatstring,'\t%d'];
end
formatstring = [formatstring,'\n'];

outputfile = fopen(filename,'w');
fprintf(outputfile,formatstring,matrix');
fclose(outputfile);




% Written by Hajime Hirase
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu