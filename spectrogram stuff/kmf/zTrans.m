% helper function zTrans computes Fisher's z-transform
function z = zTrans(x)
z = 0.5 * log((1+x) ./ (1-x));


% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu