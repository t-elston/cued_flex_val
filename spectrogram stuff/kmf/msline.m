function h = msline(Bins)
% h=msline(Bins) Add median smooth line to scatter plot.
%   MSLINE superimposes the running median line on each line object
%   in the current axis (Except LineStyles '-','--','.-'.)
% 
%   H = MSLINE returns the handle to the line object(s) in H.
%   
%   See also MedianSmooth
%
%  Bins is an optional argument passed to MedianSmooth.  It may
%  be a number of bins, or the specification of the left edges.
%  Defautl = 20.

if (nargin<1)
	Bins = 20;
end;

holdstate = ishold;
hold on;
lh = findobj(get(gca,'Children'),'Type','line');
if nargout == 1, 
   h = [];
end
count = 0;
for k = 1:length(lh)
    xdat = get(lh(k),'Xdata');
    ydat = get(lh(k),'Ydata');
    datacolor = get(lh(k),'Color');
    style = get(lh(k),'LineStyle');
    if ~strcmp(style,'-') & ~strcmp(style,'--') & ~strcmp(style,'-.')
       count = count + 1;
       
       [m BinsO] = MedianSmooth(xdat, ydat, Bins);
       
       newline = plot(BinsO, m, 'Color', datacolor);
       
       if nargout == 1
           h(count) = newline;    
       end
   end
end

if count == 0
   disp('No allowed line types found. Nothing done.');
end

if (holdstate==0)
	hold off
end;

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu