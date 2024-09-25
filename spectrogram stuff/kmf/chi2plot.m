function h = chi2plot(x, V)
%CHI2PLOT Displays a chi squared probability plot.
%   H = CHI2PLOT(X, V) makes Chi Squared with V degrees of freedom
%   probability plot of the data in X. For matrix, X, CHI2PLOT displays
%   a plot for each column. %   H is a handle to the plotted lines.
%   
%   The purpose of a chi squared probability plot is to graphically assess
%   whether the data in X could come from a chi squared distribution. If the
%   data are chi squared the plot will be linear. Other distribution types
%   will introduce curvature in the plot.  

%   by Ken Harris, 2000.

[n, m] = size(x);
if n == 1
   x = x';
   n = m;
end

[sx i]= sort(x);
minx  = min(sx(1,:));
maxx  = max(sx(n,:));
range = maxx-minx;

if range>0
  minxaxis  = 0; % minx-0.025*range;
  maxxaxis  = maxx+0.025*range;
else
  minxaxis  = 0; % minx - 1;
  maxxaxis  = maxx + 1;
end

eprob = [0.5./n:1./n:(n - 0.5)./n];
y  = chi2inv(eprob,V)';

minyaxis  = minxaxis; chi2inv(0.25 ./n,V);
maxyaxis  = chi2inv((n-0.25) ./n,V);


p     = [0.01 0.10 0.25 0.5...
         0.75 0.90 0.95 0.98 0.99 0.997 0.999];

label1= str2mat('0.01 ','0.10','0.25','0.50');
label2= str2mat('0.75','0.90','0.95','0.98','0.99','0.997', '0.999');
label = [label1;label2];

tick  = chi2inv(p,V);

q1x = prctile(x,25);
q3x = prctile(x,75);
q1y = prctile(y,25);
q3y = prctile(y,75);
qx = [q1x; q3x];
qy = [q1y; q3y];


dx = q3x - q1x;
dy = q3y - q1y;
slope = dy./dx;
%centerx = (q1x + q3x)/2;
%centery = (q1y + q3y)/2;
maxx = max(x);
minx = min(x);
%maxy = centery + slope.*(maxx - centerx);
%miny = centery - slope.*(centerx - minx);

maxy = maxx;
miny = minx;

%mx = [minx; maxx];
%my = [miny; maxy];
mx = [0; maxx];
my = [0; maxy];

if nargout == 1
   h = plot(sx,y,'+',qx,qy,'-',mx,my,'-.');
%    h = plot(chi2cdf(sx),chi2cdf(y),'+',chi2cdf(qx),chi2cdf(qy),'-',chi2cdf(mx),chi2cdf(my),'-.');
else
   plot(sx,y,'+',qx,qy,'-',mx,my,'-.');
end

%set(gca,'YTick',tick,'YTickLabel',label);
%ylabel('Probability');
set(gca,'YLim',[minyaxis maxyaxis],'XLim',[minxaxis maxxaxis]);
ylabel('Predicted');
xlabel('\chi^2');
%title('Chi Squared Probability Plot');

grid on;
%whitebg(gcf,[0 0 0])



% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu