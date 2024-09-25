% [Reg Residuals] = MultiReg(Dependent, Independent, Axes)
% Do a multiple linear regression fit
%
% Dependent is a matrix.  each column corresponds to an individual linear
% regression so if you are using it to regress entire waveshapes, sample
% number is the second coordinate.  Rows correspond to individual data
% points. so spike number is the first coordinate.
%
% Independent is also a matrix.  Each column corresponds to an independent
% variable.  You don't need to include a column of 1's - the program will do
% that for you.  Rows correspond to individual data points. so spike number
% is again the first coordinate.
%
% Set Axes to override default axes.
%
% Reg (output variable) will contain coefficients of 1 and of the independent
% variable, with sample number along the 2nd coordinate.
%
% If no output argument is requested, you will get a lot of fits to individual cases.
%
% Residuals (output variable) is the same shape as Dependent, containing
% resiuals after the fit.

function [Reg, Residuals] = MultiReg(Dependent, Independent, Axes);

% Check sizes
nInstances = size(Dependent, 1);
if nInstances ~= size(Independent, 1)
	error('Dependent and Independent must have the same number of rows');
end;

nSamples =size(Dependent, 2);
nIndependentVars = size(Independent, 2);

% set Axes if not set
if nargin <3 Axes = [0,121,-800,200]; end;

% zero means of Independent variables columnwise
Independent = Independent - ones(nInstances, 1) * mean(Independent, 1);

% find top and bottom values for independent variables
TopInd = max(Independent, [], 1);
BottomInd = min(Independent, [], 1);

% add extra column of 1's
Independent = [ones(nInstances, 1) , Independent];

% Do the least squares fit
Reg = Independent \ Dependent;

% Reconstructed waveforms
Reconstructed = Independent * Reg;

% Residuals
Residuals =Dependent - Reconstructed;

% Residual sum sq error
SumSqErrors = sum(Residuals.*Residuals, 2);
TotalError = sum(SumSqErrors)

% plot templates

plot([Reg(1,:)', Reg(1, :)' + TopInd(1)*Reg(2,:)', Reg(1, :)' + BottomInd(1)*Reg(2,:)']);

% only continue if there was _no_ outupt argument requested
if (nargout >=1) return; end;

pause
for i=1:nIndependentVars+1
	plot(Reg(i,:));
	title(sprintf('Template %d', i));
	pause;
end;


for i=1:nSamples
	plot([Dependent(i,:)', Reconstructed(i,:)', Reg(1,:)', Residuals(i,:)', zeros(nSamples, 1)]);
	axis(Axes);
	title(sprintf('Spike %d: Sum squared residuals %f', i, SumSqErrors(i)));
	pause;
end


% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu