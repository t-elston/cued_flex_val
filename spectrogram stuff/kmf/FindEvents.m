% [t chi2 p DoF] = FindEvents(x, WinLen, fS, fMax);
%
% Goes through a time series x and looks for outlying events as determined by
% their likelihood under the assumption that the time series is Gaussian
% in the range 0...fMax.  fMax defaults to Nyquist (fS/2)
%
% The output t gives the times of the start of periods of length WinLen, sorted
% starting with the most outlying.
%
% WinLen is the length of the window used for detecting outlying sequences. default 128
% fS is the sampling frequency. Defaults to 2.
%
% chi2 returns chi-squared statistic for each event - p returns its p-value
% DoF returns the number of degrees of freedom expected.

function [t, chi2, p, DoF] = FindEvents(x, WinLen, fS, fMax);

if (nargin<2 | isempty(WinLen)) WinLen = 128; end;
if (nargin<3 | isempty(fS)) fS = 2; end;
if (nargin<4 | isempty(fMax)) fMax = fS/2; end;

% deal with the fact that x might be multivariate
if size(x,1) == 1
	x = x(:);
end
nDims = size(x,2);

% calculate spectrograms
for d=1:nDims
	[Spec(:,:,d) f t] = specgram(x(:,d), WinLen, fS, hanning(WinLen), WinLen/2);
end
nTimeBins = length(t);
topF = max(find(f<fMax)); % frequency bin for top frequency


% calculate cross-spectral matrix
CsdInv = zeros(topF,nDims,nDims);
psgsc = zeros(topF,nTimeBins);
for fb=1:topF
	tmp = reshape(Spec(fb,:,:),[nTimeBins nDims]);
	CsdMH = cov(tmp)^-0.5;
	% calculate surprise-gram
	psgsc(fb,:) = sum(abs(tmp * CsdMH).^2,2)';
end

% MeanPower = mean(abs(Spec).^2,2);
% psgsc = abs(Spec).^2 ./ repmat(MeanPower, 1, size(Spec,2));

% calculate chi^2
if (topF < length(f))
	chi2 = psgsc(1,:) + 2*sum(psgsc(2:topF,:), 1); % bottom freq bin has only 1 dof
	DoF = nDims*(2*topF-1);
else
	chi2 = psgsc(1,:) + 2*sum(psgsc(2:topF-1,:), 1) +psgsc(topF,:); % same for top bin
	DoF = nDims*(2*topF-2);
end;

[sorted index] = sort(chi2);

% sort t and chi2 into descending order of chi2
order = flipud(index(:));
t = t(order);
chi2 = chi2(order);

if (nargout == 0 | nargout >=3)
	% calculate p-values
	p = 1-chi2cdf(chi2,DoF);
end

if nargout==0
	% display outlying epochs

	% make bookmark file
	%msave('bookmk', (order-1)*WinLen/2); % save bookmarks for .eeg file
	
	% calculate range for y axis
	Minx = prctile(x(:),1);
	Maxx = prctile(x(:),99);

	for i = 1:nTimeBins
	   subplot(3,1,1)
	   xr = (order(i)-1)*WinLen/2 + (1:WinLen); % time range
	   plot(xr/fS, x(xr,:));
	   xlim([xr(1) xr(end)]/fS);	
%	   ylim([Minx Maxx]);
%	   subplot(3,1,2)
%	   plot(f,20*log10([abs(Spec(:,order(i))), sqrt(MeanPower)]));
%	   title('power spectrum');
	   subplot(3,1,3);	
	   plot(f(1:topF),psgsc(:,order(i)));
	   title('Power spectrum scaled by mean');
	   fprintf('bin %d time %f: chi^2 %f p=%g q=%g\n', order(i), t(i), chi2(i),p(i),i/nTimeBins);
	   pause
	end
%end


%return

%if 0
	% coefficient of variation
	MeanAmp = mean(abs(Spec),2);
	figure; %(1)
	subplot(3,1,1);
	plot(f, 20*log10([MeanAmp, sqrt(MeanPower)]));
	legend('Mean Amp', 'RMS Amp');
	ylabel('dB');
	grid on;
	subplot(3,1,2);
	plot(f, MeanPower./MeanAmp.^2);
	title('coefficient of variation');
	title('outlyingness');
	subplot(3,1,3)
	hist(chi2, 100);
end



% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu