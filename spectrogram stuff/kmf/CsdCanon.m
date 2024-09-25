% CsdCanon(Csd, Grp, f, nCanons, AnchorChannel, Gamma)
%
% Does a display based on canonical coherence analysis of a Csd.
% The input Csd should be a 3d array, with the last
% two dimensions being channel number
%
% Grp is a vector of 1's and 2's specifying which channels belong
% to which group.
%
% f is an array of frequencies, as produced by mtcsd.
% default is 0...1
%
% nCanons is the number of canonical coherences to consider
%
% AnchorChannel is the channel which phases are measured
% relative to.  It's an array with one entry per canonical var.
% if you only give a single number, they will all be the same
%
% gamma is a contrast modifier for the Complex Color plots
% (see ComplexImage).  Default 0.8

function CsdEig(Csd, Grp, f, nCanons, AnchorChannel, Gamma)

if (nargin<3 | isempty(f))	f = ((1:size(Csd,1)) - 1) / (size(Csd,1)); end
if (nargin<4 | isempty(nCanons)) nEigs = 1; end;
if (nargin<5 | isempty(AnchorChannel)) AnchorChannel = ones(nCanons,1); end
if (nargin<6 | isempty(Gamma)) Gamma = 0.8; end;

nCh = size(Csd,3);
nFreqs = length(f);
nSubPlots = 2*nCanons+1;

if any(Grp ~=1 & Grp ~=2)
	error('Grp should consist of 1s and 2s');
end

if (length(AnchorChannel) == 1)
	AnchorChannel = repmat(AnchorChannel, nCanons, 1);
elseif (length(AnchorChannel) ~= nCanons)
	error('AnchorChannel should be nCanons elements long');
end

Grp1 = find(Grp==1);
nGrp1 = length(Grp1);
Grp2 = find(Grp==2);
nGrp2 = length(Grp2);

if (nCanons > min(nGrp1,nGrp2))
	error('Cannot have more canonical covariates than dimensions of the smaller group');
end;


CohereStore = zeros(nFreqs,nCanons);
VecStore1 = zeros(nFreqs, nGrp1, nCanons);
VecStore2 = zeros(nFreqs, nGrp2, nCanons);

% find out which group the anchor channel is in, and which index in that group
for i=1:nCanons
	AnchorGrp(i) = Grp(AnchorChannel(i));
	AnchorIndex(i) = sum(Grp==AnchorGrp(i) & 1:nCh <= AnchorChannel(i));
end;

for f0=1:length(f)
	CsMat = squeeze(Csd(f0,:,:));
	
	CXX = CsMat(Grp1,Grp1);
	CXY = CsMat(Grp1,Grp2);
	CYX = CsMat(Grp2,Grp1);
	CYY = CsMat(Grp2,Grp2);
	
	CXXMH = CXX ^ -0.5;
	CYYMH = CYY ^ -0.5;
	
	% matrix to do svd on ...
	M = CXXMH * CXY * CYYMH;
	% do svd
	[c s d] = svd(M);
	% then calculate a and b (the projections)
	a = CXXMH * c;
	b = CYYMH * d;
	
	%R2 is diagonal elements of s --- but sometimes s is already a vector which will mess up the diag function
	if (min(size(s)) == 1) 
		R2 = s(:);
	else
		R2 = diag(s);
	end
	
	CohereStore(f0,:) = R2(1:nCanons)'; % store cannonical coefficient
	
	SD1 = sqrt(diag(CXX));
	SD2 = sqrt(diag(CYY));
	
	% store correlations of canonical variables with original variables
	VecStore1(f0,:,:) = NormVector(CXX*a(:,1:nCanons)) ;
	VecStore2(f0,:,:) = NormVector(CYY*b(:,1:nCanons)) ;
			
%	fprintf('%f %f %f %f\n', CYY, b, SD2, CYY*b/SD2);
%	pause
	
end

% Now make the phase be relative to the anchor channel...

for i=1:nCanons
	if AnchorGrp(i) == 1
		Anchor = VecStore1(:,AnchorIndex(i), i);
	else
		Anchor = VecStore2(:,AnchorIndex(i), i);
	end
	
	AnchorPhase = Anchor ./ abs(Anchor);
	
	VecStore1(:,:,i) = VecStore1(:,:,i) ./ repmat(AnchorPhase,[1 nGrp1]);
	VecStore2(:,:,i) = VecStore2(:,:,i) ./ repmat(AnchorPhase,[1 nGrp2]);
end
%Plot canonical coherences 
subplot(nSubPlots,1,1)
cla
plot(f, real(CohereStore));
axis([0 max(f) 0 1]);
	
% Plot Canonical vectors
for e=1:nCanons
	subplot(nSubPlots,1,2*e);
	ComplexImage(f,1:nGrp1,VecStore1(:,:,e)',Gamma) 
	
	subplot(nSubPlots,1,1+2*e);
	ComplexImage(f,1:nGrp2,VecStore2(:,:,e)',Gamma) 

end	


% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu