function [RLspecxy]=LinearizeSpectra(LRRunsSpectra,numColumns,RightOrLeft);
[r c]=size(LRRunsSpectra);
if numColumns==12
    xCoord=13;
    yCoord=14;
elseif numColumns==2
    xCoord=3;
    yCoord=4;
elseif numColumns==3 %this is for the CFCs? has been modified to work with CFC outputs, but unsure if this affects anything else.
    xCoord=4;
    yCoord=5;
end
for a=1:r
    %%%This bit to get linearized trajectory for each run%%
    runX=LRRunsSpectra{a,xCoord};
    runY=LRRunsSpectra{a,yCoord};
    %This checks if there are artifacts to do with the 'end' event trigger
    if ~isscalar(runX)
        %Get the reference point and transformed circular trajectory
        ref=[(min(runX)+max(runX))/2 (min(runY)+max(runY))/2];
        theta=atan2(runY-ref(:,2),runX-ref(:,1))+pi;

        RefLin=linspace(0,2*pi,(min(abs(runX))+max(abs(runX)))*2+(min(abs(runY))+max(abs(runY)))*2)';
        RefLin(:,2)=1:1:length(RefLin);
        RefLin(:,2)=circshift(RefLin(:,2),-235); %%%%%%%%%%%%%%%%%%%%%%%%%%These two lines, and the ones below, are added to correct the orientation
        RefLin(:,2)=flipud(RefLin(:,2));
        %If it's left runs, do this extra step to correct mirror imaging of
        %the trajectories
        if RightOrLeft==2
            LinRef=linspace(0,2*pi,(min(abs(runX))+max(abs(runX)))*2+(min(abs(runY))+max(abs(runY)))*2)'; 
            LinRef=circshift(LinRef,round(length(RefLin)/2));
            RefLin(:,2)=1:1:length(LinRef); %Note this is an awkward way of calculating linearized coordinates for left trajectories
            RefLin(:,2)=circshift(RefLin(:,2),-235); %%%%%%%%%%%%%%%%%%%%%%%%%%These two lines, and the ones below, are added to correct the orientation
            RefLin(:,2)=flipud(RefLin(:,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%check if the above actually does what it is
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%supposed to do!!!!%%%%%%%%%%%%%%%%%%%%%%%%%
        else
        end
        %Create new linearized trajectory
        newS=[];
        for i=1:length(theta)
            dist=(RefLin(:,1)-theta(i)).^2;
            [m,ind]=min(dist);
            newS(i,:)=RefLin(ind(1),2);
        end

        LinS=[];
        %For each spectrum (i.e. HPC/PFCxtheta/LG/MG/HG), get sampled
        %points in space (linearized trajecotry) and an estimate spectral
        %value. This bit also averages the points in space that have been
        %duplicated.
        for c=1:numColumns;
            for i=1:length(newS)
                LinS(i)=mean(LRRunsSpectra{a,c}(find(newS==newS(i))));
            end
        %Construct the vector pairs of linearized trajectory x spectral
        %estimate
        zz=[newS LinS'];
        [C,ia,ic] = unique(zz,'rows');

    %     t=(1:length(C));
    %     sigma=3;
    %     dt=diff(t);
    %     dt=dt(1);
    %     a=1/(sqrt(2*pi)*sigma);
    %     gfilter = dt*a*exp(-0.5*((t - mean(t)).^2)/(sigma^2));
    %     i = gfilter < dt*a*1.e-6;
    %     gfilter(i) = [];
    % 
    %     filtC=conv(C(:,2),gfilter,'same');
        %Collate the space x spectrum pairs
        RLspecxy{a,c}=C;
        LinS=[];
        end
    else
    end
end

    %     figure;
    %     subplot(2,1,1);
    %     bar(C(:,1),C(:,2));
    %     subplot(2,1,2);
    %     bar(C(:,1),filtC);
    %     figure;
    %     imagesc(LRRunsSpectra{a,18});
    %     axis xy;
