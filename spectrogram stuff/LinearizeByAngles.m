function [RRuns LRuns RNormLinRun LNormLinRun]=LinearizeByAngles(PX,PY,P_ts,AllRuns);

%Get trajectories for left and right runs
k=1;
l=1;
for i=1:length(AllRuns)
        if max(AllRuns{1,i}(:,1))>175 %%%This should be changed to 175/50
            RRuns{1,k}=AllRuns{1,i};
            k=k+1;
        elseif min(AllRuns{1,i}(:,1))<-175  %%%This should be changed to -175/50
            LRuns{1,l}=AllRuns{1,i};
            l=l+1;
        else
            
        end        
end



% %Create a cell array containing all the angles from the reference point to
% %the coodinate in each run. This bit does it for each individual runs,
% %which is scrapped in favour of a "true" reference point based on all runs
% %on the same side for the given session (below).
% 
% %Right
% for i=1:length(RRuns)
%     [r c]=size(RRuns{1,i});
%     for j=1:r
%         ref=[(min(RRuns{1,i}(:,1))+max(RRuns{1,i}(:,1)))/2 (min(RRuns{1,i}(:,2))+max(RRuns{1,i}(:,2)))/2];
%         Rtheta{1,i}(j,1)=atan2(RRuns{1,i}(j,2)-ref(:,2),RRuns{1,i}(j,1)-ref(:,1))+pi;
%     end
% end
% 
% %Left
% for i=1:length(LRuns)
%     [r c]=size(LRuns{1,i});
%     for j=1:r
%         ref=[(min(LRuns{1,i}(:,1))+max(LRuns{1,i}(:,1)))/2 (min(LRuns{1,i}(:,2))+max(LRuns{1,i}(:,2)))/2];
%         Ltheta{1,i}(j,1)=atan2(LRuns{1,i}(j,2)-ref(:,2),LRuns{1,i}(j,1)-ref(:,1))*-1+pi;
%     end
% end

%This bit calculates the reference point for the angle calculation along
%the track based on all runs on a given side 

overlaidR=cell2mat(RRuns');
overlaidL=cell2mat(LRuns');

% %Right
% 
%         ref=[(min(overlaidR(:,1))+max(overlaidR(:,1)))/2 (min(overlaidR(:,2))+max(overlaidR(:,2)))/2];
%         Rtheta{1,i}(j,1)=atan2(RRuns{1,i}(j,2)-ref(:,2),RRuns{1,i}(j,1)-ref(:,1))+pi;
% 
% 
% %Left
% 
%         ref=[(min(overlaidL(:,1))+max(overlaidL(:,1)))/2 (min(overlaidL(:,2))+max(overlaidL(:,2)))/2];
%         Ltheta{1,i}(j,1)=atan2(LRuns{1,i}(j,2)-ref(:,2),LRuns{1,i}(j,1)-ref(:,1))*-1+pi;


%Right
for i=1:length(RRuns)
    [r c]=size(RRuns{1,i});
    for j=1:r
        ref=[(min(overlaidR(:,1))+max(overlaidR(:,1)))/2 (min(overlaidR(:,2))+max(overlaidR(:,2)))/2];
        Rtheta{1,i}(j,1)=atan2(RRuns{1,i}(j,2)-ref(:,2),RRuns{1,i}(j,1)-ref(:,1))+pi;
    end
end

%Left
for i=1:length(LRuns)
    [r c]=size(LRuns{1,i});
    for j=1:r
        ref=[(min(overlaidL(:,1))+max(overlaidL(:,1)))/2 (min(overlaidL(:,2))+max(overlaidL(:,2)))/2];
        Ltheta{1,i}(j,1)=atan2(LRuns{1,i}(j,2)-ref(:,2),LRuns{1,i}(j,1)-ref(:,1))*-1+pi;
    end
end

%Linearizing trajectories by normalizing angles into an arbitary position
%system
RefLin=linspace(0,2*pi,(min(abs(PX))+max(abs(PX)))*2+(min(abs(PY))+max(abs(PY)))*2)';
RefLin(:,2)=1:1:length(RefLin);
RefLin(:,2)=circshift(RefLin(:,2),-235);
RefLin(:,2)=flipud(RefLin(:,2));

LinRef=linspace(0,2*pi,(min(abs(PX))+max(abs(PX)))*2+(min(abs(PY))+max(abs(PY)))*2)'; 
LinRef=circshift(LinRef,round(length(RefLin)/2));
LinRef(:,2)=1:1:length(LinRef);
LinRef(:,2)=circshift(LinRef(:,2),-235);
LinRef(:,2)=flipud(LinRef(:,2));
%LinRef(2,:)=length(LinRef):-1:1;


%Fetch the closest normalized arbitary value for a given angle/coordinate

%Right
for i=1:length(Rtheta)
    for j=1:length(Rtheta{1,i})
            dist=(RefLin(:,1)-Rtheta{1,i}(j,:)).^2;
            [m,ind]=min(dist);
            RNormLinRun{1,i}(j,:)=RefLin(ind(1),2);
    end
end

%Left
for i=1:length(Ltheta)
    for j=1:length(Ltheta{1,i})
            dist=(LinRef(:,1)-Ltheta{1,i}(j,:)).^2;
            [m,ind]=min(dist);
            LNormLinRun{1,i}(j,:)=LinRef(ind(1),2);
    end
end

%This bit removes the single point "trajectories" leftover by introducing
%the "end" trigger that is not associated with the actual tracking.
RNormLinRun(cellfun(@isscalar, RNormLinRun)) = [];
LNormLinRun(cellfun(@isscalar, LNormLinRun)) = [];

%Fetch spike time-stamps corresponding to the xy coordinates and the
%linearized vector

% UnitID=UnitID;
% 
% %Right
% for i=1:length(RRuns)
%     Unitind=find(UnitID>=min(RRuns{1,i}(:,3)) & UnitID<=max(RRuns{1,i}(:,3)));
%     ts=UnitID(Unitind);
%     M=length(ts);
%     for ii=1:M
%         tdiff=(RRuns{1,i}(:,3)-ts(ii)).^2;
%         [m,ind]=min(tdiff);
%         Rspkxy{1,i}(ii,1)=RRuns{1,i}(ind(1),1);
%         Rspkxy{1,i}(ii,2)=RRuns{1,i}(ind(1),2);
%         Rspkxy{1,i}(ii,3)=RNormLinRun{1,i}(ind(1));
%     end
% end
% 
% %Left
% for i=1:length(LRuns)
%     Unitind=find(UnitID>=min(LRuns{1,i}(:,3)) & UnitID<=max(LRuns{1,i}(:,3)));
%     ts=UnitID(Unitind);
%     M=length(ts);
%     for ii=1:M
%         tdiff=(LRuns{1,i}(:,3)-ts(ii)).^2;
%         [m,ind]=min(tdiff);
%         Lspkxy{1,i}(ii,1)=LRuns{1,i}(ind(1),1);
%         Lspkxy{1,i}(ii,2)=LRuns{1,i}(ind(1),2);
%         Lspkxy{1,i}(ii,3)=LNormLinRun{1,i}(ind(1));
%     end
% end
% 
% % figure;
% % subplot(2,1,1);
% %     hold all;
% %     for i=1:length(Rspkxy)
% %         hist(Rspkxy{1,i}(:,3),1:10:1252);     
% %     end
% % subplot(2,1,2);
% %     hold all;
% %     for i=1:length(Lspkxy)
% %         hist(Lspkxy{1,i}(:,3),1:10:1252);     
% %     end
% 
% figure;
% subplot(2,1,1);
%         hist(cell2mat(Rspkxy'),1:10:1252);     
% subplot(2,1,2);
%         hist(cell2mat(Lspkxy'),1:10:1252);     

    
%% %Sanity check that this is actually working (i.e. example code):
% load 20140206_C1.mat
% SpkPath;
% k=1;
% l=1;
% for i=1:length(AllRuns)
%         if max(AllRuns{1,i}(:,1))>200
%             RRuns{1,k}=AllRuns{1,i};
%             k=k+1;
%         elseif min(AllRuns{1,i}(:,1))<-200
%             LRuns{1,l}=AllRuns{1,i};
%             l=l+1;
%         end        
% end
% i=1;
% j=24;
% ref=[(min(RRuns{1,i}(:,1))+max(RRuns{1,i}(:,1)))/2 (min(RRuns{1,i}(:,2))+max(RRuns{1,i}(:,2)))/2];
% figure;
% scatter(getcolumn(RRuns{1,1}(:,1:2),1),getcolumn(RRuns{1,1}(:,1:2),2),'DisplayName','RRuns{1,1}(:,1:2)(:,2) vs. RRuns{1,1}(:,1:2)(:,1)','YDataSource','RRuns{1,1}(:,1:2)(:,2)');figure(gcf)
% hold on
% scatter(ref(:,1),ref(:,2));
% Rtheta{1,i}=atan2(RRuns{1,i}(j,2)-ref(:,2),RRuns{1,i}(j,1)-ref(:,1));
% scatter(RRuns{1,i}(j,1),RRuns{1,i}(j,2));
% theta = Rtheta{1,1};
% r = 300; % magnitude (length) of arrow to plot
% x = ref(:,1); y = ref(:,2);
% u = r * cos(theta); % convert polar (theta,r) to cartesian
% v = r * sin(theta);
% quiver(x,y,u,v);

% % For the left side
% i=1;
% j=24;
% ref=[(min(LRuns{1,i}(:,1))+max(LRuns{1,i}(:,1)))/2 (min(LRuns{1,i}(:,2))+max(LRuns{1,i}(:,2)))/2];
% figure;
% scatter(getcolumn(LRuns{1,1}(:,1:2),1),getcolumn(LRuns{1,1}(:,1:2),2),'DisplayName','LRuns{1,1}(:,1:2)(:,2) vs. LRuns{1,1}(:,1:2)(:,1)','YDataSource','RRuns{1,1}(:,1:2)(:,2)');figure(gcf)
% hold on
% scatter(ref(:,1),ref(:,2));
% Ltheta{1,i}=atan2(LRuns{1,i}(j,2)-ref(:,2),LRuns{1,i}(j,1)-ref(:,1));
% scatter(LRuns{1,i}(j,1),LRuns{1,i}(j,2));
% theta = Ltheta{1,1};
% r = 300; % magnitude (length) of arrow to plot
% x = ref(:,1); y = ref(:,2);
% u = r * cos(theta); % convert polar (theta,r) to cartesian
% v = r * sin(theta);
% quiver(x,y,u,v);
