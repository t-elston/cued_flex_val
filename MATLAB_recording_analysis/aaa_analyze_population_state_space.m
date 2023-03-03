
% load the data
load('C:\Users\Thomas Elston\Documents\MATLAB\Projects\CuedFlexVal\population_state_space_data.mat');

HPC_cue = nanmean(HPC_cue_PCs, 3);
OFC_cue = nanmean(OFC_cue_PCs, 3);

HPC_choice = nanmean(HPC_choice_PCs, 3);
OFC_choice = nanmean(OFC_choice_PCs, 3);

n_boots = size(HPC_cue_PCs,3);

hpc_cue_dist=[];
ofc_cue_dist=[];

cue_PCs = [1,2];

% let's get some distance stats for the cue period
for b = 1:n_boots
    
    hpc_cue_dist(b,1) = sqrt(sum((HPC_cue_PCs(1,cue_PCs,b) - HPC_cue_PCs(2,cue_PCs,b)).^2, 2)) +...
        sqrt(sum((HPC_cue_PCs(3,cue_PCs,b) - HPC_cue_PCs(4,cue_PCs,b)).^2, 2));
    
    hpc_cue_dist(b,2) = sqrt(sum((HPC_cue_PCs(1,cue_PCs,b) - HPC_cue_PCs(4,cue_PCs,b)).^2, 2)) +...
        sqrt(sum((HPC_cue_PCs(2,cue_PCs,b) - HPC_cue_PCs(3,cue_PCs,b)).^2, 2));
    
    ofc_cue_dist(b,1) = sqrt(sum((OFC_cue_PCs(1,cue_PCs,b) - OFC_cue_PCs(2,cue_PCs,b)).^2, 2)) +...
        sqrt(sum((OFC_cue_PCs(3,cue_PCs,b) - OFC_cue_PCs(4,cue_PCs,b)).^2, 2));
    
    ofc_cue_dist(b,2) = sqrt(sum((OFC_cue_PCs(1,cue_PCs,b) - OFC_cue_PCs(4,cue_PCs,b)).^2, 2)) +...
        sqrt(sum((OFC_cue_PCs(2,cue_PCs,b) - OFC_cue_PCs(3,cue_PCs,b)).^2, 2));
    
end % of looping over bootstraps

cmap = cbrewer('qual','Paired',12);

% plotting for CUE PHASE
cue_figure = figure; 
set(cue_figure,'Position',[400 400 800 600])
%define axes
HPC_PC_ax  = axes('Position',[.25   .6  .3  .3]);
OFC_PC_ax  = axes('Position',[.65   .6  .3  .3]);
Hist_ax    = axes('Position',[.25   .2  .7  .3]);

axes(HPC_PC_ax);
hold on
plot3(HPC_cue(1,1), HPC_cue(1,2), HPC_cue(1,3), 'Marker','s','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(6,:));
plot3(HPC_cue(2,1), HPC_cue(2,2), HPC_cue(2,3), 'Marker','o','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(6,:));
plot3(HPC_cue(3,1), HPC_cue(3,2), HPC_cue(3,3), 'Marker','s','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(2,:));
plot3(HPC_cue(4,1), HPC_cue(4,2), HPC_cue(4,3), 'Marker','o','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(2,:));
xlim([-15, 15]);
ylim([-15, 15]);
xlabel('PC1');
ylabel('PC2');
legend({'Context A, Cue 1','Context A, Cue 2', 'Context B, Cue 1','Context B, Cue 2'}, 'Location', [.08 .8 .1 .1]);
title('HPC');
HPC_cue_pval = 1- sum((hpc_cue_dist(:,2) - hpc_cue_dist(:,1)) > 0) / n_boots;

axes(OFC_PC_ax);
hold on
plot3(OFC_cue(1,1), OFC_cue(1,2), OFC_cue(1,3), 'Marker','s','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(5,:));
plot3(OFC_cue(2,1), OFC_cue(2,2), OFC_cue(2,3), 'Marker','o','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(5,:));
plot3(OFC_cue(3,1), OFC_cue(3,2), OFC_cue(3,3), 'Marker','s','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(1,:));
plot3(OFC_cue(4,1), OFC_cue(4,2), OFC_cue(4,3), 'Marker','o','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(1,:));
xlim([-15, 15]);
ylim([-15, 15]);
title('OFC');
OFC_cue_pval = 1 - sum((ofc_cue_dist(:,2) - ofc_cue_dist(:,1)) > 0) / n_boots;

axes(Hist_ax);
hold on
histogram(hpc_cue_dist(:,1), 'Normalization','Probability','BinWidth',.1,'FaceColor',cmap(4,:),'EdgeColor','none');
histogram(hpc_cue_dist(:,2), 'Normalization','Probability','BinWidth',.1,'FaceColor',[.2 .2 .2],'EdgeColor','none');
histogram(ofc_cue_dist(:,1), 'Normalization','Probability','BinWidth',.1,'FaceColor',cmap(3,:),'EdgeColor','none');
histogram(ofc_cue_dist(:,2), 'Normalization','Probability','BinWidth',.1,'FaceColor',[.6 .6 .6],'EdgeColor','none');
legend({'HPC Within-Context','HPC Across-Context', 'OFC Within-Context', 'OFC Across-Context'},'Location', [.08 .4 .1 .1]);
h = gca;
h.YAxis.Visible = 'off';
xlabel('Euclidean Distance (a.u.)');



figure; 
hold on
plot3(HPC_choice(1,1),HPC_choice(1,2), HPC_choice(1,4), 'Marker','o','color','r', 'MarkerSize',15, 'MarkerFaceColor','r');
plot3(HPC_choice(2,1),HPC_choice(2,2), HPC_choice(2,4), 'Marker','^','color','r', 'MarkerSize',15, 'MarkerFaceColor','r');
plot3(HPC_choice(3,1),HPC_choice(3,2), HPC_choice(3,4), 'Marker','s','color','r', 'MarkerSize',15, 'MarkerFaceColor','r');
plot3(HPC_choice(4,1),HPC_choice(4,2), HPC_choice(4,4), 'Marker','pentagram','color','r', 'MarkerSize',15, 'MarkerFaceColor','r');
plot3(HPC_choice(5,1),HPC_choice(5,2), HPC_choice(5,4), 'Marker','o','color','b', 'MarkerSize',15, 'MarkerFaceColor','b');
plot3(HPC_choice(6,1),HPC_choice(6,2), HPC_choice(6,4), 'Marker','^','color','b', 'MarkerSize',15, 'MarkerFaceColor','b');
plot3(HPC_choice(7,1),HPC_choice(7,2), HPC_choice(7,4), 'Marker','s','color','b', 'MarkerSize',15, 'MarkerFaceColor','b');
plot3(HPC_choice(8,1),HPC_choice(8,2), HPC_choice(8,4), 'Marker','pentagram','color','b', 'MarkerSize',15, 'MarkerFaceColor','b');
legend({'C1, V1','C1, V2','C1, V3', 'C1, V4','C2, V1','C2, V2','C2, V3', 'C2, V4'});
xlabel('PC1');
ylabel('PC2');
zlabel('PC4');
view(-142, -47);