% load the data
% load('C:\Users\Thomas Elston\Documents\MATLAB\Projects\CuedFlexVal\preprocessed data\Don\population\boot_PCs/popdata.mat');

load('C:\Users\Thomas Elston\Documents\MATLAB\Projects\CuedFlexVal\preprocessed data\King\population\boot_PCs/popdata.mat');

HPC_choice = nanmean(HPC_choice_PCs, 3);
OFC_choice = nanmean(OFC_choice_PCs, 3);

n_boots = size(HPC_choice,3);

hpc_val_dist=[];
ofc_val_dist=[];
hpc_slopes=[];
ofc_slopes=[];
OFC_angles=[];
HPC_angles=[];
OFC_rot_C1=[];
HPC_rot_C1=[];

d_type = 'euclidean';

x=1;
y=2;
z=3;

% let's get some distance stats
for b = 1:n_boots
    
    % CHOICE PERIOD
    % look at distance between same value reps across the contexts
    for v = 1:4
        hpc_val_dist(b,v) = get_dist_v01(HPC_choice_PCs(v,:,b), HPC_choice_PCs(v+4,:,b), d_type);
        ofc_val_dist(b,v) = get_dist_v01(OFC_choice_PCs(v,:,b), OFC_choice_PCs(v+4,:,b), d_type);
    end
    
    % compute slopes of value distances for each bootstrap
    hpc_slopes(b,:) = regress(hpc_val_dist(b,:)', [ones(size(hpc_val_dist(b,:)))', [1:4]']);
    ofc_slopes(b,:) = regress(ofc_val_dist(b,:)', [ones(size(ofc_val_dist(b,:)))', [1:4]']);
    
    % look for rotations
    OFC_ctx1 = [OFC_choice_PCs(1:4,x,b), OFC_choice_PCs(1:4,y,b), OFC_choice_PCs(1:4,z,b)];
    OFC_ctx2 = [OFC_choice_PCs(5:8,x,b), OFC_choice_PCs(5:8,y,b), OFC_choice_PCs(5:8,z,b)];
    
    HPC_ctx1 = [HPC_choice_PCs(1:4,x,b), HPC_choice_PCs(1:4,y,b), HPC_choice_PCs(1:4,z,b)];
    HPC_ctx2 = [HPC_choice_PCs(5:8,x,b), HPC_choice_PCs(5:8,y,b), HPC_choice_PCs(5:8,z,b)];
    
    [OFC_angles(b), OFC_rot_C1(:,:,b)] = check_best_rotation(OFC_ctx1, OFC_ctx2);
    [HPC_angles(b), HPC_rot_C1(:,:,b)] = check_best_rotation(HPC_ctx1, HPC_ctx2);     
    
end % of looping over bootstraps




cmap = cbrewer('qual','Paired',12);

%-----------------------------------------------------
choice_figure = figure; 
set(choice_figure,'Position',[100 100 1000 800], 'renderer','Painters');
%define axes
HPC_gen_PC_ax  = axes('Position',[.15   .7  .2  .2]);
HPC_rot_PC_ax  = axes('Position',[.75   .7  .2  .2]);
OFC_gen_PC_ax  = axes('Position',[.15   .4  .2  .2]);
OFC_rot_PC_ax  = axes('Position',[.75   .4  .2  .2]);

% Choice_dist_ax = axes('Position',[.15   .1  .2  .2]);
% Choice_slope_ax= axes('Position',[.45   .1  .2  .2]);
% angle_ax       = axes('Position',[.75   .1  .2  .2]);

% x = 1; 
% y = 2; 
% z = 3;

axes(HPC_gen_PC_ax);
hold on
plot3(HPC_choice(1,x),HPC_choice(1,y), HPC_choice(1,z), 'Marker','o','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(6,:));
plot3(HPC_choice(2,x),HPC_choice(2,y), HPC_choice(2,z), 'Marker','^','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(6,:));
plot3(HPC_choice(3,x),HPC_choice(3,y), HPC_choice(3,z), 'Marker','s','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(6,:));
plot3(HPC_choice(4,x),HPC_choice(4,y), HPC_choice(4,z), 'Marker','pentagram','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(6,:));
plot3(HPC_choice(1:4,x),HPC_choice(1:4,y), HPC_choice(1:4,z), 'color',cmap(6,:), 'LineWidth', 3);

plot3(HPC_choice(5,x),HPC_choice(5,y), HPC_choice(5,z), 'Marker','o','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(2,:));
plot3(HPC_choice(6,x),HPC_choice(6,y), HPC_choice(6,z), 'Marker','^','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(2,:));
plot3(HPC_choice(7,x),HPC_choice(7,y), HPC_choice(7,z), 'Marker','s','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(2,:));
plot3(HPC_choice(8,x),HPC_choice(8,y), HPC_choice(8,z), 'Marker','pentagram','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(2,:));
plot3(HPC_choice(5:8,x),HPC_choice(5:8,y), HPC_choice(5:8,z), 'color',cmap(2,:), 'LineWidth', 3);

% legend({'C1, V1','C1, V2','C1, V3', 'C1, V4','C2, V1','C2, V2','C2, V3', 'C2, V4'},'Location', [.08 .8 .1 .1]);
title('HPC');
xlabel(['PC' num2str(x)]);
ylabel(['PC' num2str(y)]);
zlabel(['PC' num2str(z)]);
xlim([-18, 18]);
ylim([-18, 18]);
zlim([-18, 18]);
view(10,20); % King
% view(64,68); % Don



rot_HPC = nanmean(HPC_rot_C1,3);
axes(HPC_rot_PC_ax);
hold on
plot3(HPC_choice(5,x),HPC_choice(5,y), HPC_choice(5,z), 'Marker','o','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(2,:));
plot3(HPC_choice(6,x),HPC_choice(6,y), HPC_choice(6,z), 'Marker','^','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(2,:));
plot3(HPC_choice(7,x),HPC_choice(7,y), HPC_choice(7,z), 'Marker','s','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(2,:));
plot3(HPC_choice(8,x),HPC_choice(8,y), HPC_choice(8,z), 'Marker','pentagram','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(2,:));
plot3(HPC_choice(5:8,x),HPC_choice(5:8,y), HPC_choice(5:8,z), 'color',cmap(2,:), 'LineWidth', 3);

plot3(rot_HPC(1,1),rot_HPC(1,2), rot_HPC(1,3), 'Marker','o','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(10,:));
plot3(rot_HPC(2,1),rot_HPC(2,2), rot_HPC(2,3), 'Marker','^','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(10,:));
plot3(rot_HPC(3,1),rot_HPC(3,2), rot_HPC(3,3), 'Marker','s','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(10,:));
plot3(rot_HPC(4,1),rot_HPC(4,2), rot_HPC(4,3), 'Marker','pentagram','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(10,:));
plot3(rot_HPC(1:4,1),rot_HPC(1:4,2), rot_HPC(1:4,3), 'color',cmap(10,:), 'LineWidth', 3);
xlabel(['PC' num2str(x)]);
ylabel(['PC' num2str(y)]);
zlabel(['PC' num2str(z)]);
xlim([-18, 18]);
ylim([-18, 18]);
zlim([-18, 18]);
view(10,20); % King
% view(64,68); % Don
%-----


axes(OFC_gen_PC_ax); 
hold on
plot3(OFC_choice(1,x),OFC_choice(1,y), OFC_choice(1,z), 'Marker','o','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(5,:));
plot3(OFC_choice(2,x),OFC_choice(2,y), OFC_choice(2,z), 'Marker','^','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(5,:));
plot3(OFC_choice(3,x),OFC_choice(3,y), OFC_choice(3,z), 'Marker','s','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(5,:));
plot3(OFC_choice(4,x),OFC_choice(4,y), OFC_choice(4,z), 'Marker','pentagram','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(5,:));
plot3(OFC_choice(1:4,x),OFC_choice(1:4,y), OFC_choice(1:4,z), 'color',cmap(5,:), 'LineWidth', 3);

plot3(OFC_choice(5,x),OFC_choice(5,y), OFC_choice(5,z), 'Marker','o','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(1,:));
plot3(OFC_choice(6,x),OFC_choice(6,y), OFC_choice(6,z), 'Marker','^','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(1,:));
plot3(OFC_choice(7,x),OFC_choice(7,y), OFC_choice(7,z), 'Marker','s','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(1,:));
plot3(OFC_choice(8,x),OFC_choice(8,y), OFC_choice(8,z), 'Marker','pentagram','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(1,:));
plot3(OFC_choice(5:8,x),OFC_choice(5:8,y), OFC_choice(5:8,z), 'color',cmap(1,:), 'LineWidth', 3);
xlabel(['PC' num2str(x)]);
ylabel(['PC' num2str(y)]);
zlabel(['PC' num2str(z)]);
xlim([-20, 20]);
ylim([-20, 20]);
zlim([-20, 20]);
xticks([-10, 0, 10])
yticks([-10, 0, 10])
zticks([-10, 0, 10])
view(-40,10); % King
% view(70,15); % Don



rot_OFC = nanmean(OFC_rot_C1,3);
axes(OFC_rot_PC_ax);
hold on
plot3(OFC_choice(5,x),OFC_choice(5,y), OFC_choice(5,z), 'Marker','o','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(1,:));
plot3(OFC_choice(6,x),OFC_choice(6,y), OFC_choice(6,z), 'Marker','^','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(1,:));
plot3(OFC_choice(7,x),OFC_choice(7,y), OFC_choice(7,z), 'Marker','s','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(1,:));
plot3(OFC_choice(8,x),OFC_choice(8,y), OFC_choice(8,z), 'Marker','pentagram','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(1,:));
plot3(OFC_choice(5:8,x),OFC_choice(5:8,y), OFC_choice(5:8,z), 'color',cmap(1,:), 'LineWidth', 3);

plot3(rot_OFC(1,1),rot_OFC(1,2), rot_OFC(1,3), 'Marker','o','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(9,:));
plot3(rot_OFC(2,1),rot_OFC(2,2), rot_OFC(2,3), 'Marker','^','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(9,:));
plot3(rot_OFC(3,1),rot_OFC(3,2), rot_OFC(3,3), 'Marker','s','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(9,:));
plot3(rot_OFC(4,1),rot_OFC(4,2), rot_OFC(4,3), 'Marker','pentagram','color','k', 'MarkerSize',15, 'MarkerFaceColor',cmap(9,:));
plot3(rot_OFC(1:4,1),rot_OFC(1:4,2), rot_OFC(1:4,3), 'color',cmap(9,:), 'LineWidth', 3);
xlabel(['PC' num2str(x)]);
ylabel(['PC' num2str(y)]);
zlabel(['PC' num2str(z)]);
xlim([-20, 20]);
ylim([-20, 20]);
zlim([-20, 20]);
xticks([-10, 0, 10])
yticks([-10, 0, 10])
zticks([-10, 0, 10])
view(-40,10); % King
% view(70,15); % Don

% % Prepare data to plot value-related distances
% val1_x = make_jittered_xvals(hpc_val_dist(:,1), [.7, 1.3]);
% val2_x = make_jittered_xvals(hpc_val_dist(:,1), [1.7, 2.3]);
% val3_x = make_jittered_xvals(hpc_val_dist(:,1), [2.7, 3.3]);
% val4_x = make_jittered_xvals(hpc_val_dist(:,1), [3.7, 4.3]);
% 
% 
% axes(Choice_slope_ax);
% hold on
% histogram(hpc_slopes(:,2), 'Normalization','Probability','BinWidth',.1,'FaceColor',cmap(10,:),'EdgeColor','none');
% histogram(ofc_slopes(:,2), 'Normalization','Probability','BinWidth',.1,'FaceColor',cmap(9,:),'EdgeColor','none');
% h = gca;
% h.YAxis.Visible = 'off';
% xlabel('Beta(value) (Distance ~ Value + int)');
% % legend({'HPC', 'OFC'},'Location', [.08 .4 .1 .1]);
% 
% axes(angle_ax);
% hold on
% histogram(HPC_angles, 'Normalization','Probability','BinWidth',1,'FaceColor',cmap(10,:),'EdgeColor','none');
% histogram(OFC_angles, 'Normalization','Probability','BinWidth',1,'FaceColor',cmap(9,:),'EdgeColor','none');
% h = gca;
% h.YAxis.Visible = 'off';
% xlabel('Mean Dimension-Wise Rotation (Degrees)');