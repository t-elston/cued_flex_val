function PSTH_raster_plot(bhv,firing_rates,ts,  rasters, raster_ts, notes, siginfo)

raster_ts = raster_ts(1:numel(rasters(1,:)));
state1_cmap = cbrewer('seq','Blues',5);
state2_cmap = cbrewer('seq','Reds',5);
val_cmap = cbrewer('seq','Purples',8);

unit_info{1,1} = ['Region: ' notes.Target{1}];
unit_info{2,1} = [notes.Filename{1}(1:12) '_' num2str(notes.Channel(1)) notes.Unit{1}];

siginfo(siginfo==0) = NaN;

% collect the trial indices
s1_ix = bhv.state == 1;
s2_ix = bhv.state == 2;
v1_ix = bhv.chosenval == 1;
v2_ix = bhv.chosenval == 2;
v3_ix = bhv.chosenval == 3;
v4_ix = bhv.chosenval == 4;

% let's get the means for each group
s1_mean = nanmean(firing_rates(s1_ix,:),1);
s2_mean = nanmean(firing_rates(s2_ix,:),1);

s1_v1_mean = nanmean(firing_rates(s1_ix & v1_ix,:),1);
s1_v2_mean = nanmean(firing_rates(s1_ix & v2_ix,:),1);
s1_v3_mean = nanmean(firing_rates(s1_ix & v3_ix,:),1);
s1_v4_mean = nanmean(firing_rates(s1_ix & v4_ix,:),1);

s2_v1_mean = nanmean(firing_rates(s2_ix & v1_ix,:),1);
s2_v2_mean = nanmean(firing_rates(s2_ix & v2_ix,:),1);
s2_v3_mean = nanmean(firing_rates(s2_ix & v3_ix,:),1);
s2_v4_mean = nanmean(firing_rates(s2_ix & v4_ix,:),1);

v1_mean = nanmean(firing_rates(v1_ix,:),1);
v2_mean = nanmean(firing_rates(v2_ix,:),1);
v3_mean = nanmean(firing_rates(v3_ix,:),1);
v4_mean = nanmean(firing_rates(v4_ix,:),1);

% collect the rasters
n_raster_trials = 20;
s1_rasters = get_raster_trials(n_raster_trials*4, rasters, s1_ix);
s2_rasters = get_raster_trials(n_raster_trials*4, rasters, s2_ix);
v1_rasters = get_raster_trials(n_raster_trials*2, rasters, v1_ix);
v2_rasters = get_raster_trials(n_raster_trials*2, rasters, v2_ix);
v3_rasters = get_raster_trials(n_raster_trials*2, rasters, v3_ix);
v4_rasters = get_raster_trials(n_raster_trials*2, rasters, v4_ix);
s1_v1_rasters = get_raster_trials(n_raster_trials, rasters, s1_ix & v1_ix);
s1_v2_rasters = get_raster_trials(n_raster_trials, rasters, s1_ix & v2_ix);
s1_v3_rasters = get_raster_trials(n_raster_trials, rasters, s1_ix & v3_ix);
s1_v4_rasters = get_raster_trials(n_raster_trials, rasters, s1_ix & v4_ix);
s2_v1_rasters = get_raster_trials(n_raster_trials, rasters, s2_ix & v1_ix);
s2_v2_rasters = get_raster_trials(n_raster_trials, rasters, s2_ix & v2_ix);
s2_v3_rasters = get_raster_trials(n_raster_trials, rasters, s2_ix & v3_ix);
s2_v4_rasters = get_raster_trials(n_raster_trials, rasters, s2_ix & v4_ix);

% create the figure and axes
fig = figure;
set(fig,'Position',[100 100 1500 600]);

state_raster_ax = axes('Position',[.1 .55 .25 .3]);
state_mean_ax   = axes('Position',[.1 .15 .25 .4]);

val_raster_ax = axes('Position',[.4 .55 .25 .3]);
val_mean_ax   = axes('Position',[.4 .15 .25 .4]);

stateXval_raster_ax = axes('Position',[.7 .55 .25 .3]);
stateXval_mean_ax   = axes('Position',[.7 .15 .25 .4]);


axes(state_raster_ax);
hold on
plot(raster_ts,s1_rasters','.','color',state1_cmap(5,:));
s2_rasters(~isnan(s2_rasters)) = s2_rasters(~isnan(s2_rasters)) + max(s1_rasters,[],'all');
plot(raster_ts,s2_rasters','.','color',state2_cmap(5,:));
axis tight
title('State');


axes(state_mean_ax);
hold on
plot(ts,s1_mean,'LineWidth',2,'color',state1_cmap(5,:));
plot(ts,s2_mean,'LineWidth',2,'color',state2_cmap(5,:));
plot(ts,siginfo(1,:).*2,'k','Marker','.','MarkerSize',15);



axes(val_raster_ax);
hold on
plot(raster_ts,v1_rasters','.','color',val_cmap(5,:));

v2_rasters(~isnan(v2_rasters)) = v2_rasters(~isnan(v2_rasters)) + max(v1_rasters,[],'all');
v3_rasters(~isnan(v3_rasters)) = v3_rasters(~isnan(v3_rasters)) + max(v2_rasters,[],'all');
v4_rasters(~isnan(v4_rasters)) = v4_rasters(~isnan(v4_rasters)) + max(v3_rasters,[],'all');
plot(raster_ts,v2_rasters','.','color',val_cmap(6,:));
plot(raster_ts,v3_rasters','.','color',val_cmap(7,:));
plot(raster_ts,v4_rasters','.','color',val_cmap(8,:));
axis tight
title('Value');


axes(val_mean_ax);
hold on
plot(ts,v1_mean,'LineWidth',2,'color',val_cmap(5,:));
plot(ts,v2_mean,'LineWidth',2,'color',val_cmap(6,:));
plot(ts,v3_mean,'LineWidth',2,'color',val_cmap(7,:));
plot(ts,v4_mean,'LineWidth',2,'color',val_cmap(8,:));
plot(ts,siginfo(2,:).*2,'k','Marker','.','MarkerSize',15);


axes(stateXval_raster_ax);
hold on
s1_v2_rasters(~isnan(s1_v2_rasters)) = s1_v2_rasters(~isnan(s1_v2_rasters)) + max(s1_v1_rasters,[],'all');
s1_v3_rasters(~isnan(s1_v3_rasters)) = s1_v3_rasters(~isnan(s1_v3_rasters)) + max(s1_v2_rasters,[],'all');
s1_v4_rasters(~isnan(s1_v4_rasters)) = s1_v4_rasters(~isnan(s1_v4_rasters)) + max(s1_v3_rasters,[],'all');
plot(raster_ts,s1_v1_rasters','.','color',state1_cmap(2,:));
plot(raster_ts,s1_v2_rasters','.','color',state1_cmap(3,:));
plot(raster_ts,s1_v3_rasters','.','color',state1_cmap(4,:));
plot(raster_ts,s1_v4_rasters','.','color',state1_cmap(5,:));

s2_v1_rasters(~isnan(s2_v1_rasters)) = s2_v1_rasters(~isnan(s2_v1_rasters)) + max(s1_v4_rasters,[],'all');
s2_v2_rasters(~isnan(s2_v2_rasters)) = s2_v2_rasters(~isnan(s2_v2_rasters)) + max(s2_v1_rasters,[],'all');
s2_v3_rasters(~isnan(s2_v3_rasters)) = s2_v3_rasters(~isnan(s2_v3_rasters)) + max(s2_v2_rasters,[],'all');
s2_v4_rasters(~isnan(s2_v4_rasters)) = s2_v4_rasters(~isnan(s2_v4_rasters)) + max(s2_v3_rasters,[],'all');
plot(raster_ts,s2_v1_rasters','.','color',state2_cmap(2,:));
plot(raster_ts,s2_v2_rasters','.','color',state2_cmap(3,:));
plot(raster_ts,s2_v3_rasters','.','color',state2_cmap(4,:));
plot(raster_ts,s2_v4_rasters','.','color',state2_cmap(5,:));
axis tight
title('State x Value');



axes(stateXval_mean_ax);
hold on
plot(ts,s1_v1_mean,'LineWidth',2,'color',state1_cmap(2,:));
plot(ts,s1_v2_mean,'LineWidth',2,'color',state1_cmap(3,:));
plot(ts,s1_v3_mean,'LineWidth',2,'color',state1_cmap(4,:));
plot(ts,s1_v4_mean,'LineWidth',2,'color',state1_cmap(5,:));
plot(ts,s2_v1_mean,'LineWidth',2,'color',state2_cmap(2,:));
plot(ts,s2_v2_mean,'LineWidth',2,'color',state2_cmap(3,:));
plot(ts,s2_v3_mean,'LineWidth',2,'color',state2_cmap(4,:));
plot(ts,s2_v4_mean,'LineWidth',2,'color',state2_cmap(5,:));
plot(ts,siginfo(3,:).*2,'k','Marker','.','MarkerSize',15);

a = annotation('textbox',[.1,.75,.3,.2],'String',unit_info,'FitBoxToText','on','Interpreter', 'none');




end % of function 