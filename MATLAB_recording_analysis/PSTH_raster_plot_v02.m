function PSTH_raster_plot_v02(bhv,firing_rates,ts,  rasters, raster_ts, notes, cue_pval, choice_pvals)

raster_ts = raster_ts(1:numel(rasters(1,:)));
state1_cmap = cbrewer('seq','Blues',7);
state2_cmap = cbrewer('seq','Reds',7);
val_cmap = cbrewer('seq','BuPu',10);

unit_info{1,1} = ['Region: ' notes.Target{1}];
unit_info{2,1} = [notes.Filename{1}(1:12) '_' num2str(notes.Channel(1)) notes.Unit{1}];

% make some dummy timing arrays
cue_times_ix = double(ts >= -650 & ts <= -300);
choice_time_ix = double(ts >= 10 & ts <= 400);


cue_times_ix(cue_times_ix==0) = NaN;
choice_time_ix(choice_time_ix==0) = NaN;

cue_state_sig = double(cue_times_ix * cue_pval(1) < .05);
choice_state_sig = double(choice_time_ix * choice_pvals(1) < .05);
choice_val_sig = double(choice_time_ix * choice_pvals(2) < .05);
choice_stateXval_sig = double(choice_time_ix * choice_pvals(3) < .05);

cue_state_sig(cue_state_sig==0) = NaN;
choice_state_sig(choice_state_sig==0) = NaN;
choice_val_sig(choice_val_sig==0) = NaN;
choice_stateXval_sig(choice_stateXval_sig==0) = NaN;




% collect the trial indices
s1_ix = bhv.state == 1;
s2_ix = bhv.state == 2;
type1_ix = bhv.state_type == 1;
type2_ix = bhv.state_type == 2;
v1_ix = bhv.chosenval == 1;
v2_ix = bhv.chosenval == 2;
v3_ix = bhv.chosenval == 3;
v4_ix = bhv.chosenval == 4;

% let's get the means for each group
s1_t1_mean = nanmean(firing_rates(s1_ix & type1_ix,:),1);
s1_t2_mean = nanmean(firing_rates(s1_ix & type2_ix,:),1);
s2_t1_mean = nanmean(firing_rates(s2_ix & type1_ix,:),1);
s2_t2_mean = nanmean(firing_rates(s2_ix & type2_ix,:),1);

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
s1_t1_rasters = get_raster_trials(n_raster_trials*2, rasters, s1_ix & type1_ix);
s1_t2_rasters = get_raster_trials(n_raster_trials*2, rasters, s1_ix & type2_ix);
s2_t1_rasters = get_raster_trials(n_raster_trials*2, rasters, s2_ix & type1_ix);
s2_t2_rasters = get_raster_trials(n_raster_trials*2, rasters, s2_ix & type2_ix);

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
s1_t2_rasters(~isnan(s1_t2_rasters)) = s1_t2_rasters(~isnan(s1_t2_rasters)) + max(s1_t1_rasters,[],'all');
s2_t1_rasters(~isnan(s2_t1_rasters)) = s2_t1_rasters(~isnan(s2_t1_rasters)) + max(s1_t2_rasters,[],'all');
s2_t2_rasters(~isnan(s2_t2_rasters)) = s2_t2_rasters(~isnan(s2_t2_rasters)) + max(s2_t1_rasters,[],'all');

plot(raster_ts,s1_t1_rasters','.','color',state1_cmap(3,:));
plot(raster_ts,s1_t2_rasters','.','color',state1_cmap(6,:));
plot(raster_ts,s2_t1_rasters','.','color',state2_cmap(3,:));
plot(raster_ts,s2_t2_rasters','.','color',state2_cmap(6,:));
axis tight
plot([0 0], ylim,'k','LineWidth',2);
plot([-700 -700], ylim,'k','LineWidth',2);
xlim([-1000 1500]);
set(gca,'XColor','none','YColor','none');
title('State');


axes(state_mean_ax);
hold on
plot(ts,s1_t1_mean,'LineWidth',2,'color',state1_cmap(3,:));
plot(ts,s1_t2_mean,'LineWidth',2,'color',state1_cmap(6,:));
plot(ts,s2_t1_mean,'LineWidth',2,'color',state2_cmap(3,:));
plot(ts,s2_t2_mean,'LineWidth',2,'color',state2_cmap(6,:));
plot(ts,cue_state_sig.*min(ylim) + .2*min(ylim),'k','Marker','.','MarkerSize',15);
plot(ts,choice_state_sig.*min(ylim) + .2*min(ylim),'k','Marker','.','MarkerSize',15);
axis tight
xlim([-1000 1500]);
plot([0 0], ylim,'k','LineWidth',2);
plot([-700 -700], ylim,'k','LineWidth',2);
ylabel('firing rate (Hz)');
xlabel('time from choice');
legend({'S1, C1','S1, C2', 'S2, C1','S2, C2'},'NumColumns',2);



axes(val_raster_ax);
hold on
plot(raster_ts,v1_rasters','.','color',val_cmap(4,:));

v2_rasters(~isnan(v2_rasters)) = v2_rasters(~isnan(v2_rasters)) + max(v1_rasters,[],'all');
v3_rasters(~isnan(v3_rasters)) = v3_rasters(~isnan(v3_rasters)) + max(v2_rasters,[],'all');
v4_rasters(~isnan(v4_rasters)) = v4_rasters(~isnan(v4_rasters)) + max(v3_rasters,[],'all');
plot(raster_ts,v2_rasters','.','color',val_cmap(6,:));
plot(raster_ts,v3_rasters','.','color',val_cmap(8,:));
plot(raster_ts,v4_rasters','.','color',val_cmap(10,:));
axis tight
set(gca,'XColor','none','YColor','none');
plot([0 0], ylim,'k','LineWidth',2);
plot([-700 -700], ylim,'k','LineWidth',2);
xlim([-1000 1500]);
title('Value');


axes(val_mean_ax);
hold on
plot(ts,v1_mean,'LineWidth',2,'color',val_cmap(4,:));
plot(ts,v2_mean,'LineWidth',2,'color',val_cmap(6,:));
plot(ts,v3_mean,'LineWidth',2,'color',val_cmap(8,:));
plot(ts,v4_mean,'LineWidth',2,'color',val_cmap(10,:));
plot(ts,choice_val_sig.*min(ylim) + .2*min(ylim),'k','Marker','.','MarkerSize',15);
axis tight
xlim([-1000 1500]);
plot([0 0], ylim,'k','LineWidth',2);
plot([-700 -700], ylim,'k','LineWidth',2);
legend('Val 1','Val 2', 'Val 3','Val 4');



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
set(gca,'XColor','none','YColor','none');
plot([0 0], ylim,'k','LineWidth',2);
plot([-700 -700], ylim,'k','LineWidth',2);
xlim([-1000 1500]);
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
plot(ts,choice_stateXval_sig.*min(ylim) + .2*min(ylim),'k','Marker','.','MarkerSize',15);
axis tight
xlim([-1000 1500]);
plot([0 0], ylim,'k','LineWidth',2);
plot([-700 -700], ylim,'k','LineWidth',2);
legend({'s1 v1','s1 v2', 's1 v3','s1 v4','s2 v1','s2 v2', 's2 v3',' s2 v4'},'NumColumns',2);


a = annotation('textbox',[.1,.75,.3,.2],'String',unit_info,'FitBoxToText','on','Interpreter', 'none');




end % of function 