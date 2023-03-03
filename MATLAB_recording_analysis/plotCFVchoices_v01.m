function plotCFVchoices_v01(bhv)


bhv(bhv.forcedchoice==1,:)=[];

% define color scheme
ct = cbrewer('qual','Set1',9);




% make a panel looking at comparing each thing against each other
lvals = unique(bhv.Lval);
rvals = unique(bhv.Rval);


choicegrid = NaN(numel(rvals), numel(lvals), 2);


for s_ix = 1:2
for r_ix = 1:numel(rvals)
    
    rval = rvals(r_ix);
    
    for l_ix = 1:numel(lvals)
        
        lval = lvals(l_ix);
        
        choicegrid(r_ix,l_ix,s_ix) = nanmean(bhv.pickedbest(bhv.Rval == rval & bhv.Lval == lval & bhv.state==s_ix));
        
    end % of cycling through lvals
    
end % of cycling through rvals
end % of cycling through states

cmap = viridis(100);

figure; 
subplot(1, 2, 1); 
imagesc(choicegrid(:,:,1));
colormap(cmap);
caxis([0,1]);
xlabel('Right Value');
ylabel('Left Value');
xticks([1:4]);
yticks([1:4]);
cb = colorbar;
cb.FontSize = 12;
cb.Position = [.1 .73 .01 .2];
set(gca,'FontSize',12,'LineWidth',1);
title('Context A');
axis square

subplot(1, 2, 2); 
imagesc(choicegrid(:,:,2));
colormap(cmap);
caxis([0,1]);
xlabel('Right Value');
ylabel('Left Value');
xticks([1:4]);
yticks([1:4]);
set(gca,'FontSize',12,'LineWidth',1);
title('Context B');
axis square

figure;
% get group level summary
[r1_means,r1_sems] = grpstats(bhv.pickedleft(bhv.state==1),bhv.optdiff(bhv.state==1),{'mean','sem'});
[r2_means,r2_sems] = grpstats(bhv.pickedleft(bhv.state==2),bhv.optdiff(bhv.state==2),{'mean','sem'});
hold on
errorbar(unique(bhv.optdiff),r1_means,r1_sems,'color',ct(1,:),'marker','.',...
    'MarkerSize',30,'LineWidth',3);
errorbar(unique(bhv.optdiff),r2_means,r2_sems,'color',ct(2,:),'marker','.',...
    'MarkerSize',30,'LineWidth',3);

xticks([-3:3]);
set(gca,'FontSize',12,'LineWidth',1);
xlabel('Lval - Rval');
ylabel('p(Choose Left)');
legend('Context A','Context B');
axis square










hits = bhv(bhv.pickedbest ==1,:);


figure;
hold on
histogram(hits.rt((hits.Lval + hits.Rval) > 4),'BinWidth',20,'Normalization','probability')
histogram(hits.rt((hits.Lval + hits.Rval) < 4),'BinWidth',20,'Normalization','probability')
xlabel('rt (ms)');
legend('+ Best Opt','- Best Opt');
set(gca,'LineWidth',1,'FontSize',12);












end % of function