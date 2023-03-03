% aaa_CFV_pop_trajectory_main_v01.m

datadir = 'C:\Users/Thomas Elston/Documents/PYTHON/CuedFlexVal/population_data/';

f_data = dir([datadir '*.mat']);

n_boots = 100;

% load data
load([datadir 'HPC_popdata']);

bhv = array2table(bhv,'VariableNames',bhv_varnames);
u_ids = unique(u_names);
n_units = numel(u_ids);
n_times = size(zFRs,2);

% create label vector
ctx_label = [ones(n_times,1) ; ones(n_times,1)  ; ones(n_times,1)  ; ones(n_times,1)  ;
             ones(n_times,1) *-1 ; ones(n_times,1) *-1 ; ones(n_times,1) *-1 ; ones(n_times,1) *-1];    
         
val_label = [ones(n_times,1) ; ones(n_times,1)*2 ; ones(n_times,1)*3 ; ones(n_times,1)*4  ;
             ones(n_times,1) ; ones(n_times,1)*2 ; ones(n_times,1)*3 ; ones(n_times,1)*4];   
         
for b = 1:n_boots
    
    all_FRs = [];
    
    u_ctr = 0;
    
    for u = 1:numel(u_ids)
        
        u_ix = contains(u_names, u_ids(u));
        
        % grab 80% of the trials for each condition
        [c1_v1] = get_random_proportion_v02(u_ix & bhv.state == 1 & bhv.chosenval == 1,.8);
        [c1_v2] = get_random_proportion_v02(u_ix & bhv.state == 1 & bhv.chosenval == 2,.8);
        [c1_v3] = get_random_proportion_v02(u_ix & bhv.state == 1 & bhv.chosenval == 3,.8);
        [c1_v4] = get_random_proportion_v02(u_ix & bhv.state == 1 & bhv.chosenval == 4,.8);
        [c2_v1] = get_random_proportion_v02(u_ix & bhv.state == -1 & bhv.chosenval == 1,.8);
        [c2_v2] = get_random_proportion_v02(u_ix & bhv.state == -1 & bhv.chosenval == 2,.8);
        [c2_v3] = get_random_proportion_v02(u_ix & bhv.state == -1 & bhv.chosenval == 3,.8);
        [c2_v4] = get_random_proportion_v02(u_ix & bhv.state == -1 & bhv.chosenval == 4,.8);
        
        if ~isempty(c1_v1) & ~isempty(c1_v2) & ~isempty(c1_v3) & ~isempty(c1_v4) & ...
           ~isempty(c2_v1) & ~isempty(c2_v2) & ~isempty(c2_v3) & ~isempty(c2_v4)
       
            u_ctr = u_ctr+1;
            all_FRs(:, u_ctr) = [nanmean(zFRs(c1_v1, :),1)' ; 
                                 nanmean(zFRs(c1_v2, :),1)' ;
                             nanmean(zFRs(c1_v3, :),1)' ;
                             nanmean(zFRs(c1_v4, :),1)' ;
                             nanmean(zFRs(c2_v1, :),1)' ; 
                             nanmean(zFRs(c2_v2, :),1)' ;
                             nanmean(zFRs(c2_v3, :),1)' ;
                             nanmean(zFRs(c2_v4, :),1)'];
                         
        end % of ensuring each neuron was represented in each condition
                     
    end % of looping over neurons
    
    % compute the PCs
    %[~,PCs] = pca(zscore(all_FRs,[],1), 'NumComponents',10);  
    [~,PCs] = pca(all_FRs, 'NumComponents',10);  


    ctx_val_distances(b,:) = mean(reshape(sqrt(sum((PCs(ctx_label==1,:) - PCs(ctx_label==-1,:)).^2, 2))',n_times,[]),2);

end % of looping over bootstraps

ctx1_PCs = mean(reshape(PCs(ctx_label==1,:),n_times,[],4),3);
ctx2_PCs = mean(reshape(PCs(ctx_label==-1,:),n_times,[],4),3);

x = 1;
y = 4;
z = 3; 
figure;
hold on
plot3(ctx1_PCs(:,x), ctx1_PCs(:,y), ctx1_PCs(:,z), 'LineWidth',3)
plot3(ctx2_PCs(:,x), ctx2_PCs(:,y), ctx2_PCs(:,z), 'LineWidth',3)
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');


[ctx_mean, ctx_ci] = GetMeanCI(ctx_val_distances,'percentile');

figure;
plot(ts, ctx_mean, 'LineWidth',3)


