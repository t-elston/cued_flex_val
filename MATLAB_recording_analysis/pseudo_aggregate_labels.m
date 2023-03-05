function [ix_labels, ix_frs] = pseudo_aggregate_labels(ix, bhv, firing_rates)

ix_labels=[];
ix_frs=[];

ix_trials = bhv(ix, :);
ix_firing_rates = firing_rates(ix, :,:);

states = unique(bhv.state);
state_types = unique(bhv.state_type);
vals = unique(bhv.chosenval);

for s = 1:numel(states)
        for v = 1:numel(vals)
            
            % find the trials which match this current constellations of
            % conditions
            
            s_t_v_ix = (ix_trials.state == states(s)) & (ix_trials.chosenval == vals(v));
            
            n_trials = sum(s_t_v_ix);
            
            %TODO 
            % make chunks of labels and firing rates that reflect this
            % specific constellation of factors
            % aggregate those over constellations
            this_label_constellation = [ix_trials.state(s_t_v_ix) , ix_trials.state_type(s_t_v_ix) , ix_trials.chosenval(s_t_v_ix) ] ;
            
            ix_labels = [ix_labels ; this_label_constellation];
            ix_frs  = cat(1,ix_frs, ix_firing_rates(s_t_v_ix,:,:));
            
            
        end % of looping over values
end % of looping over states



end % of function