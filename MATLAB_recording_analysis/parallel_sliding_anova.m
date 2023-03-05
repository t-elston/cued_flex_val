function [sig_factors, shuff_sig_factors] = parallel_sliding_anova(firing_rates, factors,fnames)

[n_trials, n_times] = size(firing_rates);

for t = 1:n_times
                        
                
            anova_frs = firing_rates(:,t);
            
            [p,tbl] = anovan(anova_frs,factors,'model','interaction','display','off', 'varnames',fnames);
                    
            [shuffled_p, shuffled_tbl] =  anovan(shuffle(anova_frs),factors,'model','interaction','display','off','varnames',fnames);       
             
            % extract the significant bins
            sig_factors(1,t) = p(1) < .001;
            sig_factors(2,t) = p(2) < .001;
            sig_factors(3,t) = p(4) < .001;
            
            % extract the significant bins
            shuff_sig_factors(1,t) = shuffled_p(1) < .01;
            shuff_sig_factors(2,t) = shuffled_p(2) < .01;
            shuff_sig_factors(3,t) = shuffled_p(4) < .01;
            
  end % of looping over times


end