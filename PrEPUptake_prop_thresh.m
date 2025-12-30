function [Pop_hiv, PrEP_inds] = PrEPUptake_prop_thresh(Population, Pop_hiv, hiv_index, HIV_st, alive_ind, infect_alive, cfg, pop_index, PrEP_inds,cfg_PrEP, year_counter)
    % this function simulates enrollment of susceptible individuals into
    % PrEP programme

    % it will need to be updated for the forward simulations: testing
    % should be less frequent

    % Create a logical array indicating whether each element of alive_ind.Data is in infect_alive.Data
    isInfected = ismember(alive_ind.Data, infect_alive.Data);

    % Select elements from alive_ind.Data that are not in infect_alive.Data
    ind_non_infect = alive_ind.Data(~isInfected);

    % Direct logical indexing to find indices where the condition is true
    ind = ind_non_infect(Pop_hiv.Data(hiv_index.status, ind_non_infect) == HIV_st.susc);

    % determine how many people will take up PrEP
    if numel(ind)>0
        year_week = 52;
        prep_up_w = cfg.prep_up/year_week;
        num_new = sum(rand(1,numel(ind))<prep_up_w);
        if num_new>0
           
            % determine who will take PrEP
            % schema: ind_s, age
            
            % find all these whose number of casual partners in the last 6
            % months exceeds the predefined threshold
            temp_ind = find(Population.Data(pop_index.casual_prop,ind)>=cfg.prep_thresh);
            
            eligible_ind = ind(temp_ind);
            % eligible = [eligible_ind; Population.Data(pop_index.age_bin, eligible_ind)]';

            if numel(eligible_ind)>0 
                
                if num_new < numel(eligible_ind)
                    ind_sel = datasample(eligible_ind,num_new);
                else
                    ind_sel = eligible_ind;
                end

                ind_pr = ind_sel;
                Pop_hiv.Data(hiv_index.status,ind_pr) = HIV_st.suscPrep;
                PrEP_inds.Data = [PrEP_inds.Data, ind_pr];

                % set testing frequency
                if year_counter <2025
                    Pop_hiv.Data(hiv_index.test_rate, ind_pr) = cfg_PrEP.old_rate;
                else
                    Pop_hiv.Data(hiv_index.test_rate, ind_pr) = cfg_PrEP.new_rate;
                end
            end
        end
    end

end