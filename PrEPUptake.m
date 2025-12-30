function [Pop_hiv, PrEP_inds]=PrEPUptake(Population,Pop_hiv,hiv_index,HIV_st,alive_ind,infect_alive,cfg,prep_uptake,pop_index, rec_casual_parts, PrEP_inds)
    % this function simulates enrollment of susceptible individuals into
    % PrEP programme

    
    % Create a logical array indicating whether each element of alive_ind.Data is in infect_alive.Data
    isInfected = ismember(alive_ind.Data, infect_alive.Data);

    % Select elements from alive_ind.Data that are not in infect_alive.Data
    ind_non_infect = alive_ind.Data(~isInfected);

    % Direct logical indexing to find indices where the condition is true
    ind = ind_non_infect(Pop_hiv.Data(hiv_index.status, ind_non_infect) == HIV_st.susc);

    if numel(ind)>0
        % duration of PrEP in weeks
        year_week = 52;
        prep_up_d=cfg.prep_up/year_week;
        % determine how many people will take up PrEP
        num_new=sum(rand(1,numel(ind))<prep_up_d);
        if num_new>0
            % determine who will take PrEP
            % schema: ind_s, age
            eligible=[];

            counter=1;
            while counter<=numel(ind) & floor(numel(eligible)/2)<num_new
                
                % determine whether this person is eligible by accessing
                % their number of partners within  the last year
                id=Population.Data(pop_index.id,ind(counter));

                % num_casual_parts=numel(retrievePartners(rec_casual_parts,id));

                %if num_casual_parts>=cfg.prep_thresh % number of casual partners exceeds rate
                if Population.Data(pop_index.casual_prop,ind(counter))>=cfg.prep_thresh % propensity to acquire casual partners exceeds threshold
                    % get age to seed the probability of taking prep
                    age_bin=Population.Data(pop_index.age_bin,ind(counter));
                    pre_up_d=prep_uptake(age_bin);
                    eligible=[eligible; ind(counter) pre_up_d];
                end
                counter=counter+1;
            end

            if floor(numel(eligible)/2)>0 
                
                if numel(eligible(:,1))>num_new
                    ind_pr = randsample(eligible(:,1),num_new,true,eligible(:,2)./(sum(eligible(:,2))));
                else
                    ind_pr = eligible(:,1);
                end
                
                Pop_hiv.Data(hiv_index.status,ind_pr)=HIV_st.suscPrep;
                Pop_hiv.Data(hiv_index.test_rate, ind_pr) = cfg.prep_test_rate;
                Pop_hiv.Data(hiv_index.testing_bin, ind_pr) = 5;
                PrEP_inds.Data = [PrEP_inds.Data, ind_pr'];
            end
        end
    end

end