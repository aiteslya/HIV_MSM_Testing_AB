function [Population, pop_hiv,fl_collect,num_died,alive_ind,steady_dur,casual_dur,Rels_steady,ind_died,infect_alive,reset_casual_prop,total_infect_undiagn_ind,total_diagn_ind,total_treat_ind,total_sup_ind, PrEP_inds, all_casual_parts] = Death_Load(Population,pop_hiv,fl_collect,pop_index,hiv_index,DeathProbs_w,tcounter,alive_ind,steady_dur,casual_dur,Rels_steady,infect_alive,total_infect_undiagn_ind,total_diagn_ind,total_treat_ind,total_sup_ind, all_casual_parts, PrEP_inds)

    % retiring individuals from the population if the yearly mark has been
    % crossed and they have turned 75 or due to the background mortality

    % death due to the background mortality
    % this if close is debugged
    if numel(alive_ind.Data)>0
        mu_arr=DeathProbs_w(Population.Data(pop_index.age_bin,alive_ind.Data)); % assess death rates for everyone, age-dependent
        r1_arr=rand(1,numel(alive_ind.Data));% roll the die
        ind_died=alive_ind.Data(r1_arr<(mu_arr));
    end

    year_week = 52;

    % death due to the being too old
    if mod(tcounter + year_week/2, year_week)==0 % a year has passed
        old_age_remove=alive_ind.Data(find(Population.Data(pop_index.age,alive_ind.Data)>=75));
        ind_died=[ind_died,old_age_remove];
    end

    
    if numel(ind_died)>0
        ind_died = unique(ind_died);
        num_died = numel(ind_died);
        rel_dur=[];
        % handle casual partnerships in the hash map
        for ind=ind_died
            deceased_id = Population.Data(pop_index.id,ind);
            [all_casual_parts, rel_dur_ind] = handleDeath(all_casual_parts, deceased_id, tcounter);
            rel_dur = [rel_dur; rel_dur_ind'];
        end

        casual_dur.Data = [casual_dur.Data; rel_dur];

        [Population,pop_hiv,Rels_steady,steady_dur,alive_ind,infect_alive,reset_casual_prop,total_infect_undiagn_ind,total_diagn_ind,total_treat_ind,total_sup_ind, PrEP_inds] = Remove_List(Population,pop_hiv,Rels_steady,steady_dur, alive_ind,pop_index,hiv_index,ind_died,tcounter,infect_alive,total_infect_undiagn_ind,total_diagn_ind,total_treat_ind,total_sup_ind, PrEP_inds);

    else
        num_died = 0;
        reset_casual_prop=[];
    end

end