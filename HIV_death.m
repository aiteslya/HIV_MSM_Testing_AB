function [Population,fl_collect,Pop_hiv,Rels_steady,steady_dur,casual_dur,alive_ind,num_died_HIV,ind_dead,infect_alive,total_infect_undiagn_ind,total_diagn_ind,total_treat_ind,total_sup_ind,reset_casual_prop,all_casual_parts, PrEP_inds]=HIV_death(Population,fl_collect,Pop_hiv,Rels_steady,steady_dur,casual_dur,alive_ind,pop_index,hiv_index,cfg,tcounter,infect_alive,HIV_st,total_infect_undiagn_ind,total_diagn_ind,total_treat_ind,total_sup_ind,all_casual_parts, PrEP_inds)
% this function performs the module that commits HIV-related death events
% the parameters passed in cfg structure are probabilities per year, need to
% be converted to probabilities per week
    year_week=52;
    % conversion of probabilities from per year to per day
    mu_eaD=cfg.mu_ea/year_week;
    mu_laD=cfg.mu_la/year_week;

    % find individuals in early AIDS stage
    % this is already calculated once prior to execution, so it just can be
    % passed.
    ind_ea1=alive_ind.Data(find(Pop_hiv.Data(hiv_index.status,alive_ind.Data)==HIV_st.EA));
    ind_ea2=alive_ind.Data(find(Pop_hiv.Data(hiv_index.status,alive_ind.Data)==HIV_st.EAD));
    ind_ea3=alive_ind.Data(find(Pop_hiv.Data(hiv_index.status,alive_ind.Data)==HIV_st.EAT));
    ind_ea=sort([ind_ea1,ind_ea2,ind_ea3]);

    ind_dead1=ind_ea(rand(1,numel(ind_ea))<mu_eaD);

    % find individuals in late hiv stage
    ind_la1=alive_ind.Data(find(Pop_hiv.Data(hiv_index.status,alive_ind.Data)==HIV_st.LA));
    ind_la2=alive_ind.Data(find(Pop_hiv.Data(hiv_index.status,alive_ind.Data)==HIV_st.LAD));
    ind_la3=alive_ind.Data(find(Pop_hiv.Data(hiv_index.status,alive_ind.Data)==HIV_st.LAT));
    ind_la=sort([ind_la1,ind_la2,ind_la3]);

    ind_dead2=ind_la(rand(1,numel(ind_la))<mu_laD);

    ind_dead=sort([ind_dead1,ind_dead2]);

    num_died_HIV=numel(ind_dead);
    if num_died_HIV>0
        % handle casual partnerships in the hash map
        rel_dur = [];
        for ind=ind_dead
            deceased_id = Population.Data(pop_index.id,ind);
            [all_casual_parts, rel_dur_ind] = handleDeath(all_casual_parts, deceased_id, tcounter);
            rel_dur = [rel_dur; rel_dur_ind'];
        end

        casual_dur.Data = [casual_dur.Data; rel_dur];

        [Population,Pop_hiv,Rels_steady,steady_dur,alive_ind,infect_alive,reset_casual_prop,total_infect_undiagn_ind,total_diagn_ind,total_treat_ind,total_sup_ind, PrEP_inds] = Remove_List(Population,Pop_hiv,Rels_steady,steady_dur, alive_ind,pop_index,hiv_index,ind_dead,tcounter,infect_alive,total_infect_undiagn_ind,total_diagn_ind,total_treat_ind,total_sup_ind, PrEP_inds);
    else
        reset_casual_prop=[];
    end
end
