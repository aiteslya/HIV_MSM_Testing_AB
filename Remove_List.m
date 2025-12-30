function [Population,pop_hiv,Rels_steady,steady_dur, alive_ind,infect_alive,reset_casual_prop,total_infect_undiagn_ind,total_diagn_ind,total_treat_ind,total_sup_ind, PrEP_inds] = Remove_List(Population,pop_hiv,Rels_steady,steady_dur,alive_ind,pop_index,hiv_index,ind_died,tcounter,infect_alive,total_infect_undiagn_ind,total_diagn_ind,total_treat_ind,total_sup_ind, PrEP_inds)
                                                                                                                                                                      
                                                                                                                                                                      
    % this function removes individuals who died (indices given in ind_died list) from tables, records
    % pertinent statistics and updates the relationship of their partners

    year_week = 52;% duration of year in days
    % remove people from the population and update their records
    alive_ind.Data = setdiff(alive_ind.Data,ind_died);
    PrEP_inds.Data = setdiff(PrEP_inds.Data, ind_died);
    infect_alive.Data=setdiff(infect_alive.Data,ind_died);
    reset_casual_prop=[];

    if numel(ind_died)>0
        for ind=ind_died
            % clean up HIV containers
            if ~isempty(total_infect_undiagn_ind.Data)
                total_infect_undiagn_ind.Data=setdiff(total_infect_undiagn_ind.Data,ind);
            end
            if ~isempty(total_diagn_ind)
                total_diagn_ind=setdiff(total_diagn_ind,ind);
            end
            if ~isempty(total_treat_ind)
                total_treat_ind=setdiff(total_treat_ind,ind);
            end
            if ~isempty(total_sup_ind)
                total_sup_ind=setdiff(total_sup_ind,ind);
            end

            % clean up Population index
            Population.Data(pop_index.death,ind)=(tcounter-mod(tcounter,year_week))/year_week;% update the year of death
            % retrieve their id
            id=Population.Data(pop_index.id,ind);
            % update steady partners
            % retrieve partners ids
            sp1_id=Population.Data(pop_index.sp1_id,ind);
            if sp1_id~=-1 % if there is a partner in slot 1
                % find index of the partner
                ind_sp1=find(Population.Data(pop_index.id,:)==sp1_id);
                reset_casual_prop=[reset_casual_prop, ind_sp1];
                % update the number of steady partnerships of the partner
                Population.Data(pop_index.nsteady,ind_sp1)=Population.Data(pop_index.nsteady,ind_sp1)-1;

                ids=sort([id sp1_id]);
                rel_ind=find(Rels_steady.Data(:,1)==ids(1) & Rels_steady.Data(:,2)==ids(2));
                rec = Rels_steady.Data(rel_ind,:);
                Rels_steady.Data(rel_ind,:)=[];
                 

                steady_dur.Data=[steady_dur.Data; min(tcounter-Population.Data(pop_index.sp1_stdate,ind), rec(:,4)-rec(:,3)) ]; % record relationship duration
                % update the partner's steady relationship record
                Population.Data(pop_index.sp1_id,ind_sp1)=-1;
                Population.Data(pop_index.sp1_stdate,ind_sp1)=-40000;
                % remove this partnership from the the register of the
                % currently ongoing partnerships
                
            end

            % reset entry of the Population table to be used for newly
            % created individual
            Population.Data(pop_index.id, ind) = 0;
            Population.Data(pop_index.birth, ind) = -300;
            Population.Data(pop_index.death, ind) = -1;
            Population.Data(pop_index.age, ind) = -1;
            Population.Data(pop_index.sp1_id, ind) = -1;
            Population.Data(pop_index.sp1_stdate, ind) = -40000;
            Population.Data(pop_index.nsteady, ind) = 0;
            Population.Data(pop_index.casual_prop, ind) = -1;
            Population.Data(pop_index.age_bin, ind) = 0;
            Population.Data(pop_index.steady_prop, ind) = -1;
                                  
            % clean up HIV table
            pop_hiv.Data(hiv_index.id,ind)=0;
            pop_hiv.Data(hiv_index.status,ind)=0;
            pop_hiv.Data(hiv_index.date_inf,ind)=-40000;
            pop_hiv.Data(hiv_index.date_diagn,ind)=-40000;
            pop_hiv.Data(hiv_index.date_art,ind)=-40000;
            pop_hiv.Data(hiv_index.date_sup,ind)=-40000;
            pop_hiv.Data(hiv_index.date_new_stage,ind)=-40000;
            pop_hiv.Data(hiv_index.ever_tested, ind) = -40000;
            pop_hiv.Data(hiv_index.test_rate, ind) = 0;
            pop_hiv.Data(hiv_index.testing_bin, ind) = 0;
        end
    end

    % in the case both partners removed on the same turn no need to update
    % their casual propensities going forward
    reset_casual_prop=setdiff(reset_casual_prop,ind_died);
end