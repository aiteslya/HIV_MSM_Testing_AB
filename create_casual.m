function [Pop_hiv, infect_alive, infect_source_data,infect_target_data, all_casual_parts, rec_casual_parts, new_rels_casual, newly_infect1, infectionLog]=create_casual(Population,Pop_hiv,pop_index,tcounter,alive_ind,casual_dur_distr,age_attr,sero_attr,condom,burn_fl,infect_source_data,infect_target_data,infect_alive,hiv_index,cfg,HIV_st, all_casual_parts, rec_casual_parts, infectionLog)
                                                                                                                              
% create new casual pairs (per day)
ind_avail = alive_ind.Data(find(Population.Data(pop_index.casual_prop,alive_ind.Data)>0));

% determine how many pairs can be created

num_pairs=numel(alive_ind.Data);%= numel(ind_avail);

year_week = 52; % duration of the year in weeks

sigma_w=cfg.sigma_cas/year_week; % convert from probability per year to probability per week

% non-dimensionalize the distribution of relationships durations
casual_dur_distr = casual_dur_distr/sum(casual_dur_distr);
cum_casual_dur_distr = cumsum(casual_dur_distr);

% determine how many of the couples form

r1_arr=rand(1,num_pairs);
num_pairs_real=sum(r1_arr<sigma_w);

% array to keep track of couples created in this step for whom the first
% sexual contact has already taken place
new_rels_casual = [];

newly_infect1 = [];
% disp('create of new casual partnerships')
% tic

[all_casual_parts, rec_casual_parts, list_ids] = create_casual_pair(Population,pop_index,tcounter,ind_avail,age_attr,sero_attr,Pop_hiv,hiv_index, cum_casual_dur_distr, all_casual_parts, rec_casual_parts,num_pairs_real);

%toc

if ~isempty(list_ids)
    new_rels_casual = [new_rels_casual; list_ids];

    if burn_fl % burn in is over simulate the infection
       % simulate infection
        %n_pairs = size(list_ids,1);

        % disp('search for couples with transmission potential')
        % % 
        % tic

        list_can_trans = [HIV_st.A1 HIV_st.A23 HIV_st.A45 HIV_st.C HIV_st.EA HIV_st.LA HIV_st.A1D HIV_st.A23D HIV_st.A45D HIV_st.CD HIV_st.EAD HIV_st.LAD HIV_st.A1T HIV_st.A23T HIV_st.A45T HIV_st.CT HIV_st.EAT HIV_st.LAT];
      % Extract the column of IDs from Pop_hiv.Data for comparison
        id_data = Pop_hiv.Data(hiv_index.id,:);
        
        % Find indices of list_ids in id_data using ismember
        [~, loc1] = ismember(list_ids(:,1), id_data);
        
        % The 'loc' output from ismember gives the indices of list_ids in id_data
        % where there is a match, if no match exists, loc will be 0 (since list_ids are unique and
        % are guaranteed to be in id_data, this shouldn't happen unless there's a data integrity issue)
        
        % Store these indices in list_inds
        list_inds(:,1) = loc1;

         % Find indices of list_ids in id_data using ismember
        [~, loc2] = ismember(list_ids(:,2), id_data);
        
        % The 'loc' output from ismember gives the indices of list_ids in id_data
        % where there is a match, if no match exists, loc will be 0 (since list_ids are unique and
        % are guaranteed to be in id_data, this shouldn't happen unless there's a data integrity issue)
        
        % Store these indices in list_inds
        list_inds(:,2) = loc2;

        %is_relevant_id1 = ismember(Pop_hiv.Data(hiv_index.id,:), relevant_ids1);
        relevant_indices1 = list_inds(:,1);
        relevant_indices2 = list_inds(:,2);
        

        % Logical vectors for transmission capability and susceptibility using relevant indices
        can_transmit1 = ismember(Pop_hiv.Data(hiv_index.status, relevant_indices1), list_can_trans);
        is_susceptible2 = (Pop_hiv.Data(hiv_index.status, relevant_indices2) == HIV_st.susc) | (Pop_hiv.Data(hiv_index.status, relevant_indices2) == HIV_st.suscPrep);
        
        % Logical vectors for transmission capability and susceptibility using relevant indices
        can_transmit2 = ismember(Pop_hiv.Data(hiv_index.status, relevant_indices2), list_can_trans);
        is_susceptible1 = (Pop_hiv.Data(hiv_index.status, relevant_indices1) == HIV_st.susc) | (Pop_hiv.Data(hiv_index.status, relevant_indices1) == HIV_st.suscPrep);

        % Apply logical AND for each condition pair
        cond1 = can_transmit1 & is_susceptible2; % Vector representing transmission from 1 to 2
        cond2 = can_transmit2 & is_susceptible1; % Vector representing transmission from 2 to 1
        
        % Determine if each pair meets either condition
        valid_pairs = cond1 | cond2;
        transmit_pairs = list_ids(valid_pairs,:);
        % toc
        % 
        % disp('Infection in new casual partnerships')
        % tic
        if height(transmit_pairs)>0

            for counter_pairs = 1:height(transmit_pairs)
                id1 = transmit_pairs(counter_pairs,1);
                id2 = transmit_pairs(counter_pairs,2);
                [Pop_hiv,infect_alive,infect_source_data,infect_target_data, newly_infect1, infectionLog] = infect_new_casual(Population,Pop_hiv,pop_index,hiv_index,cfg,infect_alive,tcounter,HIV_st,id1,id2,condom,infect_source_data,infect_target_data, newly_infect1, infectionLog, rec_casual_parts);
            end
            infect_alive.Data=[infect_alive.Data, newly_infect1];
        end

       %toc
    end
else
    warning('Create_casual: could not create a single pair');
end

end