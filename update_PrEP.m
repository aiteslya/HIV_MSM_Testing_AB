function [Pop_hiv, inds_changed_PreP, PrEP_inds] = update_PrEP(Pop_hiv, Population, alive_ind, infect_alive, pop_index, hiv_index, HIV_st, ind_list, old_prop, cfg_prep, test_edge, PrEP_inds, test_distr, rec_casual_parts, test_rates, year_counter)
% this function checks new propensities of individuals who use PrEP. If the
% propensity is lower than it used to be, then the function attempts to
% find someone with similar or larger propensity and re-assign the PrEP

% a version of this function will need to be created for the forward
% scenaro, in order to set the correct testing rate

% input verification dimensions of ind_list and old props should be 1 x n
% and should be the same

% Input verification: Check if both inputs are vectors and have the same length
if ~isvector(ind_list) || ~isvector(old_prop) || numel(ind_list) ~= numel(old_prop)
    error('Inputs ind_list and old_props must be vectors of the same length.');
end

isInfected = ismember(alive_ind.Data, infect_alive.Data);

% Select elements from alive_ind.Data that are not in infect_alive.Data
ind_non_infect = alive_ind.Data(~isInfected);
% select entries where individual does not take PrEP
ind_non_infect = ind_non_infect(Pop_hiv.Data(hiv_index.status, ind_non_infect) == HIV_st.susc);

all_propensities = Population.Data(pop_index.casual_prop, ind_non_infect);
% new (updated) propensities of individuals who take PrEP
new_prop = Population.Data(pop_index.casual_prop, ind_list);

indices = ind_list(new_prop<old_prop);
current_props = new_prop(new_prop<old_prop);

inds_changed_PreP = [];

if numel(indices) > 0
    counter = 1;
    for ind = indices
        max_prop = max(all_propensities(all_propensities>current_props(counter)));

        if numel(max_prop)>0
            ind_repl = ind_non_infect(find(all_propensities == max_prop, 1));
            if max_prop > current_props(counter)
                % switch out PrEP status
                if Pop_hiv.Data(hiv_index.status,ind) == HIV_st.suscPrep % has not been infected while on PrEP
                    Pop_hiv.Data(hiv_index.status,ind) = HIV_st.susc;
                else
                    disp('Someone was infected while on PrEP');
                end
                Pop_hiv.Data(hiv_index.status, ind_repl) = HIV_st.suscPrep;
                ind_non_infect(:, find(all_propensities == max_prop, 1)) = [];
                all_propensities(:, find(all_propensities == max_prop, 1)) = [];

                % assign testing rates
                % new PrEP user
                if year_counter < 2025
                    Pop_hiv.Data(hiv_index.test_rate, ind_repl) = cfg_prep.old_rate;
                else
                    Pop_hiv.Data(hiv_index.test_rate, ind_repl) = cfg_prep.new_rate;
                end

                Pop_hiv.Data(hiv_index.testing_bin, ind_repl) = 5; % hard-coded: PrEP testing bin

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % de-enrolled user
                
                id = Population.Data(pop_index.id, ind);
                

                % retrieve number of partners
               
                if isKey(rec_casual_parts, id)
                    num_parts = numel(rec_casual_parts(id));
                    num_parts_bin = discretize(num_parts,test_edge);
                else
                    num_parts_bin = 1;
                end
                Pop_hiv.Data(hiv_index.test_rate, ind) = randsample(test_rates,1, true,test_distr(:, num_parts_bin));
                Pop_hiv.Data(hiv_index.testing_bin, ind) = num_parts_bin;

                %%%%%%


                inds_changed_PreP = [inds_changed_PreP; ind_repl ind];

                % update PrEP inds
                PrEP_inds.Data = setdiff(PrEP_inds.Data, ind);
                PrEP_inds.Data = [PrEP_inds.Data, ind_repl];
                
            end
        end
        counter = counter + 1;
    end
end
end