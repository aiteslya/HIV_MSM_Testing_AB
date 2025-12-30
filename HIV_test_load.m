function [Pop_hiv, newly_diagn,prop_diagn_AHI,new_total_diagn_ind,new_AHI_diagn_ind, PrEP_inds, newly_tested, infect_diagn_time, testingLog] = HIV_test_load(Pop_hiv,hiv_index,HIV_st,infect_old,prop_diagn_AHI,tcounter,new_total_diagn_ind,new_AHI_diagn_ind,PrEP_inds, alive_ind, infect_diagn_time, testingLog, Population, rec_casual_parts, pop_index, sim_type, partner_counts)
% PrEP_inds: contains people who were infected while on PrEP
% diagnosis, but not people who were just infected on the same day

    % myname=mfilename;
    year_week = 52;
    newly_diagn = [];% list of indices who was just newly diagnosed
    newly_tested = [];


    inds_can_test = alive_ind.Data(Pop_hiv.Data(hiv_index.status,alive_ind.Data)<=HIV_st.suscPrep);


    % determine who will obtain diagnosis as the result of the test
    hiv_statuses = Pop_hiv.Data(hiv_index.status,inds_can_test);
    % this is a hardcoded condition that people in Fiebig stages 1 and 2
    % cannot be diagnosed
    
    cond1 = (hiv_statuses ~= HIV_st.A23) & (hiv_statuses ~= HIV_st.A1);
    cond2 = hiv_statuses < HIV_st.suscPrep; %

    % the below is a problem
    cond3 = ismember(inds_can_test, infect_old);

    %can_be_diagnosed = cond1 & cond2 & cond3;
    can_be_diagnosed = cond2 & cond3;    

    testing_rates_can_test = Pop_hiv.Data(hiv_index.test_rate, inds_can_test);
    testing_rates_can_test_realized = Pop_hiv.Data(hiv_index.ever_tested, inds_can_test);
    
    hiv_st_can_test = Pop_hiv.Data(hiv_index.status, inds_can_test);
   
    
   if tcounter>2*year_week 
       [scenarioName, pct1, pct2, pct3] = parseSimType(sim_type);
    
        switch scenarioName
            case 'ahi' % increased testing rate in people with early HIV infection
                inds_changed_test = find(hiv_st_can_test<HIV_st.C & hiv_st_can_test>HIV_st.susc);
                
                testing_rates_can_test(inds_changed_test) = max(pct1, testing_rates_can_test(inds_changed_test));
    
                can_be_diagnosed = cond2 & cond3;
            case 'late_inc'    % increase testing rate in people with late HIV infection

                    inds_pot_changed_test = find((tcounter-testing_rates_can_test_realized)/year_week >= pct3); % period of time. NOTE: CHANGED FOR PCT1 TO BE CONSISTENT WITH OTHER INTERVENTIONS
                    if abs(pct2-1.000)<1e-4
                        inds_changed_test = inds_pot_changed_test;
                    else
                        error('HIV_test_load has not been coded for a portion of individuals to increase their testing rate');
                    end
                    new_test_rate = 1/pct3; 
                    testing_rates_can_test(inds_changed_test) = max(new_test_rate, testing_rates_can_test(inds_changed_test));
                
            case 'immigr_late'    % increase testing rate in people with late HIV infection

                % new May 9 update
                inds_pot_changed_test = find((tcounter-testing_rates_can_test_realized)/year_week > pct2); % period of time
                new_test_rate = 1/pct2;
                testing_rates_can_test(inds_pot_changed_test) = max(new_test_rate, testing_rates_can_test(inds_pot_changed_test));
                   
            case 'parts_inc'
    
                % deliniate the population by their status, identifying only these who
                % can get tested
                additional_inds = setdiff(alive_ind.Data, partner_counts(1,:));
                remove_inds = setdiff(partner_counts(1,:), alive_ind.Data);
            
                new_recs = [additional_inds; zeros(1,numel(additional_inds))];
                
            
                partner_counts_updated = [partner_counts'; new_recs']';
            
                % Logical index for columns to keep
                keep_cols = ~ismember(partner_counts_updated(1,:), remove_inds);
                
                % Apply the mask
                partner_counts_updated = partner_counts_updated(:, keep_cols);
            
               % check the correspondes between partner_counts and alive_ind.Data
            
            
               if sum(abs(partner_counts_updated(1, :)-alive_ind.Data))>0
                   % Extract current indices from first row of partner_counts_updated
                    current_order = partner_counts_updated(1, :);
                    
                    % Get the order that matches alive_ind.Data
                    [~, reorder_idx] = ismember(alive_ind.Data, current_order);
                    
                    % Safety check: if any index is not found, that's an error
                    if any(reorder_idx == 0)
                        error('Some alive individuals are missing in partner_counts_updated.');
                    end
                
                    % Reorder the columns
                    partner_counts_updated = partner_counts_updated(:, reorder_idx);
               end
    
                % Extract indices (row 1) and values (row 2) from partner_counts_updated
                all_inds = partner_counts_updated(1, :);
                all_counts = partner_counts_updated(2, :);
                
                % Find positions of inds_can_test in all_inds
                [~, pos] = ismember(inds_can_test, all_inds);
                
                % Safety check
                if any(pos == 0)
                    warning('Some inds_can_test are missing from partner_counts_updated — assigning 0 to those.');
                end
                
                % Get partner counts, defaulting to 0 for missing
                partner_counts_can_test = zeros(size(inds_can_test));
                partner_counts_can_test(pos > 0) = all_counts(pos(pos > 0));
    
                % Old code
                % inds_changed_test = find(partner_counts_can_test >= pct1);
                % end of old code
                
                % new June 17 update
                dur = 1/pct2;
                inds_pot_changed_test = find((tcounter-testing_rates_can_test_realized)/year_week > dur & partner_counts_can_test >= pct1); % period of time

                inds_changed_test = inds_pot_changed_test;

    
                testing_rates_can_test(inds_changed_test) = max(pct2, testing_rates_can_test(inds_changed_test));
            case 'immigr_parts'
    
                % deliniate the population by their status, identifying only these who
                % can get tested
                additional_inds = setdiff(alive_ind.Data, partner_counts(1,:));
                remove_inds = setdiff(partner_counts(1,:), alive_ind.Data);
            
                new_recs = [additional_inds; zeros(1,numel(additional_inds))];
                
            
                partner_counts_updated = [partner_counts'; new_recs']';
            
                % Logical index for columns to keep
                keep_cols = ~ismember(partner_counts_updated(1,:), remove_inds);
                
                % Apply the mask
                partner_counts_updated = partner_counts_updated(:, keep_cols);
            
               % check the correspondes between partner_counts and alive_ind.Data
            
            
               if sum(abs(partner_counts_updated(1, :)-alive_ind.Data))>0
                   % Extract current indices from first row of partner_counts_updated
                    current_order = partner_counts_updated(1, :);
                    
                    % Get the order that matches alive_ind.Data
                    [~, reorder_idx] = ismember(alive_ind.Data, current_order);
                    
                    % Safety check: if any index is not found, that's an error
                    if any(reorder_idx == 0)
                        error('Some alive individuals are missing in partner_counts_updated.');
                    end
                
                    % Reorder the columns
                    partner_counts_updated = partner_counts_updated(:, reorder_idx);
               end
    
                % Extract indices (row 1) and values (row 2) from partner_counts_updated
                all_inds = partner_counts_updated(1, :);
                all_counts = partner_counts_updated(2, :);
                
                % Find positions of inds_can_test in all_inds
                [~, pos] = ismember(inds_can_test, all_inds);
                
                % Safety check
                if any(pos == 0)
                    warning('Some inds_can_test are missing from partner_counts_updated — assigning 0 to those.');
                end
                
                % Get partner counts, defaulting to 0 for missing
                partner_counts_can_test = zeros(size(inds_can_test));
                partner_counts_can_test(pos > 0) = all_counts(pos(pos > 0));
    
                % Old code
                % inds_changed_test = find(partner_counts_can_test >= pct1);
                % end of old code
                
                % new June 17 update
                dur = 1/pct3;
                inds_pot_changed_test = find((tcounter-testing_rates_can_test_realized)/year_week > dur & partner_counts_can_test >= pct2); % period of time

                inds_changed_test = inds_pot_changed_test;

    
                testing_rates_can_test(inds_changed_test) = max(pct3, testing_rates_can_test(inds_changed_test));
        end
   end
    % determine who will be tested
    r_arr = rand(1, numel(inds_can_test));
    tested = r_arr<(testing_rates_can_test/year_week);

    
    ind_can_be_diagnosed =  inds_can_test(can_be_diagnosed);

    hiv_statuses_can_be_diagnosed = Pop_hiv.Data(hiv_index.status, ind_can_be_diagnosed);


   assert((sum(hiv_statuses_can_be_diagnosed<1)+sum(hiv_statuses_can_be_diagnosed>7))==0);

   if sum(tested)>0
       %sum(tested)
        newly_tested = inds_can_test(tested);

        for ind = inds_can_test(tested)
            id = Population.Data(pop_index.id,ind);
            HIV_st_id = Pop_hiv.Data(hiv_index.status, ind);
            num_parts_id = numel(retrievePartners(rec_casual_parts, id));
            record = {id ind tcounter HIV_st_id, Population.Data(pop_index.age,ind) ...
                                Population.Data(pop_index.casual_prop,ind) ...
                                tcounter - Pop_hiv.Data(hiv_index.date_inf,ind) num_parts_id }; 
        testingLog.logEvent(record{:});
        end
    end

    was_diagnosed_log_index = tested & can_be_diagnosed;

    was_diagnosed_inds = inds_can_test(was_diagnosed_log_index);


    assert(sum(Pop_hiv.Data(hiv_index.status, was_diagnosed_inds)<0 & Pop_hiv.Data(hiv_index.status, was_diagnosed_inds)>=7)==0);

    was_diagnosed_PrEP_inds = PrEP_inds.Data(ismember(PrEP_inds.Data,was_diagnosed_inds));

    PrEP_inds.Data = setdiff(PrEP_inds.Data,was_diagnosed_PrEP_inds);

    % record that people have been tested at least once
    Pop_hiv.Data(hiv_index.ever_tested, inds_can_test(tested)) = tcounter;

    % for individuals who were diagnosed update: 1. HIV status in the
    % cascade of care 2. time of detection
    Pop_hiv.Data(hiv_index.status, was_diagnosed_inds) = Pop_hiv.Data(hiv_index.status, was_diagnosed_inds) + 7;
    Pop_hiv.Data(hiv_index.date_diagn, was_diagnosed_inds) = tcounter;

    % 

    
    
    new_total_diagn_ind=[new_total_diagn_ind, was_diagnosed_inds];

    % record the age of infection at diagnosis

    infect_diagn_time.Data = [infect_diagn_time.Data tcounter - Pop_hiv.Data(hiv_index.date_inf, new_total_diagn_ind)];

    newly_diagn=[newly_diagn, was_diagnosed_inds];

end