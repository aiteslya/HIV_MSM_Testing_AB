function [Pop_hiv, newly_diagn,prop_diagn_AHI,new_total_diagn_ind,new_AHI_diagn_ind, PrEP_inds, newly_tested, infect_diagn_time, testingLog] = HIV_test(Pop_hiv,hiv_index,HIV_st,infect_old,prop_diagn_AHI,tcounter,new_total_diagn_ind,new_AHI_diagn_ind,PrEP_inds, alive_ind, infect_diagn_time, testingLog, Population, rec_casual_parts, pop_index)
% PrEP_inds: contains people who were infected while on PrEP
% diagnosis, but not people who were just infected on the same day

    % myname=mfilename;
    year_week = 52;
    newly_diagn = [];% list of indices who was just newly diagnosed
    newly_tested = [];

    % deliniate the population by their status, identifying only these who
    % can get tested
    inds_can_test = alive_ind.Data(Pop_hiv.Data(hiv_index.status,alive_ind.Data)<=HIV_st.suscPrep);

    % determine who will be tested
    r_arr = rand(1, numel(inds_can_test));
    tested = r_arr<(Pop_hiv.Data(hiv_index.test_rate, inds_can_test)/year_week);

    % determine who will obtain diagnosis as the result of the test
    hiv_statuses = Pop_hiv.Data(hiv_index.status,inds_can_test);
    % this is a hardcoded condition that people in Fiebig stages 1 and 2
    % cannot be diagnosed
    
    cond1 = (hiv_statuses ~= HIV_st.A23) & (hiv_statuses ~= HIV_st.A1);
    cond2 = hiv_statuses < HIV_st.suscPrep; %

    % the below is a problem
    cond3 = ismember(inds_can_test, infect_old);

    can_be_diagnosed = cond1 & cond2 & cond3;

    ind_can_be_diagnosed =  inds_can_test(can_be_diagnosed);

    hiv_statuses_can_be_diagnosed = Pop_hiv.Data(hiv_index.status, ind_can_be_diagnosed);

   assert((sum(hiv_statuses_can_be_diagnosed<3)+sum(hiv_statuses_can_be_diagnosed>7))==0);

   if sum(tested)>0
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