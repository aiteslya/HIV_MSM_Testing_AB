function [Pop_hiv,infect_alive,newly_infect,infect_source_data,infect_target_data, infectionLog] = infection(Population,Pop_hiv,pop_index,hiv_index,cfg,infect_alive,tcounter,HIV_st,condom,new_rels_casual,infect_source_data,infect_target_data, all_casual_parts, Rels_steady, infectionLog, rec_casual_parts)

% this function executes infection events per day
% set up the parameters by converting probabilities per year to probability
% per day
% condom matrix: 4 by 6 : 
% 1-st row - steady no PrEP, 2-nd - steady with PrEP, 3-rd casual no PrEP,
% 4-th - casual with PrEP
% each column denotes age bracket: 1: 15-25, 2: 25-35, 3: 35-45, 4: 45-55,
% 5: 55-65, 1: 65-75
% thus each entry of condom matrix denotes probability of using condom
% during AI
    year_week=52; % number of weeks in a year
    year_day = 365;
    week = 7;
    % probability of anal intercourse in a casual partnership per week
    cont_prob_st=cfg.cont_rate_st/year_week;
    cont_prob_cas=cfg.cont_rate_cas/year_week; 
    newly_infect=[];

    for ind = infect_alive.Data
        if ismember(Pop_hiv.Data(hiv_index.status,ind),[HIV_st.A1,HIV_st.A23,HIV_st.A45,HIV_st.C,HIV_st.EA,HIV_st.LA,HIV_st.A1D,HIV_st.A23D,HIV_st.A45D,HIV_st.CD,HIV_st.EAD,HIV_st.LAD,HIV_st.A1T,HIV_st.A23T,HIV_st.A45T,HIV_st.CT,HIV_st.EAT,HIV_st.LAT]) % individual can transmit
            % determine the status
            st = Pop_hiv.Data(hiv_index.status,ind);
            
            % id of the individual
            id = Pop_hiv.Data(hiv_index.id,ind);

            % determine the multiplier of probability of transmission per
            % contact depending on the HIV stage of the index individual

            if st==HIV_st.A1 | st==HIV_st.A1D | st==HIV_st.A1T
                multiplier=cfg.a1;
            elseif st==HIV_st.A23 | st==HIV_st.A23D | st==HIV_st.A23T
                multiplier=cfg.a2;
            elseif st==HIV_st.A45 | st==HIV_st.A45D | st==HIV_st.A45T
               multiplier=cfg.a3;
            elseif st==HIV_st.C | st==HIV_st.CD | st==HIV_st.CT
                multiplier=1;
            elseif st==HIV_st.EA | st==HIV_st.EAD | st==HIV_st.EAT
                multiplier=cfg.a4;
            elseif st==HIV_st.LA | st==HIV_st.LAD | st==HIV_st.LAT
                multiplier=0;
            end
            
            CPStMult=1;
            CPCasMult=1;
            prob_trans_chr=cfg.prob_trans_chr;
            % reduction in contact rate due to being diagnosed
            if ismember(st,[HIV_st.A1D,HIV_st.A23D,HIV_st.A45D,HIV_st.CD,HIV_st.EAD,HIV_st.LAD,HIV_st.A1T,HIV_st.A23T,HIV_st.A45T,HIV_st.CT,HIV_st.EAT,HIV_st.LAT])
                if (tcounter - Pop_hiv.Data(hiv_index.date_diagn,ind))<year_week
                    CPStMult=CPStMult*cfg.diagn_cont_mult;
                    CPCasMult=CPCasMult*cfg.diagn_cont_mult;
                end
            end

            % transmission within steady partnerships
            % retrieve partners
            % ids
            
            sp1_id = Population.Data(pop_index.sp1_id,ind);

            % indices
            if sp1_id>0
                % determine which fraction of a week this partnership will last
                ind_rel = find(ismember(Rels_steady.Data(:,1:2),sort([id, sp1_id]),'rows'));
                time_end = Rels_steady.Data(ind_rel,4);
                if time_end<tcounter+1
                    dur = time_end - tcounter;
                else
                    dur = 1;
                end
                ind_p1=find(Population.Data(pop_index.id,:)==sp1_id);
                HIVst_p1=Pop_hiv.Data(hiv_index.status,ind_p1);

                if HIVst_p1==HIV_st.susc | HIVst_p1==HIV_st.suscPrep % susceptible to infection

                    age_bin1=Population.Data(pop_index.age_bin,ind);
                    age_bin2=Population.Data(pop_index.age_bin,ind_p1);
                
                    cond1=condom(1,age_bin1);
                    cond2=condom(1,age_bin2);

                    if HIVst_p1==HIV_st.susc
                        prob_infect=CPStMult*cont_prob_st*prob_trans_chr*multiplier;
                    else
                        prob_infect=CPStMult*cont_prob_st*prob_trans_chr*multiplier*cfg.prep_red;
                    end

                    % reconcile condom use
                    % probability that condom is not used
                    cond_effect=1-max(cond1,cond2);

                    prob_infect=prob_infect*cond_effect*dur;
                    
                    r1=rand;
                    if r1<prob_infect % infection event
                        newly_infect=[newly_infect,ind_p1];

                        % record the state of newly infected individual
                        % prior to infection
                        HIV_st_p1=Pop_hiv.Data(hiv_index.status,ind_p1);
                        Pop_hiv.Data(hiv_index.status,ind_p1)=1;

                        Pop_hiv.Data(hiv_index.date_inf,ind_p1)=tcounter;
                        % new addition: estimate the time until progression
                        % to I_{p,2+3}
                        % convert from years to days
                        infect_progres_rate=cfg.gamma_a1/year_day;
                        % sample form the exponential distribution
                        dur=ceil(-log(1-rand)/infect_progres_rate);
                        Pop_hiv.Data(hiv_index.date_new_stage,ind_p1)=tcounter+dur/week;
            
                        % record who infected who
                        id=Population.Data(pop_index.id,ind);
                        age=Population.Data(pop_index.age,ind);
                        prop=Population.Data(pop_index.casual_prop,ind);
                        HIV_st_id=Pop_hiv.Data(hiv_index.status,ind);
                        age_infect = Pop_hiv.Data(hiv_index.date_inf,ind);
        
                        id_p1=Population.Data(pop_index.id,ind_p1);
                        age_p1=Population.Data(pop_index.age,ind_p1);
                        prop_p1=Population.Data(pop_index.casual_prop,ind_p1);
                        
        
                        if tcounter>cfg.T_burn*year_week % HARDCODED
                            % old code, update by more efficient and
                            % expansive data record below
                            % % record the details for the source of infection
                            % infect_source_data=[infect_source_data;tcounter id age prop HIV_st_id age_infect];
                            % % record the details for the target of infection
                            % infect_target_data=[infect_target_data; tcounter id_p1 age_p1 prop_p1 HIV_st_id HIV_st_p1];

                            % retrieve recent number of casual partners
                            num_parts_id = numel(retrievePartners(rec_casual_parts, id));
                            num_parts_p1 = numel(retrievePartners(rec_casual_parts, id_p1));

                            record = {id id_p1 ind ind_p1 tcounter HIV_st_id HIV_st_p1...
                                Population.Data(pop_index.age,ind) Population.Data(pop_index.age,ind_p1)...
                                Population.Data(pop_index.casual_prop,ind) Population.Data(pop_index.casual_prop,ind_p1)...
                                tcounter - Pop_hiv.Data(hiv_index.date_inf,ind) num_parts_id num_parts_p1 1 Pop_hiv.Data(hiv_index.test_rate,ind) Pop_hiv.Data(hiv_index.test_rate,ind_p1)}; % 1 denotes the type of partnership, 1:steady, 2: non-steady
                            infectionLog.logEvent(record{:});
                        end


                    end
                end
            end

            % transmission within casual partnerships
            % retrieve casual partners
%             cp_id_arr=Pop_casual.Data(2:2:2*max_num_casual,ind);
%             cp_id_arr=cp_id_arr(cp_id_arr>0);

            % id=Population.Data(pop_index.id,ind);
  
            parts = retrievePartners(all_casual_parts, id);
            
            
            if ~isempty(parts) % there are casual partners
                cp_id_arr = zeros(numel(parts), 1);
                cp_start = zeros(numel(parts),1);
                cp_dur = zeros(numel(parts), 1);
                
                counter = 1;
                for part = parts
                    cp_id_arr(counter) = part.id;
                    cp_start(counter) = part.start_date;
                    cp_dur(counter) = part.end_date;
                    
                    counter = counter + 1;
                end

                cp_counter = 1;
                for cp_id=cp_id_arr'
                    % get id of the ego
                    if ~ismember(sort([id,cp_id]),new_rels_casual(:,1:2),'rows') % these individuals did not just get together
                        % note: the check for whether the individuals just
                        % got together was removed in a version where step
                        % is a week
                        age_bin1=Population.Data(pop_index.age_bin,ind);
                        cond1=condom(2,age_bin1);
                        
                        
                        ind_cp=find(Population.Data(pop_index.id,:)==cp_id);
                        age_bin2=Population.Data(pop_index.age_bin,ind_cp);
                        
                        HIVst_cp=Pop_hiv.Data(hiv_index.status,ind_cp);
                        
                        if HIVst_cp==HIV_st.susc | HIVst_cp==HIV_st.suscPrep % susceptible to infection
                            cond2=condom(2,age_bin2);
                            if HIVst_cp==HIV_st.susc
                                prob_infect=CPCasMult*cont_prob_cas*prob_trans_chr*multiplier;
                                
                            else
                                prob_infect=CPCasMult*cont_prob_cas*prob_trans_chr*multiplier*cfg.prep_red;
                              
                            end

                            time_end = cp_dur(cp_counter);
                            if time_end<tcounter+1
                                dur = time_end - tcounter;
                            else
                                dur = 1;
                            end

                            
                            % reconcile condom use
                            % probability that condom is not used
                            cond_effect=1-max(cond1,cond2); % probability of condom not used
                            prob_infect=prob_infect*cond_effect*dur;
                            r1=rand;

                            if r1<prob_infect % infection event
                                newly_infect=[newly_infect,ind_cp];
                                % record the state before changing it
                                HIV_st_cp=Pop_hiv.Data(hiv_index.status,ind_cp);
                                Pop_hiv.Data(hiv_index.status,ind_cp)=1;
                                Pop_hiv.Data(hiv_index.date_inf,ind_cp)=tcounter;

                                % new addition: estimate the time until progression
                                % to I_{p,2+3}
                                % convert from years to days
                                infect_progres_rate=cfg.gamma_a1/year_day;
                                % sample form the exponential distribution
                                dur=ceil(-log(1-rand)/infect_progres_rate);
                                Pop_hiv.Data(hiv_index.date_new_stage,ind_cp)=tcounter+dur/week;

                                age=Population.Data(pop_index.age,ind);
                                prop=Population.Data(pop_index.casual_prop,ind);
                                HIV_st_id=Pop_hiv.Data(hiv_index.status,ind);
                                age_infect = Pop_hiv.Data(hiv_index.date_inf,ind);
                
                                id_cp=Population.Data(pop_index.id,ind_cp);
                                age_cp=Population.Data(pop_index.age,ind_cp);
                                prop_cp=Population.Data(pop_index.casual_prop,ind_cp);
                                
                
                                % record the details for the source of
                                % infection: time point, id, biological age,
                                % propensity to acquire partners, status,
                                % age of infection of infector
                                infect_source_data=[infect_source_data; tcounter id age prop HIV_st_id age_infect];
                                % record the details for the target of infection: time point, id, biological age of infector,
                                % propensity to acquire partners of
                                % infectee, status and age of infection of infector
                                infect_target_data=[infect_target_data; tcounter id_cp age_cp prop_cp HIV_st_id HIV_st_cp];

                                % new record keeping code

                                % retrieve recent number of casual partners
                                num_parts_id = numel(retrievePartners(rec_casual_parts, id));
                                num_parts_cp = numel(retrievePartners(rec_casual_parts, id_cp));
    
                                record = {id id_cp ind ind_cp tcounter HIV_st_id HIV_st_cp...
                                    Population.Data(pop_index.age,ind) Population.Data(pop_index.age,ind_cp)...
                                    Population.Data(pop_index.casual_prop,ind) Population.Data(pop_index.casual_prop,ind_cp)...
                                    tcounter - Pop_hiv.Data(hiv_index.date_inf,ind) num_parts_id num_parts_cp 2 Pop_hiv.Data(hiv_index.test_rate,ind) Pop_hiv.Data(hiv_index.test_rate,ind_cp)};
                                infectionLog.logEvent(record{:});
                            end

                        end
                    % else
                    %     disp('This pair was created today');
                    end
                    cp_counter = cp_counter +1;
                end
            end
        end
    end
    infect_alive.Data = [infect_alive.Data , newly_infect];
end