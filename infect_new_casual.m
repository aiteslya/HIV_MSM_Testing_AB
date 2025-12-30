function [Pop_hiv,infect_alive,infect_source_data,infect_target_data, newly_infect, infectionLog] = infect_new_casual(Population,Pop_hiv,pop_index,hiv_index,cfg,infect_alive,tcounter,HIV_st,id1,id2,condom,infect_source_data,infect_target_data, newly_infect, infectionLog, rec_casual_parts)

% this function simulates event of infection that may occur upon the first
% sex act of individuals who reside in the population table under indices
% ind1 and ind2
% set up the parameters by converting probabilities per year to probability
% per day
% condom matrix: 4 by 6 : 
% 1-st row - steady no PrEP, 2-nd - steady with PrEP, 3-rd casual no PrEP,
% 4-th - casual with PrEP
% each column denotes age bracket: 1: 15-25, 2: 25-35, 3: 35-45, 4: 45-55,
% 5: 55-65, 1: 65-75
% thus each entry of condom matrix denotes probability of using condom
% during AI
    ind1 = find(Population.Data(pop_index.id,:)==id1);
    ind2 = find(Population.Data(pop_index.id,:)==id2);
    year_week = 52;
    
    % can individual residing at ind1 transmit?
    st1=ismember(Pop_hiv.Data(hiv_index.status,ind1),[HIV_st.A1,HIV_st.A23,HIV_st.A45,HIV_st.C,HIV_st.EA,HIV_st.LA,HIV_st.A1D,HIV_st.A23D,HIV_st.A45D,HIV_st.CD,HIV_st.EAD,HIV_st.LAD,HIV_st.A1T,HIV_st.A23T,HIV_st.A45T,HIV_st.CT,HIV_st.EAT,HIV_st.LAT]);
    % can individual residing at ind2 transmit?
    st2=ismember(Pop_hiv.Data(hiv_index.status,ind2),[HIV_st.A1,HIV_st.A23,HIV_st.A45,HIV_st.C,HIV_st.EA,HIV_st.LA,HIV_st.A1D,HIV_st.A23D,HIV_st.A45D,HIV_st.CD,HIV_st.EAD,HIV_st.LAD,HIV_st.A1T,HIV_st.A23T,HIV_st.A45T,HIV_st.CT,HIV_st.EAT,HIV_st.LAT]);
    
    if (st1+st2)==1 % exactly one person can transmit
        % find the index of infectious
        ind=ind1*st1+ind2*(1-st1);
        % set the index of the partner of infectious
        ind_cp=ind1*(1-st1)+ind2*st1;
        
        % determine the status
        st=Pop_hiv.Data(hiv_index.status,ind);
        % determine the multiplier of probability of transmission per
        % contact

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
        
        CPCasMult=1;
        prob_trans_chr=cfg.prob_trans_chr;
        if ismember(st,[HIV_st.A1D,HIV_st.A23D,HIV_st.A45D,HIV_st.CD,HIV_st.EAD,HIV_st.LAD,HIV_st.A1T,HIV_st.A23T,HIV_st.A45T,HIV_st.CT,HIV_st.EAT,HIV_st.LAT])
            if (tcounter - Pop_hiv.Data(hiv_index.date_diagn, ind))<year_week
                CPCasMult=CPCasMult*cfg.diagn_cont_mult;
            end
        end

        % transmission within casual partnerships
        HIVst_cp=Pop_hiv.Data(hiv_index.status,ind_cp);
        
        age_bin1 = Population.Data(pop_index.age_bin,ind);
        age_bin2 = Population.Data(pop_index.age_bin,ind_cp);
        cond1=condom(2,age_bin1);
        if HIVst_cp==HIV_st.susc | HIVst_cp==HIV_st.suscPrep % susceptible to infection
            cond2=condom(2,age_bin2);
            if HIVst_cp==HIV_st.susc
                prob_infect=CPCasMult*prob_trans_chr*multiplier;
            elseif HIVst_cp==HIV_st.suscPrep
                prob_infect=CPCasMult*prob_trans_chr*multiplier*cfg.prep_red;
            end
            
            % reconcile condom use
            % probability that condom is not used
            cond_effect=1-max(cond1,cond2);
            prob_infect=prob_infect*cond_effect;
            r1=rand;
            if r1<prob_infect % infection event
                newly_infect=[newly_infect,ind_cp];
                HIV_st_cp=Pop_hiv.Data(hiv_index.status,ind_cp);
                Pop_hiv.Data(hiv_index.status,ind_cp)=1;
                Pop_hiv.Data(hiv_index.date_inf,ind_cp)=tcounter;

                % new addition: estimate the time until progression
                % to I_{p,2+3}
                % convert from years to days
                infect_progres_rate=cfg.gamma_a1/year_week;
                % sample form the exponential distribution
                dur=ceil(-log(1-rand)/infect_progres_rate);
                Pop_hiv.Data(hiv_index.date_new_stage,ind_cp)=tcounter+dur;

                s_fl=1;
                % record the data about who infected whom
                id=Population.Data(pop_index.id,ind);
                age=Population.Data(pop_index.age,ind);
                prop=Population.Data(pop_index.casual_prop,ind);
                HIV_st=Pop_hiv.Data(hiv_index.status,ind);
                age_infect = Pop_hiv.Data(hiv_index.date_inf,ind);

                id_cp=Population.Data(pop_index.id,ind_cp);
                age_cp=Population.Data(pop_index.age,ind_cp);
                prop_cp=Population.Data(pop_index.casual_prop,ind_cp);
                                
                if tcounter>=cfg.T_burn*year_week
                    % % record the details for the source of infection
                    % infect_source_data=[infect_source_data;tcounter id age prop HIV_st age_infect];
                    % % record the details for the target of infection
                    % infect_target_data=[infect_target_data; tcounter id_cp age_cp prop_cp HIV_st HIV_st_cp];
                    % retrieve recent number of casual partners
                    num_parts_id = numel(retrievePartners(rec_casual_parts, id));
                    num_parts_cp = numel(retrievePartners(rec_casual_parts, id_cp));

                    record = {id id_cp ind ind_cp tcounter st HIV_st_cp...
                        Population.Data(pop_index.age,ind) Population.Data(pop_index.age,ind_cp)...
                        Population.Data(pop_index.casual_prop,ind) Population.Data(pop_index.casual_prop,ind_cp)...
                        tcounter - Pop_hiv.Data(hiv_index.date_inf,ind) num_parts_id num_parts_cp 2 Pop_hiv.Data(hiv_index.test_rate,ind) Pop_hiv.Data(hiv_index.test_rate,ind_cp)};
                    infectionLog.logEvent(record{:});
                end
                
            end
        end

        
    end
end