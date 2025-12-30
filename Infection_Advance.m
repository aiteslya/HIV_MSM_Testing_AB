function [Pop_hiv]=Infection_Advance(Pop_hiv,hiv_index,cfg,infect_alive,HIV_st,tcounter)
% this function goes over all individuals who are infected and determines whether they advance to the next stage

% convert probabilities from per year to per week

year_week = 52; % duration of year in weeks
gamma_a23W = cfg.gamma_a23/year_week;
gamma_a45W = cfg.gamma_a45/year_week;
gamma_cW = cfg.gamma_c/year_week;
gamma_eaW = cfg.gamma_ea/year_week;


for ind=infect_alive.Data
    st=Pop_hiv.Data(hiv_index.status,ind);
    if ismember(st,[HIV_st.A1,HIV_st.A23,HIV_st.A45,HIV_st.C,HIV_st.EA,HIV_st.A1D,HIV_st.A23D,HIV_st.A45D,HIV_st.CD,HIV_st.EAD,HIV_st.A1T,HIV_st.A23T,HIV_st.A45T,HIV_st.CT,HIV_st.EAT])
        % new version, where durations of infections were sampled from distributions  
        if Pop_hiv.Data(hiv_index.date_new_stage,ind)<=tcounter % time to move to the next stage and sample the duration of that stage
            if st==HIV_st.A1 | st==HIV_st.A1D | st==HIV_st.A1T
                Gamma = gamma_a23W;
                % update status
                Pop_hiv.Data(hiv_index.status,ind)=Pop_hiv.Data(hiv_index.status,ind)+1;
                % sample form the exponential distribution
                dur=ceil(-log(1-rand)/Gamma);
                Pop_hiv.Data(hiv_index.date_new_stage,ind) = Pop_hiv.Data(hiv_index.date_new_stage,ind)+dur;
            elseif st==HIV_st.A23 | st==HIV_st.A23D | st==HIV_st.A23T
                Gamma = gamma_a45W;
                % update status
                Pop_hiv.Data(hiv_index.status,ind)=Pop_hiv.Data(hiv_index.status,ind)+1;
                % sample form the exponential distribution
                dur=ceil(-log(1-rand)/Gamma);
                Pop_hiv.Data(hiv_index.date_new_stage,ind)=Pop_hiv.Data(hiv_index.date_new_stage,ind)+dur;
            elseif st==HIV_st.A45 | st==HIV_st.A45D | st==HIV_st.A45T
                Gamma = gamma_cW;
                mean_dur=1/Gamma;
                dur=ceil(chronic_dur_sample(mean_dur,'erlang',80));
                % update status
                Pop_hiv.Data(hiv_index.status,ind)=Pop_hiv.Data(hiv_index.status,ind)+1;
                Pop_hiv.Data(hiv_index.date_new_stage,ind)=Pop_hiv.Data(hiv_index.date_new_stage,ind)+dur;
            elseif st==HIV_st.C | st==HIV_st.CD | st==HIV_st.CT
                Gamma = gamma_eaW;
                % update status
                Pop_hiv.Data(hiv_index.status,ind)=Pop_hiv.Data(hiv_index.status,ind)+1;
                % sample form the exponential distribution
                dur=ceil(-log(1-rand)/Gamma);
                Pop_hiv.Data(hiv_index.date_new_stage,ind)=Pop_hiv.Data(hiv_index.date_new_stage,ind)+dur;
             elseif st==HIV_st.EA | st==HIV_st.EAD | st==HIV_st.EAT
                % update status
                Pop_hiv.Data(hiv_index.status,ind)=Pop_hiv.Data(hiv_index.status,ind)+1;
                Pop_hiv.Data(hiv_index.date_new_stage,ind)=cfg.T*year_week +10;
            end
        end     
    end
end

end