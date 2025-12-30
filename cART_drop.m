function [Pop_hiv,total_diagn_ind,total_sup_ind]=cART_drop(Pop_hiv,hiv_index,HIV_st,cfg,total_diagn_ind,total_sup_ind,tcounter)
    % individuals drop out of cART treatment
    
    switch ModelConstants.unit
        case 'week'
            year = 52;
    end

    myname=mfilename;
    newly_drop=[];
    xi = 1 -  exp(-cfg.xi/year);
    
    % define mean sojourn times
    gamma_a1W = 1-exp(-cfg.gamma_a1/year);
    gamma_a23W = 1-exp(-cfg.gamma_a23/year);
    gamma_a45W = 1-exp(-cfg.gamma_a45/year);
    gamma_cW = 1-exp(-cfg.gamma_c/year);
    gamma_eaW = 1-exp(-cfg.gamma_ea/year);

    for ind=total_sup_ind
        r1=rand;
        if r1<xi % drop out of treatment
            st=Pop_hiv.Data(hiv_index.status,ind);
            Pop_hiv.Data(hiv_index.date_diagn,ind) = tcounter;
            Pop_hiv.Data(hiv_index.date_art,ind)=-40000;
            Pop_hiv.Data(hiv_index.date_sup,ind)=-40000;
            newly_drop=[newly_drop,ind];
            if st==HIV_st.A1S
                Pop_hiv.Data(hiv_index.status,ind)=HIV_st.A1D;
                
                Gamma=gamma_a1W;

               % sample form the exponential distribution
                r1 = rand;
                dur=ceil(-log(1-r1)/Gamma);
                Pop_hiv.Data(hiv_index.date_new_stage,ind)=tcounter+dur;

            elseif st==HIV_st.A23S
                Pop_hiv.Data(hiv_index.status,ind)=HIV_st.A23D;
                Gamma=gamma_a23W;
                % sample form the exponential distribution
                r1 = rand;
                dur=ceil(-log(1-r1)/Gamma);
                Pop_hiv.Data(hiv_index.date_new_stage,ind)=tcounter+dur;
            elseif st==HIV_st.A45S
                Pop_hiv.Data(hiv_index.status,ind)=HIV_st.A45D;
                
                Gamma=gamma_a45W;
                % sample form the exponential distribution
                r1 = rand;
                dur=ceil(-log(1-r1)/Gamma);
                Pop_hiv.Data(hiv_index.date_new_stage,ind)=tcounter+dur;

            elseif st==HIV_st.CS
                Pop_hiv.Data(hiv_index.status,ind)=HIV_st.CD;
                
                Gamma=gamma_cW;
                mean_dur=1/Gamma;
                dur=ceil(chronic_dur_sample(mean_dur,'erlang',80));
                Pop_hiv.Data(hiv_index.date_new_stage,ind)=tcounter+dur;
            elseif st==HIV_st.EAS
                Pop_hiv.Data(hiv_index.status,ind)=HIV_st.EAD;

                Gamma=gamma_eaW;
                % sample form the exponential distribution
                r1 = rand;
                dur=ceil(-log(1-r1)/Gamma);
                Pop_hiv.Data(hiv_index.date_new_stage,ind)=tcounter+dur;
            elseif  st==HIV_st.LAS
                Pop_hiv.Data(hiv_index.status,ind)=HIV_st.LAD;
                Pop_hiv.Data(hiv_index.date_new_stage,ind)=cfg.T*year +10;
            else
                error([myname,': someone who should not be be able to drop from treatment is attempting to be dropped']);
            end
        end
    end

    total_diagn_ind = [total_diagn_ind,newly_drop];
    total_sup_ind = setdiff(total_sup_ind,newly_drop);
end