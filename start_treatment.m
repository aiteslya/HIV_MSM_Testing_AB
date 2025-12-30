function [Pop_hiv,newly_treat,total_diagn_ind]=start_treatment(Pop_hiv,hiv_index,HIV_st, total_diagn_ind,tcounter, theta_year, cfg)
% start treatment
    newly_treat=[];
    unit = ModelConstants.unit;
    switch unit
        case 'week'
            year = 52;
    end

    if any(Pop_hiv.Data(hiv_index.status,total_diagn_ind)<HIV_st.A1D | Pop_hiv.Data(hiv_index.status,total_diagn_ind)>HIV_st.LAD)
        error([myname,': someone who should not be treated is treated']);
    end

    theta = 1-exp(-theta_year*cfg.treat_red_fact/year);

    for ind=total_diagn_ind
        st=Pop_hiv.Data(hiv_index.status,ind);
              
        r1=rand;
        
        if r1<theta % started treatment
            newly_treat=[newly_treat, ind];
            Pop_hiv.Data(hiv_index.date_art,ind)=tcounter;
            if st==HIV_st.A1D
                Pop_hiv.Data(hiv_index.status,ind)=HIV_st.A1T;
            elseif st==HIV_st.A23D
                Pop_hiv.Data(hiv_index.status,ind)=HIV_st.A23T;
            elseif st==HIV_st.A45D
                Pop_hiv.Data(hiv_index.status,ind)=HIV_st.A45T;
            elseif st==HIV_st.CD
                Pop_hiv.Data(hiv_index.status,ind)=HIV_st.CT;
            elseif st==HIV_st.EAD
                Pop_hiv.Data(hiv_index.status,ind)=HIV_st.EAT;
            elseif  st==HIV_st.LAD
                Pop_hiv.Data(hiv_index.status,ind)=HIV_st.LAT;
            else
                error([myname,': someone who should not be treated is treated']);
            end
        end
    end
    
    total_diagn_ind=setdiff(total_diagn_ind,newly_treat);
        
end
