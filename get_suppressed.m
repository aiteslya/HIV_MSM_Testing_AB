function [Pop_hiv,newly_sup,total_treat_ind]=get_suppressed(Pop_hiv,hiv_index,HIV_st,cfg,total_treat_ind,tcounter)
% individuals become suppressed
myname=mfilename;
if isequal(ModelConstants.unit,'week')
    year = 52;
end

newly_sup=[];

if numel(total_treat_ind)>0

    eta = 1-exp(-cfg.eta/year);
    
    for ind=total_treat_ind
        st=Pop_hiv.Data(hiv_index.status,ind);

        r1=rand;
        if r1<eta % got suppressed
            newly_sup=[newly_sup, ind];
            Pop_hiv.Data(hiv_index.date_sup,ind)=tcounter;
            if st==HIV_st.A1T
                Pop_hiv.Data(hiv_index.status,ind)=HIV_st.A1S;
                Pop_hiv.Data(hiv_index.date_sup,ind) = tcounter;
            elseif st==HIV_st.A23T
                Pop_hiv.Data(hiv_index.status,ind)=HIV_st.A23S;
                Pop_hiv.Data(hiv_index.date_sup,ind) = tcounter;
            elseif st==HIV_st.A45T
                Pop_hiv.Data(hiv_index.status,ind)=HIV_st.A45S;
                Pop_hiv.Data(hiv_index.date_sup,ind) = tcounter;
            elseif st==HIV_st.CT
                Pop_hiv.Data(hiv_index.status,ind)=HIV_st.CS;
                Pop_hiv.Data(hiv_index.date_sup,ind) = tcounter;
            elseif st==HIV_st.EAT
                Pop_hiv.Data(hiv_index.status,ind)=HIV_st.EAS;
                Pop_hiv.Data(hiv_index.date_sup,ind) = tcounter;
            elseif  st==HIV_st.LAT
                Pop_hiv.Data(hiv_index.status,ind)=HIV_st.LAS;
                Pop_hiv.Data(hiv_index.date_sup,ind) = tcounter;
            else
                error([myname,': someone who should not be suppressed is suppressed']);
            end
        end
    end

    total_treat_ind=setdiff(total_treat_ind,newly_sup);
end
end