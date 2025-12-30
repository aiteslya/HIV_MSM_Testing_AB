function [Pop_hiv, PrEP_inds, inds_drop] = PrEPDrop(Pop_hiv,hiv_index,HIV_st,cfg, PrEP_inds)
    % this function simulates enrollment of susceptible individuals into
    % PrEP programme
    year_week = 52;% duration of year in weeks
    pre_drop_w=cfg.prep_drop/year_week; % convert probability per year into probability per day

    ind=PrEP_inds.Data;
    inds_drop = [];
    if numel(ind)>0
        for counter=1:1:numel(ind)
            r1=rand;
            if r1<pre_drop_w
                if Pop_hiv.Data(hiv_index.status, ind(counter)) == HIV_st.suscPrep % if was not infected while on PrEP
                    Pop_hiv.Data(hiv_index.status,ind(counter)) = HIV_st.susc;
                end
                % remove a person from the registry of PrEP users - this
                % affects their testing rate
                inds_drop = [inds_drop, ind(counter)];
            end
        end
    end

    PrEP_inds.Data = setdiff(PrEP_inds.Data, inds_drop);
    
end