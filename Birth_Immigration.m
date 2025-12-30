function [Population, Pop_hiv, Id_counter, alive_ind, ind_new_bracket, new_infect3, new_diagn, new_suppr, infect_alive] = Birth_Immigration(Population,Pop_hiv,pop_index,hiv_index,lambda_w,tcounter,Id_counter,alive_ind, cfg_main, lambda_diagn_w, infect_alive, HIV_st)
                                                                                                                                                                    
% Function which adds new individuals to the population. These individuals
% are distributed using age of sexual debut data. The importation of
% infections from outside the Netherlands is taken into account through the
% addition of indviduals who are infected. 


% insert new individuals which are born according to poisson process with
% average lambda_w per week
% set their age to 15 years, set their date of birth 
    new_infect3 = [];
    new_diagn = [];
    new_suppr = [];

    year_w = 52;
    year_day = 365;
    age_edges = 15:10:75;
    week = 7;
    
    N_real = 200000;
    mult_fact = cfg_main.N/N_real;

    num_born=poissrnd(lambda_w);
    
    if num_born>0
        % find the first available num_born slates (accounting for
        % possibility new individuals may take place of individuals who
        % already left the populaiton
        ind = find(Population.Data(pop_index.id,:)==0, num_born);  
        % Generate probabilities to sample age of sexual debut from a distribution
        r_arr = rand(1, num_born);
        %Loop for each new person
        age_min = 31;
        age_max = 75;
        for people_count = 1:num_born
            index_pos = sample_from_distribution_v2(AgeBirthProp, r_arr(people_count));
            if index_pos < 17
                Population.Data(pop_index.age,ind(people_count))= index_pos + 14;
            else
                Population.Data(pop_index.age,ind(people_count)) = round(age_min + (age_max - age_min)*rand(1,1)-1); % linear distribution (assumed) for age above 30
            end

        end
        % Put into bins
        % Use discretize() function
        Population.Data(pop_index.age_bin, ind) = discretize(Population.Data(pop_index.age, ind), age_edges);
        Population.Data(pop_index.birth,ind) = (tcounter-mod(tcounter,year_w))/year_w-Population.Data(pop_index.age,ind); % tcounter is in weeks, needs to converted to years  
        new_ids=Id_counter:1:(Id_counter+num_born-1);
        Id_counter=Id_counter+num_born;
        Population.Data(pop_index.id,ind)=new_ids;
        Population.Data(pop_index.death,ind)=-1;
        Population.Data(pop_index.sp1_id,ind)=-1;
        Population.Data(pop_index.sp1_stdate,ind)=-40000;
        Population.Data(pop_index.nsteady,ind)=0;
        % list new individuals as having to sample propensity to acquire
        % casual partners
        ind_new_bracket=ind;
   
        alive_ind.Data=[alive_ind.Data, ind];

        %add hiv table record
        Pop_hiv.Data(hiv_index.id,ind)=new_ids;

        if tcounter>=(cfg_main.T_burn+0.5)*year_w
        % set up an array of entrance rates
            lambda_inf = mult_fact*cfg_main.import_undiagn_rate/year_w;
            lambda_diagn = lambda_diagn_w;
            lambda_art = lambda_diagn*cfg_main.prop_immigr_suppr;
            lambda_sus = lambda_d - lambda_inf - lambda_diagn - lambda_art;
            Lambda_arr = [lambda_inf lambda_diagn-lambda_art lambda_art lambda_sus];
            Lambda_arr = Lambda_arr./sum(Lambda_arr);
            % determine the statuses
            HIV_st_new = randsample(1:4, num_born,true, Lambda_arr);
            n_infect = sum(HIV_st_new == 1);
            n_diagn = sum(HIV_st_new == 2);
            n_suppr = sum(HIV_st_new == 3);
            n_susc = sum(HIV_st_new ==4);
        else
            n_infect = 0;
            n_diagn = 0;
            n_suppr = 0;
        end

        if (n_infect + n_diagn + n_suppr) == 0 % no individuals with pre-existing infection entered the population

            Pop_hiv.Data(hiv_index.status,ind)=0;
            Pop_hiv.Data(hiv_index.date_inf,ind)=-40000;
            Pop_hiv.Data(hiv_index.date_diagn,ind)=-40000;
            Pop_hiv.Data(hiv_index.date_art,ind)=-40000;
            Pop_hiv.Data(hiv_index.date_sup,ind)=-40000;
            Pop_hiv.Data(hiv_index.date_new_stage,ind)=-40000;
        else
            % decide who of new-entered people is in what stage
            % infected
            ind_select = ind;
            init_infect_seed_weight = [1/cfg_main.gamma_a1 1/cfg_main.gamma_a23 1/cfg_main.gamma_a45 1/cfg_main.gamma_c 1/cfg_main.gamma_ea 1/cfg_main.mu_la];
            % normalization of weights
            init_infect_seed_weight = init_infect_seed_weight./sum(init_infect_seed_weight);

            % define mean sojourn times
            gamma_a1W=cfg_main.gamma_a1/year_week;
            gamma_a23W=cfg_main.gamma_a23/year_week;
            gamma_a45W=cfg_main.gamma_a45/year_week;
            gamma_cW=cfg_main.gamma_c/year_week;
            gamma_eaW=cfg_main.gamma_ea/year_week;

            if n_infect>0
                ind_infect = datasample(ind_select,n_infect,'Replace',false);
                ind_select = setdiff(ind_select, ind_infect);

                new_infect3 = ind_infect;

                % set HIV status
                Pop_hiv.Data(hiv_index.status, ind_infect) = randsample(1:6,n_infect,true,init_infect_seed_weight);

                Pop_hiv.Data(hiv_index.date_diagn,ind_infect)=-40000;
                Pop_hiv.Data(hiv_index.date_art,ind_infect)=-40000;
                Pop_hiv.Data(hiv_index.date_sup,ind_infect)=-40000;
                
                for ind=ind_infect
                    % retrieve HIV status
                    st = Pop_hiv.Data(hiv_index.status,ind);
                    if st==HIV_st.A1 
                        Pop_hiv.Data(hiv_index.date_inf,ind) = (tcounter*week-1)/week; % note that -4 in each branch is added so that people do not get infected or diagnosed in this step
                    elseif st==HIV_st.A23 
                        Pop_hiv.Data(hiv_index.date_inf,ind) = (tcounter*week - year_day/cfg_main.gamma_a1 - randi(floor(1+year_day/cfg_main.gamma_a23)))/week; 
                    elseif st==HIV_st.A45 
                        Pop_hiv.Data(hiv_index.date_inf, ind) = (tcounter*week - year_day/cfg_main.gamma_a1 - year_day/cfg_main.gamma_a23 - randi(floor(1+year_day/cfg_main.gamma_a45)))/week;
                    elseif st==HIV_st.C 
                        Pop_hiv.Data(hiv_index.date_inf, ind) = (tcounter*week - year_day/cfg_main.gamma_a1 - year_day/cfg_main.gamma_a23 - year_day/cfg_main.gamma_a45 - randi(floor(1+year_day/cfg_main.gamma_c)))/week;
                    elseif st==HIV_st.EA 
                        Pop_hiv.Data(hiv_index.date_inf, ind) = (tcounter*week - year_day/cfg_main.gamma_a1 - year_day/cfg_main.gamma_a23 - year_day/cfg_main.gamma_a45 - year_day/cfg_main.gamma_c - randi(floor(1+year_day/cfg_main.gamma_ea)))/week;
                    elseif st==HIV_st.LA 
                        Pop_hiv.Data(hiv_index.date_inf, ind) = (tcounter*week - year_day/cfg_main.gamma_a1 - year_day/cfg_main.gamma_a23 - year_day/cfg_main.gamma_a45 - year_day/cfg_main.gamma_c - year_day/cfg_main.gamma_ea - randi(floor(1+year_day/cfg_main.mu_la)))/week;
                    else % this individual should not be taken care off in this branch
                        myname = mfilename;
                        error([myname, 'Time of infection is assigned to someone who was diagnosed']);
                    end
               
                    if st==HIV_st.A1 | st==HIV_st.A1D | st==HIV_st.A1T
                        Gamma=gamma_a1W;
                        % sample from an exponential distribution
                        r1 = rand;
                        dur = -log(1-r1)/Gamma;
                        Pop_hiv.Data(hiv_index.date_new_stage,ind) = tcounter+dur;
                    elseif st==HIV_st.A23 | st==HIV_st.A23D | st == HIV_st.A23T
                        Gamma=gamma_a23W;
                        % sample from an exponential distribution
                        r1 = rand;
                        dur = -log(1-r1)/Gamma;
                        Pop_hiv.Data(hiv_index.date_new_stage,ind)=tcounter+dur;
                    elseif st==HIV_st.A45 | st==HIV_st.A45D | st==HIV_st.A45T
                        Gamma=gamma_a45W;
                        % sample from an exponential distribution
                        r1 = rand;
                        dur= -log(1-r1)/Gamma;
                        Pop_hiv.Data(hiv_index.date_new_stage,ind)=tcounter+dur;
                    elseif st==HIV_st.C | st==HIV_st.CD | st==HIV_st.CT
                        Gamma=gamma_cW;
                        % sample from an Erlang distribution
                        mean_dur=1/Gamma;
                        dur=chronic_dur_sample(mean_dur,'erlang',80);
                        Pop_hiv.Data(hiv_index.date_new_stage,ind)=tcounter+dur;
                    elseif st==HIV_st.EA | st==HIV_st.EAD | st==HIV_st.EAT
                        Gamma=gamma_eaW;
                        % sample from an exponential distribution
                        r1 = rand;
                        dur=-log(1-r1)/Gamma;
                        Pop_hiv.Data(hiv_index.date_new_stage,ind)=tcounter+dur;
                     elseif st==HIV_st.LA | st==HIV_st.LAD | st==HIV_st.LAT
                        Pop_hiv.Data(hiv_index.date_new_stage,ind)=cfg_main.T*year +10;
                    end
                end
            else
                ind_infect = [];
            end

            % diagnosed

            if n_diagn>0
                ind_diagn =  datasample(ind_select,n_diagn,'Replace',false);%ind_select(unidrnd(numel(ind_select),[1, n_diagn]));
                ind_select = setdiff(ind_select, ind_diagn);

                new_diagn = ind_diagn;

                % set HIV status
                Pop_hiv.Data(hiv_index.status, ind_diagn) = randsample(HIV_st.A1D:1:HIV_st.LAD ,n_diagn,true,init_infect_seed_weight);
                Pop_hiv.Data(hiv_index.date_diagn,ind_diagn) = (tcounter*week - 2)/week;

                Pop_hiv.Data(hiv_index.date_art,ind_diagn)=-40000;
                Pop_hiv.Data(hiv_index.date_sup,ind_diagn)=-40000;
                
                for ind=ind_diagn
                    % retrieve HIV status
                    st = Pop_hiv.Data(hiv_index.status,ind);
                    if st==HIV_st.A1D 
                        Pop_hiv.Data(hiv_index.date_inf,ind) = (tcounter*week-1)/week; % note that -4 in each branch is added so that people do not get infected or diagnosed in this step
                    elseif st==HIV_st.A23D 
                        Pop_hiv.Data(hiv_index.date_inf,ind) = (tcounter*week - year_day/cfg_main.gamma_a1 - randi(floor(1+year_day/cfg_main.gamma_a23)))/week; 
                    elseif st==HIV_st.A45D 
                        Pop_hiv.Data(hiv_index.date_inf, ind) = (tcounter*week - year_day/cfg_main.gamma_a1 - year_day/cfg_main.gamma_a23 - randi(floor(1+year_day/cfg_main.gamma_a45)))/week;
                    elseif st==HIV_st.CD 
                        Pop_hiv.Data(hiv_index.date_inf, ind) = (tcounter*week - year_day/cfg_main.gamma_a1 - year_day/cfg_main.gamma_a23 - year_day/cfg_main.gamma_a45 - randi(floor(1+year_day/cfg_main.gamma_c)))/week;
                    elseif st==HIV_st.EAD 
                        Pop_hiv.Data(hiv_index.date_inf, ind) = (tcounter*week - year_day/cfg_main.gamma_a1 - year_day/cfg_main.gamma_a23 - year_day/cfg_main.gamma_a45 - year_day/cfg_main.gamma_c - randi(floor(1+year_day/cfg_main.gamma_ea)))/week;
                    elseif st==HIV_st.LAD 
                        Pop_hiv.Data(hiv_index.date_inf, ind) = (tcounter*week - year_day/cfg_main.gamma_a1 - year_day/cfg_main.gamma_a23 - year_day/cfg_main.gamma_a45 - year_day/cfg_main.gamma_c - year_day/cfg_main.gamma_ea - randi(floor(1+year_day/cfg_main.mu_la)))/week;
                    else % this individual should not be taken care off in this branch
                        myname = mfilename;
                        error([myname, 'Time of infection is assigned to someone who was diagnosed']);
                    end

                     if st==HIV_st.A1 | st==HIV_st.A1D | st==HIV_st.A1T
                        Gamma=gamma_a1W;
                        % sample from an exponential distribution
                        r1 = rand;
                        dur=ceil(-log(1-r1)/Gamma);
                        Pop_hiv.Data(hiv_index.date_new_stage,ind) = tcounter+dur;
                    elseif st==HIV_st.A23 | st==HIV_st.A23D | st==HIV_st.A23T
                        Gamma=gamma_a23W;
                        % sample from an exponential distribution
                        r1 = rand;
                        dur=ceil(-log(1-r1)/Gamma);
                        Pop_hiv.Data(hiv_index.date_new_stage,ind) = tcounter+dur;
                    elseif st==HIV_st.A45 | st==HIV_st.A45D | st==HIV_st.A45T
                        Gamma=gamma_a45W;
                        % sample from an exponential distribution
                        r1 = rand;
                        dur=ceil(-log(1-r1)/Gamma);
                        Pop_hiv.Data(hiv_index.date_new_stage,ind)=tcounter+dur;
                    elseif st==HIV_st.C | st==HIV_st.CD | st==HIV_st.CT
                        Gamma=gamma_cW;
                        % sample from an Erlang distribution
                        mean_dur=1/Gamma;
                        dur=ceil(chronic_dur_sample(mean_dur,'erlang',80));
                        Pop_hiv.Data(hiv_index.date_new_stage,ind)=tcounter+dur;
                    elseif st==HIV_st.EA | st==HIV_st.EAD | st==HIV_st.EAT
                        Gamma=gamma_eaW;
                        % sample from an exponential distribution
                        r1 = rand;
                        dur=ceil(-log(1-r1)/Gamma);
                        Pop_hiv.Data(hiv_index.date_new_stage,ind)=tcounter+dur;
                     elseif st==HIV_st.LA | st==HIV_st.LAD | st==HIV_st.LAT
                        Pop_hiv.Data(hiv_index.date_new_stage,ind)=cfg_main.T*year +10;
                    end
                end

               
            else
                ind_diagn = [];
            end

            % suppressed

            if n_suppr>0
                ind_sup = datasample(ind_select,n_suppr,'Replace',false);%ind_select(unidrnd(numel(ind_select),[1, n_suppr]));
                new_suppr = ind_sup;
                ind_select = setdiff(ind_select, ind_sup);

                % set HIV status
                Pop_hiv.Data(hiv_index.status, ind_sup) = HIV_st.CS;
                Pop_hiv.Data(hiv_index.date_inf,ind_sup) = (week*tcounter - 5)/week;
                Pop_hiv.Data(hiv_index.date_diagn,ind_sup) = (week*tcounter - 4)/week;

                Pop_hiv.Data(hiv_index.date_art,ind_sup) = (week*tcounter - 3)/week;
                Pop_hiv.Data(hiv_index.date_sup,ind_sup) = (week*tcounter - 2)/week;
                Pop_hiv.Data(hiv_index.date_new_stage,ind_sup)=-40000;
                           
            else
                ind_sup = [];
            end


            % susceptible

            ind_susc = ind_select;

            if numel(ind_susc)>0
                Pop_hiv.Data(hiv_index.status,ind_susc)=0;
                Pop_hiv.Data(hiv_index.date_inf,ind_susc)=-40000;
                Pop_hiv.Data(hiv_index.date_diagn,ind_susc)=-40000;
                Pop_hiv.Data(hiv_index.date_art,ind_susc)=-40000;
                Pop_hiv.Data(hiv_index.date_sup,ind_susc)=-40000;
                Pop_hiv.Data(hiv_index.date_new_stage,ind_susc)=-40000;
            end

        end

    else
        ind_new_bracket = [];
    end
end