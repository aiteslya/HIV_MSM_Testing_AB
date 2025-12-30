function [Population, Pop_hiv, Id_counter, alive_ind, ind_new_bracket, new_infect3, new_diagn, new_suppr, infect_alive] = Birth_V2_Immigration(Population,Pop_hiv,pop_index,hiv_index,lambda_w,tcounter,Id_counter,alive_ind, cfg_main, lambda_diagn_w, infect_alive, HIV_st)
% updated birth function accounting for importation of new infections via
% importation

% additional parameters being passed: configuration cfg, yearly number
% of diagnosed individuals who have immigrated infected_distr, and registry
% of infected alive individuals in infect_alive.Data
% assumption about infected_distr: it is rescaled to the model population
% size in the main routine
    new_infect3 = [];
    new_diagn = [];
    new_suppr = [];    
    num_born=poissrnd(lambda_w);

    N_real = 200000;
    mult_fact = cfg.N/N_real;

    year_w = 52;
    year_day = 365;
    age_edges = 15:10:75;
    week = 7;

    if num_born>0
        % find the first available num_born slates (accounting for
        % possibility new individuals may take place of individuals who
        % already left the populaiton
        ind = find(Population.Data(pop_index.id,:)==0, num_born);  
        alive_ind.Data=[alive_ind.Data, ind];

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
   
        Pop_hiv.Data(hiv_index.id,ind_new_bracket)=new_ids;
        Pop_hiv.Data(hiv_index.ever_tested, ind_new_bracket) = 0;
        Pop_hiv.Data(hiv_index.testing_bin, ind_new_bracket) = 0;
        Pop_hiv.Data(hiv_index.test_rate, ind_new_bracket) = 0;

        % Infected births only occur when burn_fl == 1

        if burn_fl == 1

            % Run through all new births and determine if any are infected
            % Infection weights
            init_infect_seed_weight = [1/cfg_main.gamma_a1 1/cfg_main.gamma_a23 1/cfg_main.gamma_a45 1/cfg_main.gamma_c 1/cfg_main.gamma_ea 1/cfg_main.mu_la];
            % normalization of weights
            init_infect_seed_weight = init_infect_seed_weight./sum(init_infect_seed_weight);

            % Perform binomial sampling to see who is infected, and then assign
            % the infection stage using the weights related to duration of each
            % stage
            % N_infections = binornd(length(ind_new_bracket), cfg.prob_birth_infec, 1); % Number of new births with HIV
            % infected_ind = datasample(ind_new_bracket, N_infections,'Replace',false); % Randomly pick N_infections new births to have HIV
            
            infected_ind = ind_new_bracket(rand(1, num_born)<cfg_main.prop_born_infect);

            % temporary assignment to be updated later with the exact stage
            Pop_hiv.Data(hiv_index.status,infected_ind) = 1; 

            infect_alive.Data = [infect_alive.Data infected_ind];

            % Now determine stage of infection
            % Code allows for all stages to be possible for now
            % So, 1: A1, 2:A23, 3:A45, 4:C, 5:EA, 6:LA

           

            % Now we have the number of infected sorted, as well as the stage
            % of infection they are at. Now we need to determine who is
            % diagnosed and undiagnosed.
            % To do this, run through all new borns. If they are uninfected
            % proceed with old code, otherwise run code to see if they are
            % diagnosed
    
        
            for people_count = 1:length(ind_new_bracket)
                % some default properties for everyone

                if Pop_hiv.Data(hiv_index.status,ind_new_bracket(people_count)) == 0 % if susceptible
                    %add hiv table record
                    Pop_hiv.Data(hiv_index.date_inf,ind_new_bracket(people_count))=-40000;
                    Pop_hiv.Data(hiv_index.date_diagn,ind_new_bracket(people_count))=-40000;
                    Pop_hiv.Data(hiv_index.date_art,ind_new_bracket(people_count))=-40000;
                    Pop_hiv.Data(hiv_index.date_sup,ind_new_bracket(people_count))=-40000;
                    Pop_hiv.Data(hiv_index.date_new_stage,ind_new_bracket(people_count))=-40000;
                    

                else
                    % code for individuals who entered with HIV infections
                    % Assumption: those who enter the Netherlands diagnosed are
                    % suppressed
                   
                    if rand < cfg_main.prop_diagn_immig % Test to see if infected person is diagnosed
                        % Straight to suppressed
                        Pop_hiv.Data(hiv_index.status, ind_new_bracket(people_count))= 23;
                        % Determine time of infection, diagnosis, treatment,
                        % suppression. For people who enter the population
                        % diagnosed with HIV
                        Pop_hiv.Data(hiv_index.date_inf, ind_new_bracket(people_count)) = (tcounter*week - 1 - year_day/Pop_hiv.Data(hiv_index.test_rate, ind_new_bracket(people_count)) - year_day/cfg_main.theta_a1 - year_day/cfg_main.eta_a1)/week;
                        Pop_hiv.Data(hiv_index.date_diagn, ind_new_bracket(people_count)) = (tcounter*week - 1 - year_day/cfg_main.theta_a1 - year_day/cfg_main.eta_a1)/week;
                        Pop_hiv.Data(hiv_index.date_art, ind_new_bracket(people_count)) = (tcounter*week - 1 - year_day/cfg_main.eta_a1)/week;
                        Pop_hiv.Data(hiv_index.date_sup, ind_new_bracket(people_count)) = (tcounter*week - 1)/week;

                        % record in suppressed register
                        total_sup_ind = [total_sup_ind ind_new_bracket(people_count)];
                    else % Not diagnosed
                        % Determine time of infection
                        
                        % retrieve HIV status
                        st = randsample(1:6, 1,true, init_infect_seed_weight);
                        Pop_hiv.Data(hiv_index.status,ind_new_bracket(people_count)) = st;

                        if st==HIV_st.A1 
                            Pop_hiv.Data(hiv_index.date_inf, ind_new_bracket(people_count)) = (tcounter*week-1)/week; % note that -4 in each branch is added so that people do not get infected or diagnosed in this step
                        elseif st==HIV_st.A23 
                            Pop_hiv.Data(hiv_index.date_inf, ind_new_bracket(people_count)) = (tcounter*week - year_day/cfg_main.gamma_a1 - randi(floor(1+year_day/cfg_main.gamma_a23)))/week; 
                        elseif st==HIV_st.A45 
                            Pop_hiv.Data(hiv_index.date_inf, ind_new_bracket(people_count)) = (tcounter*week - year_day/cfg_main.gamma_a1 - year_day/cfg_main.gamma_a23 - randi(floor(1+year_day/cfg_main.gamma_a45)))/week;
                        elseif st==HIV_st.C 
                            Pop_hiv.Data(hiv_index.date_inf, ind_new_bracket(people_count)) = (tcounter*week - year_day/cfg_main.gamma_a1 - year_day/cfg_main.gamma_a23 - year_day/cfg_main.gamma_a45 - randi(floor(1+year_day/cfg_main.gamma_c)))/week;
                        elseif st==HIV_st.EA 
                            Pop_hiv.Data(hiv_index.date_inf, ind_new_bracket(people_count)) = (tcounter*week - year_day/cfg_main.gamma_a1 - year_day/cfg_main.gamma_a23 - year_day/cfg_main.gamma_a45 - year_day/cfg_main.gamma_c - randi(floor(1+year_day/cfg_main.gamma_ea)))/week;
                        elseif st==HIV_st.LA 
                            Pop_hiv.Data(hiv_index.date_inf, ind_new_bracket(people_count)) = (tcounter*week - year_day/cfg_main.gamma_a1 - year_day/cfg_main.gamma_a23 - year_day/cfg_main.gamma_a45 - year_day/cfg_main.gamma_c - year_day/cfg_main.gamma_ea - randi(floor(1+year_day/cfg_main.mu_la)))/week;
                        else % this individual should not be taken care off in this branch
                            myname = mfilename;
                            error([myname, 'Time of infection is assigned to someone who was diagnosed']);
                        end
                        
                        % record in the undiagnosed register
                        total_infect_undiagn_ind.Data = [total_infect_undiagn_ind.Data ind_new_bracket(people_count)];
                        
                        % Use "switch case" to determine time to progress to
                        % next HIV stage
        
                       
                        % For this section, we want Gamma values (rate of
                        % changing HIV infection stage) in [days^-1]:
                        
    
                        year_day = 365;
                        year_week = 52;
                        gamma_a1D=cfg_main.gamma_a1/year_day;
                        gamma_a23D=cfg_main.gamma_a23/year_day;
                        gamma_a45D=cfg_main.gamma_a45/year_day;
                        gamma_cD=cfg_main.gamma_c/year_day;
                        gamma_eaD=cfg_main.gamma_ea/year_day;
    
                        r1 = rand; % Need random number to sample
    
                        switch Pop_hiv.Data(hiv_index.date_inf, ind_new_bracket(people_count)) % Look at HIV status
                            case st == HIV_st.A1 % Acute 1 stage
                                duration = ceil(-log(1-r1)/gamma_a1D);
                                Pop_hiv.Data(hiv_index.date_new_stage,ind_new_bracket(people_count)) = tcounter + duration/week;
                            case st == HIV_st.A23 % Acute 2/3 stage
                                duration = ceil(-log(1-r1)/gamma_a23D);
                                Pop_hiv.Data(hiv_index.date_new_stage,ind_new_bracket(people_count)) = tcounter + duration/week;
                            
                            case st == HIV_st.A45 % Acute 4/5 stage
                                duration = ceil(-log(1-r1)/gamma_a45D);
                                Pop_hiv.Data(hiv_index.date_new_stage,ind_new_bracket(people_count)) = tcounter + duration/week;
    
                            case st == HIV_st.C % Chronic stage
                                duration = ceil(-log(1-r1)/gamma_cD);
                                Pop_hiv.Data(hiv_index.date_new_stage,ind_new_bracket(people_count)) = tcounter + duration/week;
    
                            case st == HIV_st.EA % Early aids stage
                                duration = ceil(-log(1-r1)/gamma_eaD);
                                Pop_hiv.Data(hiv_index.date_new_stage,ind_new_bracket(people_count)) = tcounter + duration/week;
    
                            case st == HIV_st.LA % Late aids stage
                                Pop_hiv.Data(hiv_index.date_new_stage,ind_new_bracket(people_count))=cfg_main.T*year_week +10; %%%%%%%%%% What is going on here? cfg_main.T*year_week = 520, so this is always 530?
    
                        end
                    end
                end
            end
        else % if burn_fl == 0
             % Still need to assign uninfected 
            
        	Pop_hiv.Data(hiv_index.status,ind_new_bracket)=0;
            Pop_hiv.Data(hiv_index.date_inf,ind_new_bracket)=-40000;
            Pop_hiv.Data(hiv_index.date_diagn,ind_new_bracket)=-40000;
            Pop_hiv.Data(hiv_index.date_art,ind_new_bracket)=-40000;
            Pop_hiv.Data(hiv_index.date_sup,ind_new_bracket)=-40000;
            Pop_hiv.Data(hiv_index.date_new_stage,ind_new_bracket)=-40000;

            
        end
    else
        ind_new_bracket = [];
    end
end