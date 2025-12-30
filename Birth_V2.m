function [Population,Pop_hiv,Id_counter,alive_ind,ind_new_bracket] = Birth_V2(Population,Pop_hiv,pop_index,hiv_index,lambda_w,tcounter,Id_counter,alive_ind, AgeBirthProp)
% insert new individuals which are born according to poisson process with
% average lambda_w per week
% set their age to 15 years, set their date of birth 
    year_w = 52;
    num_born=poissrnd(lambda_w);
    age_edges = 15:10:75;

    if num_born>0
        % find the first available num_born slates (accounting for
        % possibility new individuals may take place of individuals who
        % already left the populaiton
        ind=find(Population.Data(pop_index.id,:)==0,num_born);    
        %Generate probabilities to sample from distribution
        r_arr = rand(1, num_born);
        %Loop for each new person
        age_min = 31;
        age_max = 75;
        for people_count = 1:num_born
            index_pos = sample_from_distribution_v2(AgeBirthProp, r_arr(people_count));
            if index_pos < 17
                Population.Data(pop_index.age,ind(people_count))= index_pos + 14;
            else
                Population.Data(pop_index.age,ind(people_count)) = round(age_min + (age_max - age_min)*rand(1,1)); % linear distribution (assumed) for age above 30
            end

        end
        % Put into bins
        % Use discretize() function
        Population.Data(pop_index.age_bin, ind) = discretize(Population.Data(pop_index.age, ind), age_edges);
        Population.Data(pop_index.birth,ind)=(tcounter-mod(tcounter,year_w))/year_w-Population.Data(pop_index.age,ind); % tcounter is in days, needs to converted to years  
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
        Pop_hiv.Data(hiv_index.status,ind)=0;
        Pop_hiv.Data(hiv_index.date_inf,ind)=-40000;
        Pop_hiv.Data(hiv_index.date_diagn,ind)=-40000;
        Pop_hiv.Data(hiv_index.date_art,ind)=-40000;
        Pop_hiv.Data(hiv_index.date_sup,ind)=-40000;
        Pop_hiv.Data(hiv_index.date_new_stage,ind)=-40000;
    else
        ind_new_bracket = [];
    end
end