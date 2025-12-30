function [Population,Rels_steady]=create_steady_pair(Population,pop_index,t_counter,Rels_steady,alive_ind,age_attr,sero_attr,Pop_hiv,hiv_index,all_casual_parts, N_pairs)
% function that attempts to create a pair
% find all individuals who can have a relationship 
% recall that Id_counter points to the next slot available for population
% with newly entering individual
s_fl=0;

% in the populatio of live individuals find all individuals that are single
a_ind = alive_ind.Data;
ind_single = a_ind(Population.Data(pop_index.sp1_id, a_ind) == -1);

% pick the first individual to form the pair, ensuring there are
% individuals who can be recruited
if isempty(ind_single)
    return;
end

ind1 = ind_single(randi(numel(ind_single)));
bin_age1 = Population.Data(pop_index.age_bin,ind1);

% access respective age mixing column
pref_distr=age_attr(:,bin_age1);

% determine the age group of the partner

bin_age2 = sample_from_distribution(pref_distr);

% find all single individuals within the bracket

pop_p2_age_ind = ind_single(Population.Data(pop_index.age_bin, ind_single) == bin_age2);

% determine HIV status of the partner (negative (S or S_{P}) or positive
% but undiagnosed)
% determine the HIV status of the first individual
st_ind1=Pop_hiv.Data(hiv_index.status,ind1);

pref_sero_distr = sero_attr(:, 1 + (st_ind1 > 7));

% determine the HIV  status of the partner

hiv_st2 = sample_from_distribution(pref_sero_distr);

% find all individuals who fall within chosen HIV status
% start with locating all individuals who are susceptible within the
% necessary age bracket

ind_susc=pop_p2_age_ind(Pop_hiv.Data(hiv_index.status,pop_p2_age_ind)<=7);

if hiv_st2==2 % diagnosed individual

    mask = ~ismember(pop_p2_age_ind, ind_susc);
    ind_infect = pop_p2_age_ind(mask);

    if numel(ind_infect)>0 % there are infected individuals
        pop_p2_ind=ind_infect;
    else
        pop_p2_ind=ind_susc;
    end
elseif hiv_st2==1
    pop_p2_ind=ind_susc;
else
    error('create_steady_pair: when selecting serosorting preference invalid option was selected');
end


if numel(pop_p2_ind)>0

    index = randi(numel(pop_p2_ind));
    ind2 = pop_p2_ind(index);
    % checks
   
    if ind1~=ind2 % not the same person
        id1=Population.Data(pop_index.id,ind1);
        id2=Population.Data(pop_index.id,ind2);

       
        % check whether a person whose id=id2 is a current casual partner
        % of person whose id=id1

        cond_exist = partnershipExists(all_casual_parts, id1, id2);
        if ~cond_exist
                s_fl=1;
                % note that the starting date and the duration is updated
                % outside of this function. This is inconvenient.
                Rels_steady.Data = [Rels_steady.Data ; [sort([id1, id2]) 0 0]];
                Population.Data(pop_index.nsteady, ind1) = Population.Data(pop_index.nsteady,ind1) + 1;
                Population.Data(pop_index.nsteady,ind2) = Population.Data(pop_index.nsteady,ind2) + 1;
                % insert new partners into the Population table
                % insert id2 as a new partner of id1
                Population.Data(pop_index.sp1_id,ind1)=id2;
                Population.Data(pop_index.sp1_stdate,ind1) = t_counter;
                % insert id1 as a new partner of id2
                Population.Data(pop_index.sp1_id,ind2)=id1;
                Population.Data(pop_index.sp1_stdate,ind2) = t_counter;
        end

    end

end

end