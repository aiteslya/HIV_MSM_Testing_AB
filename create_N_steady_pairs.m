function [Population, Rels_steady, s_fl] = create_N_steady_pairs(Population,pop_index,t_counter,Rels_steady,alive_ind,age_attr,sero_attr,Pop_hiv,hiv_index,all_casual_parts, HIV_st, N_pairs)
% function that attempts to create N_pairs pairs
% find all individuals who can have a relationship 
% the flag of how successful the function was: 0 - no pairs were created, 1
% - some pairs, but not all were created, 2 exactly the number of pairs
% that was desired was created

s_fl = 0;
tol = 1e-6;

% pre-compute the deliniation
% in the population of live individuals find all individuals that are single
a_ind = alive_ind.Data;
ind_single = a_ind(Population.Data(pop_index.sp1_id, a_ind) == -1);

% deliniate the population by HIV status: diagnosed and undiagnosed
mask_non_diagn = Pop_hiv.Data(hiv_index.status,ind_single) <= HIV_st.suscPrep;
mask_diagn = Pop_hiv.Data(hiv_index.status,ind_single) > HIV_st.suscPrep;

ind_non_diagn = ind_single(mask_non_diagn);
ind_diagn = ind_single(mask_diagn);
% indicator of a presence of diagnosed individuals
fl_diagn = numel(ind_diagn)>0;
% deliniate the groups into age cathegories
num_bins = 6;
for bin_counter = 1:num_bins

    ind_non_diagn_age{bin_counter} = ind_non_diagn(Population.Data(pop_index.age_bin,ind_non_diagn) == bin_counter);
    
    if fl_diagn

        ind_diagn_age{bin_counter} = ind_diagn(Population.Data(pop_index.age_bin,ind_diagn) == bin_counter);

    end
end

num_pairs = 0;
n_single = numel(ind_single);

% counter of iterations performed to create a pair
iter_counter = 0;
MAX_ITER = 10;

r1_arr = rand(numel(ind_single), 1);
r2_arr = rand(numel(ind_single), 1);

attempts_counter = 1;
attempts_thresh = N_pairs;

while num_pairs < N_pairs & n_single>1 & iter_counter < MAX_ITER
    % pick the first partner
    % ind1 = ind_single(randi(numel(ind_single)));
    
    if numel(ind_single)>0 % check that there are people to form partnerships
        % debug February 10 2025
         if sum(Population.Data(pop_index.steady_prop, ind_single)) > tol
             distr = Population.Data(pop_index.steady_prop, ind_single);
         else
             distr = ones(1, numel(ind_single));
         end

        %distr = ones(1, numel(ind_single));

        sample_index1 = sample_from_distribution_v2(distr, r1_arr(attempts_counter));
        ind1 = ind_single(sample_index1);
    
        if attempts_counter < attempts_thresh
            attempts_counter = attempts_counter + 1;
        else
           r1_arr = rand(numel(ind_single), 1);
           r2_arr = rand(numel(ind_single), 1);
           attempts_counter = 1;
        end
    
        bin_age1 = Population.Data(pop_index.age_bin,ind1);
    
        % access respective age mixing column
        pref_distr=age_attr(:,bin_age1);
        % determine the age group of the partner
        bin_age2 = sample_from_distribution(pref_distr);
       
        % set hiv_st2 = 1 as a default
        hiv_st2 = 1;
    
        if fl_diagn
            % determine HIV status of the partner (negative (S or S_{P}) or positive
            % but undiagnosed)
            % determine the HIV status of the first individual
            st_ind1 = Pop_hiv.Data(hiv_index.status,ind1);
            
            pref_sero_distr = sero_attr(:, 1 + (st_ind1 > 7));
            
            % determine the HIV  status of the partner
            
            hiv_st2 = sample_from_distribution(pref_sero_distr);
        end
    
        if fl_diagn & hiv_st2 == 2 % there are diagnosed individuals and we sample from this population
            par2_ind_arr = ind_diagn_age{bin_age2};
        else
            par2_ind_arr = ind_non_diagn_age{bin_age2};
        end
    
        % ind2 = par2_ind_arr(randi(numel(par2_ind_arr)));
        if numel(par2_ind_arr)>0    
            % remove propensity by age
             if sum(Population.Data(pop_index.steady_prop, par2_ind_arr)) > tol
                 distr = Population.Data(pop_index.steady_prop, par2_ind_arr);
             else
                 distr = ones(1,numel(par2_ind_arr));
             end
        
            %distr = ones(1,numel(par2_ind_arr));

            sample_index2 = sample_from_distribution_v2(distr, r2_arr(attempts_counter));
            ind2 = par2_ind_arr(sample_index2);
            % checks
           
            if ind1~=ind2 % not the same person
                id1 = Population.Data(pop_index.id,ind1);
                id2 = Population.Data(pop_index.id,ind2);
               
                % check whether a person whose id=id2 is a current casual partner
                % of person whose id=id1
        
                cond_exist = partnershipExists(all_casual_parts, id1, id2);
                if ~cond_exist
                        % note that the starting date and the duration is updated
                        % outside of this function. This is inconvenient.
                        Rels_steady.Data = [Rels_steady.Data ; [sort([id1, id2]) 0 0]];
                        Population.Data(pop_index.nsteady, ind1) = Population.Data(pop_index.nsteady,ind1) + 1;
                        Population.Data(pop_index.nsteady,ind2) = Population.Data(pop_index.nsteady,ind2) + 1;
                        % insert new partners into the Population table
                        % insert id2 as a new partner of id1
                        Population.Data(pop_index.sp1_id,ind1)=id2;
                        Population.Data(pop_index.sp1_stdate,ind1)=t_counter;
                        % insert id1 as a new partner of id2
                        Population.Data(pop_index.sp1_id,ind2)=id1;
                        Population.Data(pop_index.sp1_stdate,ind2)=t_counter;
        
                        hiv_st1 = Pop_hiv.Data(hiv_index.status,ind1);
                        hiv_st2 = Pop_hiv.Data(hiv_index.status,ind2);
        
                        % update the pools of available individuals
                        ind_single = setdiff(ind_single,[ind1 ind2]);
                        if hiv_st1 <= HIV_st.suscPrep
                            % no HIV diagnosis
                            ind_non_diagn_age{bin_age1} = setdiff(ind_non_diagn_age{bin_age1}, ind1);
                        else
                            % received HIV diagnosis
                            ind_diagn_age{bin_age1} = setdiff(ind_diagn_age{bin_age1}, ind1);
                        end
        
                        if hiv_st2 <= HIV_st.suscPrep
                            % no HIV diagnosis
                            ind_non_diagn_age{bin_age2} = setdiff(ind_non_diagn_age{bin_age2}, ind2);
                        else
                            % received HIV diagnosis
                            ind_diagn_age{bin_age2} = setdiff(ind_diagn_age{bin_age2}, ind2);
                        end
        
                        % update the counters
                        iter_counter = 0;
                        num_pairs = num_pairs + 1;
                        n_single = n_single - 2;
                else
                    iter_counter = iter_counter + 1;
                end
            
            else
                iter_counter = iter_counter + 1;
            end
        else
            iter_counter = iter_counter + 1;
        end
    else % there are no single individuals, further efforts do not have a point
        if num_pairs == N_pairs
            s_fl = 2;
        elseif num_pairs>0
            s_fl = 1;
        else
            s_fl = 0;
        end
        return;
    end 


end
if num_pairs == N_pairs
    s_fl = 2;
elseif num_pairs>0
    s_fl = 1;
end

end