function [Population, Rels_steady]=create_steady(Population,pop_index, tcounter,Rels_steady,alive_ind,age_attr,sero_attr,Pop_hiv,hiv_index,all_casual_parts, HIV_st, cfg, steady_dur_distr)
% create new steady pairs (per day)

% determine how many pairs can be
% created - look for individuals who are alive and do not have a steady
% partner
num_pairs = floor(sum(Population.Data(pop_index.nsteady,alive_ind.Data)<1)/2);

year_week = 52; % duration of the year in weeks
week = 7; % duration of week in days
sigma_w = cfg.sigma/year_week;% convert from probability per year to probability per week

r1_arr=rand(1,num_pairs);

% calculate the number of pairs that will be created
num_pairs_real=sum(r1_arr<sigma_w);

if num_pairs_real>0
    
    [Population, Rels_steady, ~] = create_N_steady_pairs(Population,pop_index,tcounter,Rels_steady,alive_ind,age_attr,sero_attr,Pop_hiv,hiv_index,all_casual_parts, HIV_st, num_pairs_real);
  
    % new code
    % Generate durations of steady partnerships that were created

    j_values = randsample(1:numel(steady_dur_distr),num_pairs_real,true,steady_dur_distr)/week;

    % Assign durations
    Rels_steady.Data((end-num_pairs_real+1):end, 3) = tcounter;
    Rels_steady.Data((end-num_pairs_real+1):end, 4) = tcounter + j_values;
    
end