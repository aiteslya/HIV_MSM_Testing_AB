function [Population, ind_new_bracket] = Age(Population,pop_index,alive_ind)
    
    % ageing all individuals who are alive by 1 year
    min_age=15;
    max_age=75;
    band_width=10;
    age_edge=min_age:band_width:max_age;
    
    % Vectorized increment of ages
    Population.Data(pop_index.age, alive_ind.Data) = Population.Data(pop_index.age, alive_ind.Data) + 1;
    
    % Pre-computation of condition of being on an age edge
    has_aged_to_new_band = ismember(Population.Data(pop_index.age,alive_ind.Data), age_edge);
    
    % Identification of individuals who moved to a new age bracket (i.e.
    % all these who are not 75 yet
    below_limit = Population.Data(pop_index.age,alive_ind.Data) <= age_edge(end-1);
    eligible_for_new_bracket = and(has_aged_to_new_band, below_limit);
    
    ind_new_bracket = alive_ind.Data(eligible_for_new_bracket);
    
    % Incrementing the age bracket for these individuals
    Population.Data(pop_index.age_bin, ind_new_bracket) = Population.Data(pop_index.age_bin, ind_new_bracket) + 1;

end