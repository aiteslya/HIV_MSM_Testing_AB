function prop=set_casual_prop(distr,r1)
% Function set_casual_prop samples propensity of individual whose ID=id
% from the distribution of propensities provided by row vector distr,
% containing integer values greater or equal to 0. r1 is a random number
% between 0 and 1 supplied by the user
    % find index of individual whose ID=id
    %ind=find(pop(pop_index.id,:)==id);

    if any(distr<0)
         error('set_steady_prop: A value in the inputted distribution, "distr", is negative');
    end


    Cum_distr = cumsum(distr);
    % 
    % if abs(Cum_distr(end)-1)>1e-6
    %     error('set_steady_prop: Inputted probability distribution, "distr", does not sum to 1')
    % end

    indices = 1:numel(Cum_distr);
    j = min(indices(Cum_distr > r1*sum(distr)));
    prop = j - 1; % minus 1 is necessary because propensity can be zero

    %find(cumulative_sum >= sum(distribution)*random_number, 1, 'first');

end