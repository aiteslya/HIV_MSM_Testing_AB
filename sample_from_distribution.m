function sampled_index = sample_from_distribution(distribution)
    % Sample from a given distribution 
    %
    % Parameters:
    %   distribution - A probability distribution array where elements should sum to 1
    %
    % Returns:
    %   sampled_index - The sampled index according to the distribution
    %
    % Raises:
    %   An error if the distribution is negative or does not sum to 1

    % Validate the distribution
    if any(distribution < 0)
        error('Distribution contains negative probabilities.');
    end
    if abs(sum(distribution) - 1) > 1e-6
        error('Sum of probabilities in the distribution does not equal 1.');
    end

    % Generate a random number between 0 and 1
    random_number = rand;

    % Compute the cumulative sum of the distribution
    cumulative_sum = cumsum(distribution);

    % Find the index of the first cumulative sum value that is greater than or equal to the random number
    bins = 1:numel(distribution);
    idx = min(bins(cumulative_sum >= random_number));
    
    sampled_index = idx;  
end
