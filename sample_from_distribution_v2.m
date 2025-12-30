function sampled_index = sample_from_distribution_v2(distribution, random_number)
    % Sample from a given distribution using an external random number
    %
    % Parameters:
    %   distribution - A probability distribution array where elements should sum to 1
    %   random_number - A pre-generated random number used for sampling
    %
    % Returns:
    %   sampled_index - The sampled index according to the distribution
    %
    % Raises:
    %   An error if the distribution is negative 

    % Validate the distribution
    if any(distribution < 0)
        error('Distribution contains negative probabilities.');
    end
    if sum(distribution) < 1e-6
        error('Sum of probabilities in the distribution is too close to 0.');
    end

    % Compute the cumulative sum of the distribution
    cumulative_sum = cumsum(distribution);

    % Find the index of the first cumulative sum value that is greater than or equal to the random number
    sampled_index = find(cumulative_sum >= sum(distribution)*random_number, 1, 'first');
end