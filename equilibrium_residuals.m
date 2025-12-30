function F = equilibrium_residuals(vars, p, mu, c, N_total)
    lambda = vars(1);
    N = vars(2:end);

    num_age = numel(p);
    
    F = zeros(1, num_age+1);
    F(1) = p(1) * lambda - (mu(1) + c(1)) * N(1);
    
    for k = 2:num_age
        F(k) = p(k) * lambda - (mu(k) + c(k)) * N(k) + c(k-1) * N(k-1);
    end
    
    % Constraint: sum of N_k should equal N_total
    F(num_age+1) = sum(N) - N_total;
end