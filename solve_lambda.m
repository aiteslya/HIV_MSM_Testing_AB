function [lambda_solution, N_solutions] = solve_lambda(N_total, mu, c, p)
% Define known parameters
%N_total = 5000; % Total N
% mu = [0.000305832390861   0.000489313086286   0.000914634659514   0.002352909483975 0.006753431712328   0.018654070891175]; % mu_k values

% Define initial guess for lambda
lambda_guess = 1.0;


% Initial guess for N_k
N_guess = N_total * p; % Proportional to p

% Combine initial guesses
initial_guess = [lambda_guess, N_guess];

% Solve using fsolve
options = optimoptions('fsolve', 'Display', 'none', 'TolFun', 1e-9, 'TolX', 1e-9);

[sol, fval, exitflag] = fsolve(@(vars) equilibrium_residuals(vars, p, mu, c, N_total), initial_guess, options);

% Extract lambda and N_k from solution
lambda_solution = sol(1);
N_solutions = sol(2:end);

% Display results
% fprintf('Lambda: %f\n', lambda_solution);
% fprintf('N_k values: \n');
% disp(N_solutions);

% Define function to compute residuals for fsolve


end