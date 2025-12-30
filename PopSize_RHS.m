function dydt=PopSize_RHS(t,y,pars)
% this is the right hand side of ODE sstem tracking the changes of age
% brackets of a population 
% the population is aged 15-75 with each brackets 10 years large

lambda = pars(1); % entrance rates
c = pars(2); % age rate
mu = pars(3:8); % death rates

num_brack=6;
dydt = zeros(num_brack,1);

dydt(1,1)=lambda-(mu(1)+c)*y(1);

for k=2:num_brack
    dydt(k,1)=c*y(k-1,1)-(mu(k)+c)*y(k,1);
end

end