% Model equations for the scaled SIR model
% Marisa Eisenberg 7-29-16 (marisae@umich.edu)

function dxdt = sirODE(t,x,params)
%x(1) = susceptibles s
%x(2) = infecteds i
%x(3) = recovered r

mu = 0;
beta = params(1);
gamma = params(2);

dxdt = zeros(3,1); %column vector for the state variables - can change to (2,1) and ditch the recovereds as well

dxdt(1) = mu - beta*x(1)*x(2) - mu*x(1);
dxdt(2) = beta*x(1)*x(2) - gamma*x(2) - mu*x(2);
dxdt(3) = gamma*x(2) - mu*x(3);


