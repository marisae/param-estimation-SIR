% Cost Function for the SIR model
% Marisa Eisenberg 7-29-16 (marisae@umich.edu)

function NLL = sirCost(tspan,params,data,x0fcn,yfcn)
params = abs(params);
[t,x] = ode45(@sirODE,tspan,x0fcn(params),[],params);
y = yfcn(x,params);
NLL = (y)'*ones(length(y),1) - data'*log(y);  %Poisson ML
    % missing an additive constant term but it doesn't change the threshold
    % (and makes calulation easier/faster)

% Least squares - note these are off by a scaling factor so they'll change
% the threshold
% NLL = ((data - y))'*((data - y)); %OLS
% NLL = sum((data-y).^2); % Also OLS! :) 
% or can use something like this and be in LL units (note requires choosing a sigma!):
% NLL = -sum(log(normpdf(data,y,0.1*mean(data))));