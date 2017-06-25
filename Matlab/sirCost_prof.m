% Cost Function for the SIR model
% Marisa Eisenberg 7-29-16 (marisae@umich.edu)

function NLL = sirCost_prof(tspan,params,paramsfixedfcn,data,x0fcn,yfcn)
params = paramsfixedfcn(abs(params));
[t,x] = ode45(@sirODE,tspan,x0fcn(params),[],params);
y = yfcn(x,params);
% NLL = ((data - y))'*((data - y)); %OLS
NLL = (y)'*ones(length(y),1) - data'*log(y);  %Poisson ML