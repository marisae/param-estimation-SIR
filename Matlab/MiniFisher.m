% Simplified FIM
% Marisa Eisenberg 7-29-16 (marisae@umich.edu)

function FIM = MiniFisher(tspan,params,x0fcn,xfcn,yfcn)

delta = 0.001; %percent change we'll use, default 0.1%
X = [];

for j=1:length(params)
    params1 = params;
    params2 = params;
    params1(j) = (1+delta)*params(j);
    params2(j) = (1-delta)*params(j);
    
    [t x1] = ode45(xfcn,tspan,x0fcn(params1),[],params1);
    [t x2] = ode45(xfcn,tspan,x0fcn(params2),[],params2);
    
    X = [X, (yfcn(x1,params1) - yfcn(x2,params2))./(2*delta*params(j))];
    %this fills in the jth column of the design matrix with the sensitivities to parameter j 
    %at each time point.
end

%In case the initial conditions are fixed
% X = X(2:end,:);

%FIM (simplified w/o weighting term)
FIM = X'*X;





