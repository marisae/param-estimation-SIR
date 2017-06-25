% Profile Likelihood Generator
% Marisa Eisenberg (marisae@umich.edu) - 7/31/16 - updated 6-22-17

function profile = ProfLike(params,profparam,fitter,factor)
% Definitions
%   params = point in parameter space from which to profile (the parameter estimates)
%   profparam = index of the parameter to be profiled
%   factor = the fractional/percent range to profile the parameter over
%   fitter = this is a customized fminsearch that takes two arguments:
%     params and paramsfixedfcn, which tell it the starting parameters and
%     fixes the profiled parameter. Everythng else (data, ICs, likelihood
%     function, etc.) is fixed for the entire profile so is set outside when 
%     the fitter is defined.
%     e.g. fitter = @(params,paramfixedfcn) fminsearch(@(p) siwrML(times,p,paramfixedfcn,data,x0fcn,yfcn),params,optimset('MaxFunEvals',5000,'MaxIter',5000))


% Setup
% factor = 0.5;
numpoints = 10;

% Profile
profrangeDown = linspace(params(profparam), params(profparam)*(1-factor),numpoints)'; 
profrangeUp = linspace(params(profparam), params(profparam)*(1+factor),numpoints)';
% split into up and down so we can use last fitted value as starting value for next run
profrange = [profrangeDown profrangeUp];
currfvals = [];
currparams = [];
currflags = [];
for i=1:2
    paramstemp = params;
    for j = 1:numpoints
        [i j] %track progress
        if profparam==length(params) % got to be a nicer way to do this without an if statement, but think of it later
            paramfixedfcn = @(p) [p(1:profparam-1); profrange(j,i)];
        else
            paramfixedfcn = @(p) [p(1:profparam-1); profrange(j,i); p(profparam+1:end)];
        end
        
        %Quick note about how I'm doing the paramfixedfcn---it's not really
        %good practice to pass the optimizer a parameter that won't be
        %changing (the fixed parameter), since that direction in parameter
        %space will then have a completely flat likelhood. It's okay in
        %this case---the optimizer still runs fine for the other parameters
        %and it still converges, but be careful because some optimizers will not
        %work well with this! You could (probably should...) instead pass the profiled parameter
        %separately to the ProfLike function and then combine the profiled 
        %parameter with the other parameters in the cost function before 
        %you run the ODE.
        
        [paramstemp, fvaltemp, flagtemp] = fitter(paramstemp,paramfixedfcn);
        paramstemp = paramfixedfcn(abs(paramstemp));
        currfvals = [currfvals; fvaltemp];
        currflags = [currflags; flagtemp];
        currparams = [currparams; paramstemp'];
    end
end

profile = [flipud([profrangeDown currfvals(1:numpoints) currflags(1:numpoints) currparams(1:numpoints,:)]);...
    [profrangeUp currfvals(numpoints+1:end) currflags(numpoints+1:end) currparams(numpoints+1:end,:)]];