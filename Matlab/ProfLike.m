% Profile Likelihood Generator
% Marisa Eisenberg (marisae@umich.edu) - 7/31/16 - updated 6-22-17

function profile = ProfLike(params,profindex,costfun,factor)
% Definitions
%   params = point in parameter space from which to profile (the parameter estimates)
%   profparam = index of the parameter to be profiled
%   factor = the fractional/percent range to profile the parameter over
%   costfun = this is a cost function of only the parameters (the full vector of them)
%     Everythng else (data, ICs, etc.) is fixed for the entire profile so is set outside when 
%     costfun is defined. In ProfLike, we'll put a wrapper on costfun that
%     will fix the profiled parameter.

% Setup
% factor = 0.5;
numpoints = 10;

% Costfun wrapper
fixcostfun = @(shortparams,profparamval)...
    costfun([shortparams(1:profindex-1),profparamval,shortparams(profindex:end)]);

% Profile
profrangeDown = linspace(params(profindex), params(profindex)*(1-factor),numpoints)'; 
profrangeUp = linspace(params(profindex), params(profindex)*(1+factor),numpoints)';
% split into up and down so we can use last fitted value as starting value for next run
profrange = [profrangeDown profrangeUp];
currfvals = [];
currparams = [];
currflags = [];
for i=1:2
    paramstemp = [params(1:profindex-1), params(profindex+1:end)];
    for j = 1:numpoints
        [i j] %track progress
        [paramstemp, fvaltemp, flagtemp] = fminsearch(@(p) fixcostfun(p,profrange(j,i)),paramstemp,optimset('MaxFunEvals',5000,'MaxIter',5000));
        currfvals = [currfvals; fvaltemp];
        currflags = [currflags; flagtemp];
        currparams = [currparams; [paramstemp(1:profindex-1),profrange(j,i),paramstemp(profindex:end)]]; %storing the profiled value too, so the output parameter values are easy to run the model with
    end
end

profile = [flipud([profrangeDown currfvals(1:numpoints) currflags(1:numpoints) currparams(1:numpoints,:)]);...
    [profrangeUp currfvals(numpoints+1:end) currflags(numpoints+1:end) currparams(numpoints+1:end,:)]];
end