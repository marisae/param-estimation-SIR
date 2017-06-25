% Extra code for running profile likelihoods with the SIR model. (Part 2, Prob 4)
% Note this should be run after running SIR_Example_Main!
% Marisa Eisenberg (marisae@umich.edu) - 7/31/16 - updated 6-22-17

% This code uses: ProfLike.m, sirCost.m, sirODE.m, SIR_Example_Main.m

paramlist = {'\beta','\gamma','k'};
profiles = [];

% Wrapper function for parameter estimation
costfun = @(p) sirCost(times,p,data,x0fcn,yfcn);
    % This is just me being a little lazy---made a wrapper since that's an 
    % easy way to pass everything you need to fit the parameters (e.g. 
    % data, initial conditions, measurement eqn, etc.), all in one blob.

threshold = chi2inv(0.95,length(paramests))/2 + fval;
profrange = 0.25; %percent range for profile to run across

for i=1:length(paramests)
    %Generate a profile for parameter i, using paramests as the starting
    %value and the fitter to do the parameter estimation:
    profiles(:,:,i) = ProfLike(paramests,i,costfun,profrange);
        % each profile has columns: profiled parameter value, resulting
        % cost-function (e.g. RSS) value, any flags from the optimizer, and
        % then columns for each of the other parameter estimates.
    
    %plot profile
    figure(i)
    set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial')
    hold on
    plot(profiles(:,1,i),profiles(:,2,i),'k','LineWidth',2)
    plot(paramests(i),fval,'r*','LineWidth',2)
    plot([profiles(1,1,i) profiles(end,1,i)],[threshold threshold],'r--')
    xlabel(paramlist{i})
    ylabel('Cost Function Value')
    
    %plot parameter relationships
    figure(10+i)
    set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial')
    hold on
    plot(profiles(:,1,i),profiles(:,4:end,i),'LineWidth',2)
    plot(paramests(i),paramests,'r*')
    xlabel(paramlist{i})
    ylabel('Estimated Parameter Value')
    legend(paramlist)
end

