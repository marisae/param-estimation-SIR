% SIR Example Code Identifiability Practical
% Marisa Eisenberg (marisae@umich.edu) - 7-29-16 - updated 6-22-17

% This code uses: sirODE.m, sirCost.m, MiniFisher.m

clear

%% Load Data

times = 0:7:100;
data = [97;271;860;1995;4419;6549;6321;4763;2571;1385;615;302;159;72;34]; % simulated data


% For exploring what happens when we only have early data (Part 2, Prob 5)
% times = times(1:7);
% data = data(1:7);

%% Setup

% A note---I often find it's convenient to make small functions that take
% your parameters/data/model as input and return (for example) the initial 
% conditions or measurement equation for your model. 
% So we'll make two functions:

%   x0fcn: a function that takes in the current parameter values and
%   returns the initial conditions, i.e. I(0) = data(1)/k, S(0) = 1- I(0),
%   R(0) = 0. We'll assume the parameter order is: beta, gamma, k.

x0fcn = @(params) [1-data(1)/params(3); data(1)/params(3); 0];

%   yfcn: a function that takes in the current model output and parameter
%   values, and returns the measured cases, y = I*k. Note that k =
%   params(3), and I = x(:,2).

yfcn = @(x,params) x(:,2)*params(3);

% Note the "params" and "x" here are just dummy variables--they will be
% filled once we call these functions! However, "data(1)" is actually
% referencing the dataset we just loaded. Be careful--if you change the
% data later, this will not be reflected here unless you reload the x0 function!
% You can change this if you want, e.g. to make x0fcn = @(params,data) etc.


% Starting parameter values:
params = [0.4, 0.25, 80000]; %beta, gamma, k


%% Simulate and Plot the Model

[t,x] = ode45(@sirODE,times,x0fcn(params),[],params);
y = yfcn(x,params);

figure(1)
    set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial')
    hold on
    plot(t,y,'b','LineWidth',2);
    plot(times,data,'ko','LineWidth',2);
    %legend('Initial Simulation with Starting Parameters'); 
    ylabel('Infected Population');  
    xlabel('Time (days)');

%% Parameter Estimation

% Estimate the model parameters 
[paramests, fval] = fminsearch(@(p) sirCost(times,p,data,x0fcn,yfcn),params,optimset('Display','iter','MaxFunEvals',5000,'MaxIter',5000));
paramests = abs(paramests);

% Re-simulate the model with the final parameter estimates
[test,xest] = ode45(@sirODE,times,x0fcn(paramests),[],paramests);
yest = yfcn(xest,paramests);

% Plot results
figure(1)  % Model Fit
    set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial')
%     hold on
    plot(test,yest,'k','LineWidth',2);
    plot(times,data,'ko','LineWidth',2);
    legend('Initial Simulation with Starting Parameters','Data','Model with Parameter Estimates','Location','se'); 
    ylabel('Infected Population');  
    xlabel('Time (days)');
    
figure(2)  %Residuals of the fit
    set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial')
    hold on
    plot(test,(yest-data),'o');
    ylabel('Residuals');  
    xlabel('Time (days)'); 

    
%% Fisher Information Matrix

% Numerically approximate the FIM
FIM = MiniFisher(times,paramests,x0fcn,@sirODE,yfcn);

%Calculate & print rank of FIM
rank(FIM)

%Calculate SEs & CVs, print CVs
% SEs = sqrt(diag(((length(data) - length(paramests))/fval)*inv(FIM)));
% CVs = SEs./paramests





