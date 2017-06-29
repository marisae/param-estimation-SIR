# SIR model example
# Marisa Eisenberg (marisae@umich.edu) - 6-23-17 

#### Load all the things ####
library(deSolve)
# library('Bhat')
library(Matrix)
source('MiniFisher.R')
source('ProfLike.R')
  
#### Load Data ####
times = c(0, 7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 84, 91, 98)
data = c(97, 271, 860, 1995, 4419, 6549, 6321, 4763, 2571, 1385, 615, 302, 159, 72, 34)

#shortened version for seeing how truncated data affects the estimation
# times = times[1:7]
# data = data[1:7]
  
#### Set initial parameter values ####
params = c('beta'=0.4,'gamma'=0.25, 'kappainv'=1/80000)
  # fitting 1/kappa because kappa is so giant that it slows optim down a bunch
  # (plus kappa has a huge potential range, but 1/kappa is really just between 0 and 1)
   
#### ODE ####
SIRode <- function(t, x, params){
  S = x[1]
  I = x[2]
  R = x[3]
  
  b = params[1]
  g = params[2]
  
  dS = -b*S*I
  dI = b*S*I - g*I
  dR = g*I
  
  list(c(dS, dI, dR))
}
   
#### Function to set initial conditions ####
# makes it easier to update the initial conditions whenever we change the parameter values
x0fun = function(data,params) {
  x0 = c(1-(data[1]*params[3]), data[1]*params[3], 0)
  names(x0) = c('S0','I0','R0')
  x0}
  
#### Measurement Equation ####
yfun = function(odeSim, params){odeSim[,3]/params[3]} 
  
#### Simulate the model ####
xinit <- ode(x0fun(data,params), times, SIRode, params, method='ode45')
plot(times, yfun(xinit,params), type='l')
points(times, data)
  
#### Likelihood function ####
SIRML=function(params,times,data){
  params = abs(params)
  # Simulate model
  xcurr = ode(x0fun(data,params), times, SIRode, params, method='ode45')
  
  # Measurement equation
  y = yfun(xcurr,params)
  
  # Negative Log Likelihood (NLL)
  NLL =  sum(y) - sum(data*log(y)) # Poisson ML
    # note this is a slightly shortened version--there's an additive constant term missing but it 
    # makes calculation faster and won't alter the threshold. Alternatively, can do:
  # NLL = -sum(log(dpois(round(data),round(y)))) # the round is b/c Poisson is for (integer) count data
    # this can also barf if data and y are too far apart because the dpois will be ~0, which makes the log angry
  
  # ML using normally distributed measurement error (least squares)
  # NLL = -sum(log(dnorm(data,y,0.1*mean(data)))) # example WLS assuming sigma = 0.1*mean(data)
  # NLL = sum((y - data)^2)  # alternatively can do OLS but note this will mess with the thresholds 
  #                             for the profile! This version of OLS is off by a scaling factor from
  #                             actual LL units.
  
  # return(NLL) 
}
  
#### Estimate parameters & plot ####
res = optim(params,fn=SIRML,times=times,data=data)#,method='Nelder-Mead')
paramests = res$par

xest = ode(x0fun(data,paramests), times, SIRode, paramests, method='ode45')
plot(times, yfun(xest,paramests), type='l')
points(times, data)
 
#### Calculate the simplified Fisher Information Matrix (FIM) ####
FIM = MiniFisher(times,paramests,x0fun,SIRode,yfun)
rankMatrix(FIM)[1]
# qr(FIM)$rank
  
#### Generate profile likelihoods and confidence bounds ####
profiles = list()
threshold = qchisq(0.95,length(paramests))/2 + res$value
perrange = 0.25 #percent range for profile to run across

for (i in 1:length(paramests)){
  profiles[[i]] = ProfLike(paramests,i,SIRML,times,data,perrange=perrange)
  
  #plot profile
  plot(profiles[[i]]$profparvals, profiles[[i]]$fnvals, type='l',
       xlab=names(params)[i], ylab="Negative Log Likelihood",
       ylim=c(min(profiles[[i]]$fnvals),max(c(profiles[[i]]$fnvals,threshold))))
  points(paramests[i], res$value)
  abline(h=threshold, col="red", lty="dashed")
  
  #plot parameter relationships
  # matplot(profiles[[i]]$profparvals, profiles[[i]]$paramestvals, type='l',xlab=names(params)[i],ylab="Estimated Parameter Value")
  # points(paramests[i], paramests)
  # legend("topleft",inset=.05,names(params),col=1:length(params),cex=0.8,fill=seq_len(length(params)))
}




