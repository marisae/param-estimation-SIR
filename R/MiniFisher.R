# Simplified FIM
# Marisa Eisenberg (marisae@umich.edu) - 6-23-17 

# This code calculates a simplified form of the Fisher information matrix (FIM), for use in
# testing the number of identifiable combinations. We'll mostly use this for looking at 
# structural identifiability, although you could use it more generally. The form of the FIM
# here highlights the connection to the sensitivity matrix and comes from working out the FIM
# in the case where you have normally distributed measurement error (i.e. mean = model, known
# variances for each time point). In that case, the FIM can be written as X' W X (where X' = 
# transpose of X, and W is a weighting matrix). Here, we've simplified that even further to the
# case where we take W = I (identity matrix), so that we can just look at approximate structural
# identifiability by passing in a very finely spaced tspan.

# Input definitions
# tspan = times for measurement
# params = point in parameter space to calculate the FIM at
# x0fcn = function used to calculate the initial conditions 
#   (oops, this is using data, which we aren't actually passing---should fix that eventually...)
# xfcn = ode function
# yfcn = measurement equation
# delta = percent change we'll use for calculating the derivatives, default 0.1%

# Function Output
# - A length(params) x length(params) sized FIM.

MiniFisher = function(tspan,params,x0fcn,xfcn,yfcn,delta = 0.001){
   
  X = c()
  
  for (j in 1:length(params)) {
    params1 = params
    params2 = params
    params1[j] = (1+delta)*params[j]
    params2[j] = (1-delta)*params[j]
    
    x1 = ode(x0fun(data,params1), tspan, xfcn, params1, method='ode45')
    x2 = ode(x0fun(data,params2), tspan, xfcn, params2, method='ode45')
    
    X = cbind(X, (yfcn(x1,params1) - yfcn(x2,params2))/(2*delta*params[j]))
    #this fills in the jth column of the design matrix with the sensitivities to parameter j at each time point.
  }
  #FIM (simplified w/o weighting term)
  FIM = t(X)%*%X 
}