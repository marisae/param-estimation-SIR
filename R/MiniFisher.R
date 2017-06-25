# Simplified FIM
# Marisa Eisenberg (marisae@umich.edu) - 6-23-17 

MiniFisher = function(tspan,params,x0fcn,xfcn,yfcn){
  delta = 0.001 #percent change we'll use, default 0.1%
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