# Profile Likelihood Generator
# Marisa Eisenberg (marisae@umich.edu) - 6-24-17

# Input definitions
# params = row vector of starting parameters (all, including the one to be profiled)
# profparam = index within params for the parameter to be profiled
#   ---reminder to make this allow you to pass the name instead later on
# costfun = cost function for the model - should include params, times, and data as arguments.
#   Note costfun doesn't need to be specially set up for fixing the profiled parameter, 
#   it's just the regular function you would use to estimate all the parameters
#   (it will get reworked to fix one of them inside ProfLike)
# times, data = data set (times & values, or whatever makes sense)
#   ---possibly change this so it's included in costfun and not a separate set of inputs? Hmm.
# perrange = the percent/fraction range to profile the parameter over (default is 0.5)
# numpoints = number of points to profile at in each direction (default is 10)

# Output
# A list with:
#   - profparvals: the values of the profiled parameter that were used
#   - fnvals: the cost function value at each profiled parameter value
#   - convergence: the convergence value at each profiled parameter value
#   - paramestvals: the estimates of the other parameters at each profiled parameter value

ProfLike = function(params,profindex,costfun,times,data,perrange=0.5,numpoints=10){

  # Setup
  # split into up and down so we can use last fitted value as starting value for next run
  profrangeDown = seq(params[profindex],params[profindex]*(1-perrange),length.out = numpoints)
  profrangeUp = seq(params[profindex],params[profindex]*(1+perrange),length.out = numpoints)
  profrange = cbind(profrangeDown, profrangeUp)

  currfvals = c() #probably bad practice to grow these from an empty array but meh
  currparams = c()
  currflags = c()
  
  # Make a little wrapper around costfun so that it's only fitting the non-profiled parameters
  fixcostfun = function(shortparams,profparam,profindex,times,data,costfun){
    paramstemp = append(shortparams,profparam,after=(profindex-1))
    costfun(params=paramstemp,times=times,data=data)
  }
  # Loop over all values for the profiled parameter and fit the remaining parameter
  print("Starting profile...")
  for (i in 1:2){
    shortparams = params[-profindex]
    for (j in 1:numpoints){
      print(c(i,j)) # track progress
      res = optim(shortparams,fn=fixcostfun,profparam=profrange[j,i],profindex=profindex,times=times,data=data,costfun=costfun,method='Nelder-Mead')
      shortparams = res$par #save current fitted params as starting values for next round
      fvaltemp = res$value
      flagtemp = res$convergence
      currfvals = c(currfvals,fvaltemp)
      currflags = c(currflags, flagtemp)
      currparams = rbind(currparams, append(res$par,profrange[j,i],after=profindex-1));  #recording the  full set of params (even though it includes the profparam), to make it easier to just run the model with the profile output after ProfLike is done (i.e. you don't have to recombine the profiled param with the rest to run the model)
    }
  }
  
  # clean up the lists so they go in nice order for returning
  profrange = c(profrangeDown[numpoints:1],profrangeUp)
  currfvals = c(currfvals[numpoints:1],currfvals[(numpoints+1):(2*numpoints)])
  currflags = c(currflags[numpoints:1],currflags[(numpoints+1):(2*numpoints)])
  
  names(currparams)[profindex] = names(params)[profindex]
  currparams = rbind(currparams[numpoints:1],currparams[(numpoints+1):(2*numpoints)])
  output = list("profparvals"=profrange,"fnvals"=currfvals,"convergence"=currflags,"paramestvals"=currparams)
}
