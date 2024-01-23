# Uses the package 'here' to set the working directory to the correct one
setwd(here::here("02_mcmc"))

# Alternatively, you can do this by going on the bottom-right panel, selecting 
# the tab "Files" (which is the default/first one), then navigate to "practicals", 
# then "02_mcmc" and then "R-files". Then click on the right-most button
# "More" and select from the slider menu "Set As Working Directory"
#
# Or, you can do this in the R terminal by typing the command
# setwd("PATH-TO-YOUR-FOLDER")
# The path changes depending on the operating system. In the case of the 
# Binder Virtual Machine, this would be
# setwd("~/practical/02_mcmc/")
# or in a more verbose way
# setwd("/home/rstudio/practical/02_mcmc/")
# (note that '~' is a shortcut for the home folder '/home/rstudio', in this case)


# Creates the list of data
data=list(
a = 9.2, b = 13.8, # prior parameters
y = 15,            # number of successes
m = 20,            # number of trials
n = 40,            # future number of trials
ncrit = 25)        # critical value of future successes

# First load the 'R2OpenBUGS' package
library(R2OpenBUGS)

# Then calls the function 'bugs' to call 'BUGS' in the background
# and perform the MCMC analysis
model = bugs(
  data=data,
  parameters.to.save=c("y.pred","theta","P.crit"),
  inits=NULL,
  # The model file is define externally - using here::here ensures that we 
  # pick the correct full path to the model file
  model.file=here::here("02_mcmc","drug-MCMC.txt"),
  n.chains=2,
  n.iter=5000,
  n.burnin=0,
  DIC=TRUE
)


## Or can use a R function to define the model code
model.code=function(){
  theta ~ dbeta(a,b)                   # prior distribution
  y ~ dbin(theta,m)                    # sampling distribution
  y.pred ~ dbin(theta,n)               # predictive distribution
  P.crit <- step(y.pred - ncrit + 0.5) # =1 if y.pred >= ncrit,
}
model <- bugs(
  data=data,parameters.to.save=c("y.pred","theta","P.crit"),
  inits=NULL,n.chains=2,n.iter=5000,n.burnin=0,DIC=TRUE,
  # This specifies the model code as the function 'model.code'
  model.file=model.code,
)

# Visualise the results
print(model$summary,digit=3)

# And uses bmhe functions to graphically show convergence etc
bmhe::traceplot(model,"theta")
bmhe::traceplot(model,"y.pred")
 
# Checks for autocorrelation
acf(model$sims.list$theta,lwd=4,main="theta",col="red")
acf(model$sims.list$y.pred,lwd=4,main="y.pred",col="red")
bmhe::diagplot(model,what="Rhat")
bmhe::diagplot(model,what="n.eff")

# Plots the posterior distributions
bmhe::posteriorplot(model,plot="bar")

# Manually computes PSR (up to correction factor!) and n_eff using the simulations from the model
# Use the first node
x=model$sims.array[,,1]

# Defines the relevant quantities
# Number of simulations (rows of x)
nsims=nrow(x)
# Number of chains (columsn of x)
nchains=ncol(x)
# Vector of variance for simulated values in each chain
sigma2.hat=apply(x,2,var)
# Vector of mean for simulated values in each chain
mu.hat=apply(x,2,mean)
W=mean(sigma2.hat)
B=nsims*var(mu.hat)

# Now can compute (approximated) PSR (Rhat) and neff
Rhat_approx=sqrt( (1/W) * ( ((nsims-1)/nsims) * W + ((nchains+1)/(nsims*nchains)) * B ))
neff=nchains * nsims * min((1/B)*((nsims-1)/nsims) * W + (1/nsims) * B, 1)


# NB: If you like, you can use jags, rather than BUGS. 
# The process is *almost* identical... 
#
# First, load the package 'R2jags' (instead of 'R2OpenBUGS')
library(R2jags)
#
# In the Binder VM, this is already installed, so you are good to go.
# If you're doing this on your own machine, you must install first
# jags and then in R run the following command:
# install.packages("R2jags")
#
# Then substitute the command bugs(...) with the command jags(...), like this
model2 = jags(
  data=data,
  parameters.to.save=c("y.pred","theta","P.crit"),
  inits=NULL,
  # The model file is define externally - using here::here ensures that we 
  # pick the correct full path to the model file
  model.file=here::here("02_mcmc","drug-MCMC.txt"),
  n.chains=2,
  n.iter=5000,
  n.burnin=0,
  DIC=TRUE
)
#
# (I'm saving the output as 'model2' as I want to compare with the BUGS
# output. The name is irrelevant and you could choose whatever you want...).
#
# Notice that the two objects 'model' (saved by BUGS) and 'model2' (saved by jags)
# are *very* similar, but not identical...
#
# These commands show what's inside each of the two model outputs
names(model)
names(model2)
#
# As you can see, these are slightly different. Jags stores the results generated
# by running the MCMC algorithm in a sub-object, called 'BUGSoutput'. 
# So if you want to access the variables stored in there, you need to do it
# slightly differently than shown above for BUGS. So instead of 
acf(model$sims.list$theta,lwd=4,main="theta",col="red")
# 
# you need to call the slightly revised command
acf(model2$BUGSoutput$sims.list$theta,lwd=4,main="theta",col="red")
#
# See also https://egon.stats.ucl.ac.uk/static/stat0019/practical/02_mcmc/bugs-vs-jags.html
# for more (minor) differences between BUGS and jags
