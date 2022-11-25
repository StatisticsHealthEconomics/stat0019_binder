# Uses the package 'here' to set the working directory to the correct one
setwd(here::here("practical","02_mcmc"))

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
  model.file=here::here("practical","02_mcmc","drug-MCMC.txt"),
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
  data=dataJags,parameters.to.save=c("y.pred","theta","P.crit"),
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
