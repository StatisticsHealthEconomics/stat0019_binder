# This files essentially does the same analysis as OpenBUGS

# NB: Remember to set the working directory to the one where the files are stored.
# You can do this by going on the bottom-right panel, selecting the tab "Files"
# (which is the default/first one), then navigate to "practicals", then
# "02_mcmc" and then "R-files". Then click on the right-most button
# "More" and select from the slider menu "Set As Working Directory"
#
# Alternatively, you can do this in the R terminal by typing the command
# setwd("PATH-TO-YOUR-FOLDER")
# The path changes depending on the operating system. In the case of the 
# Binder Virtual Machine, this would be
# setwd("~/practical/02_mcmc/R-files")
# or in a more verbose way
# setwd("/home/rstudio/practical/02_mcmc/R-files")
# (note that '~' is a shortcut for the home folder '/home/rstudio', in this case)

# Utility package
library(bmhe)

# Defines the number of simulations
nsim = 11000  # NB: BUGS does 1000 burn-in + 10000 simulations!

# Then defines the prior distribution for theta
alpha = 9.2
beta = 13.8
theta = rbeta(nsim, alpha, beta)
# We can visualise the prior distribution in several ways. For instance:
# 1. Make a histogram of the simulations
hist(theta)
# 2. Plot the analytic density
plot(
  seq(0,1,.001),dbeta(seq(0,1,.001),alpha,beta),t="l",xlab="theta",
  ylab="p(theta)",main="Beta(9.2,13.8) prior"
)

# As mentioned in class, this model is *conjugated*, so even BUGS doesn't
# really do any real MCMC. Rather we *know* what is the form of the posterior
# distribution p(theta | y) and we can simulate directly from it (**in this case**!)
npats = 20 # This is the sample size of the first trial
y = 15     # This is the observed number of successes in the first trial
# Simulates the posterior, which is theta | y ~ Beta(alpha+y, beta+npats-y)
post.theta = rbeta(nsim, alpha+y, beta+npats-y)

# Compare the prior and the posterior
# NB: type '?points' to see how to add to an existing plot in "base" R
# Also type '?legend' to see how to add a legend to a plot in "base" R
plot(
  seq(0,1,.001),dbeta(seq(0,1,.001),alpha,beta),t="l",
  xlab="theta",ylab="",col="blue",ylim=c(0,5.5)
)
points(
  seq(0,1,.001),dbeta(seq(0,1,.001),alpha+y,beta+npats-y),t="l",col="red"
)
legend("topright",c("Prior","Posterior"),col=c("blue","red"),bty="n",lty=1)

# We can also simulate from the "predictive" distribution.
# y.pred is simply another Binomial variable, with the same theta parameter
# only this time its values will be simulated from the *posterior*!
npats2 = 40 # This is the sample size for the second (hypothetical) trial
y.pred = rbinom(nsim,npats2,post.theta)
# "Critical" value for the new, hypothetical trial
n.crit = 25
P.crit = y.pred >= n.crit

# Then computes the summary statistics
sims=cbind(theta,post.theta,y.pred,P.crit)
stats(sims)

# And plots the relevant posterior/predictive distributions
hist(y.pred)
abline(v=n.crit,lwd=2)
