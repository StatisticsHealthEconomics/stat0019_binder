# This files essentially does the same analysis as OpenBUGS

# NB: Remember to set the working directory to the one where the files are stored.
# You can do this by going on the bottom-right panel, selecting the tab "Files"
# (which is the default/first one), then navigate to "practicals", then
# "01_monte-carlo" and then "R-files". Then click on the right-most button
# "More" and select from the slider menu "Set As Working Directory"
#
# Alternatively, you can do this in the R terminal by typing the command
# setwd("PATH-TO-YOUR-FOLDER")
# The path changes depending on the operating system. In the case of the 
# Binder Virtual Machine, this would be
# setwd("~/practical/01_monte-carlo/R-files")
# or in a more verbose way
# setwd("/home/rstudio/practical/01_monte-carlo/R-files")
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

# Now simulate the values for the trial
npats = 20 # number of patients in the trial
y = rbinom(nsim,npats,theta)
# Computes the probability of 15 or more successes in each simulation
threshold = 14.5 # y is only integers, so if the threshold is 13.5, it means
# there are 15 or more successes in a single simulation
P.crit = y > threshold

# Now shows the "dynamic traceplot" as BUGS would do
# These are the simulated values in all the steps of the simulation
plot(y,t="l",xlim=c(10500,11000))  # NB: the R command 'xlim=c(xx,yy)'
                                   # instructs R to define the limits for the 
                                   # x-axis to be between xx and yy
plot(P.crit,t="l",xlim=c(10500,11000))

# First creates a matrix by stacking together the two columns (one for P8
# and the other for y), using the 'cbind' R built-in command (type '?cbind'
# to open the help window)
sims = cbind(P.crit,y,theta)
# Then uses the function 'stats' (which is made available by 'utils.R') to show
# summary statistics like BUGS would do
stats(sims)

# Now uses R plots to show graphical summaries of the distributions
hist(as.numeric(P.crit))  # note that R defines 'P8' as a "logical" variable
                          # ie it takes values TRUE/FALSE. In order to plot its 
                          # distribution, we need to turn it into a "numeric"
                          # variable, which we do by using the R built-in 
                          # 'as.numeric' function

hist(y)   # We don't need any transformation for y, because it is
          # already a numeric variable

plot(density(theta))   # NB: uses the built-in R function 'density(x)' to plot
                       # a "kernel" approximation to the underlying continuous
                       # distribution of theta, derived by the simulations

# Now you can modify the prior for theta implying a Uniform(0,1) distribution
theta = runif(nsim, 0, 1)
plot(
  seq(0,1,.001),dunif(seq(0,1,.001),0,1),t="l",xlab="theta",ylab="p(theta)",
  main="Uniform prior"
)

# Now simulate the values for the trial
npats = 20 # number of patients in the trial
y = rbinom(nsim,npats,theta)
# Computes the probability of 15 or more successes in each simulation
threshold = 14.5 # y is only integers, so if the threshold is 13.5, it means
# there are 15 or more successes in a single simulation
P.crit = y > threshold
# Summarises the new results
sims2 = cbind(P.crit,y,theta)
stats(sims2)
# And compares with the old results
stats(sims)
