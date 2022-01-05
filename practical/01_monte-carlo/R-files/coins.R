# This files essentially does the same analysis as OpenBUGS

# Set up
n = 10 # number of "trials"
theta = 0.5 # probability of "success" (heads) - the coin is "balanced"

# Simulates n tosses of a coin
nsim = 1000 # number of simulations
y = rbinom(nsim,n,theta)  # NB: the built-in R command 'rxxxx' (eg 'rbinom')
                          # is used to simulate random values ('r') from a
                          # given distribution ('binom'=Binomial). Type 
                          # '?Distributions' (without ''!) in your R terminal
                          # to get some help

# Computes the probability of 8 or more heads in each simulation
threshold = 7.5 # y is only integers, so if the threshold is 7.5, it means
                # there are 8 or more heads in a single simulation
P8 = y > threshold

# Now summarises the results of the simulations using a specific R function
# which is defined in the utility file 'utils.R'. This is "sourced" (eg loaded
# up on the R workspace) using the R command 'source' (type '?source' for help)
source("utils.R")

# First creates a matrix by stacking together the two columns (one for P8
# and the other for y), using the 'cbind' R built-in command (type '?cbind'
# to open the help window)
sims = cbind(P8,y)
# Then uses the function 'stats' (which is made available by 'utils.R') to show
# summary statistics like BUGS would do
stats(sims)

# Now uses R plots to show graphical summaries of the distributions
hist(as.numeric(P8))  # note that R defines 'P8' as a "logical" variable
                      # ie it takes values TRUE/FALSE. In order to plot its 
                      # distribution, we need to turn it into a "numeric"
                      # variable, which we do by using the R built-in 
                      # 'as.numeric' function
hist(y)               # We don't need any transformation for y, because it is
                      # already a numeric variable
