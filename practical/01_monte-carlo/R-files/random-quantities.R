# This files essentially does the same analysis as OpenBUGS

# First we need to simulate values for a variable y ~ Normal(0,1)
# We can easily do this using the built-in R command 'rnorm'
nsims = 10000
y = rnorm(nsims, mean=0, sd=1)

# Now we can plot the resulting density using the simulations
plot(density(y))
# Or even the analytic distribution using R built-in commands
plot(
  seq(-4,4,.01),dnorm(seq(-4,4,.01),0,1),t="l",xlab="y",ylab="p(y)"
)
# NB: The analytic plot is a lot smoother! 

# To simulate a different Normal, just change the parameters...
y = rnorm(nsims, mean=1, sd=2)

# And to construct a deterministic function of that random variable, simply 
# do some algebra...
z = y^3
mean(z)           # mean of Z
var(z)            # variance of Z
sum(z>10)/nsims   # Computes the proportion of simulations where z > 10
                  # This estimates Pr(Z>10)
# Plots the density of z
plot(density(z))
