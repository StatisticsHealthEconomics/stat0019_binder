# This files uses Gibbs sampling to simulate from the posterior distribution for
# the parameters of a semi-conjugated mode, as shown in the lecture

# Can also change working directory with
# here::here("02_mcmc")
# but it's not strictly necessary here, because we're not using any external
# file and all the variables are defined below, locally...

# Vector of observed data
y = c(1.2697,7.7637,2.2532,3.4557,4.1776,6.4320,-3.6623,7.7567,5.9032,
7.2671,-2.3447,8.0160,3.5013,2.8495,0.6467,3.2371,5.8573,-3.3749,
4.1507,4.3092,11.7327,2.6174,9.4942,-2.7639,-1.5859,3.6986,2.4544,
-0.3294,0.2329,5.2846)

# "Hyper-parameters" (ie parameters for the prior distributions)
mu_0 = 0
sigma2_0 = 10000
alpha_0 = 0.01
beta_0 = 0.01

# Sample size and sample mean of the data
n = length(y)
ybar = mean(y)

# Sets the "seed" (for reproducibility). With this command, you will
# *always* get the exact same output
set.seed(13)
# Initialises the parameters
mu = tau = numeric()
sigma2 = 1/tau
mu[1] = rnorm(1,0,3)
tau[1] = runif(1,0,3)
sigma2[1] = 1/tau[1]


# Sets the number of iterations (nsim)
nsim = 1000
# Loops over to sequentially update the parameters
for (i in 2:nsim) {
  # 1. Updates the sd of the full conditional for mu
  sigma_1 = sqrt(1/(1/sigma2_0 + n/sigma2[i-1]))
  # 2. Updates the mean of the full conditional for mu
  mu_1 = (mu_0/sigma2_0 + n*ybar/sigma2[i-1])*sigma_1^2
  # 3. Samples from the updated full conditional for mu
  mu[i] = rnorm(1,mu_1,sigma_1)
  
  # 4. Updates the 1st parameter of the full conditional for tau
  alpha_1 = alpha_0+n/2
  # 5. Updates the 2nd parameter of the full conditional for tau
  beta_1 = beta_0 + sum((y-mu[i])^2)/2
  # 6. Samples from the updated full conditional for tau
  tau[i] = rgamma(1,alpha_1,beta_1)
  # 7. Re-scales the sampled value on the variance scale
  sigma2[i] = 1/tau[i]
}

sigma = sqrt(sigma2)

# Histograms from the posterior distributions
hist(mu)
hist(sigma)

# Traceplots
plot(mu,t="l",bty="l")
plot(sigma,t="l",bty="l")

