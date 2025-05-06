#' Influenza example --- source: Cooper et al (2004); Baio (2012)
#' 
#' Defines the data
#' Number of interventions 
#'   t=0: control; 
#'   t=1: prophylactic use of Neuramidase Inhibitors (NI) 
T <- 2 					

#' Evidence synthesis on effectiveness of NIs prophylaxis vs placebo
r0 <- r1 <- n0 <- n1 <- numeric()	# defines observed cases & sample sizes
r0 <- c(34,40,9,19,6,34)
r1 <- c(11,7,3,3,3,4)
n0 <- c(554,423,144,268,251,462)
n1 <- c(553,414,144,268,252,493)
S <- length(r0)				# number of relevant studies

#' Evidence synthesis on incidence of influenza in healthy adults (under t=0)
y <- m <- numeric()			# defines observed values for baseline risk
y <- c(0,6,5,6,25,18,14,3,27)
m <- c(23,241,159,137,519,298,137,24,132)
H <- length(y)

#' Data on costs
unit.cost.drug <- 2.4		# unit (daily) cost of NI
length.treat <- 6*7		# 6 weeks course of treatment
c.gp <- 19			# cost of GP visit to administer prophylactic NI
vat <- 1.175			# VAT
c.ni <- unit.cost.drug*length.treat*vat 

# Informative prior on cost of influenza 
mu.inf <- 16.78			# mean cost of influenza episode
sigma.inf <- 2.34		# sd cost of influenza episode
tau.inf <- 1/sigma.inf^2	# precision cost of influenza episode

# Informative prior on length of influenza episodes
# Compute the value of parameters (mulog,sigmalog) for a logNormal 
# distribution to have mean and sd (m,s) - check help(lognPar)
m.l=bmhe::lognPar(8.2,2)$mulog          # original value in the paper: 8.2
s.l=bmhe::lognPar(8.2,2)$sigmalog       # original value in the paper: sqrt(2)

#' Prepares to launch `JAGS`
library(R2jags)

#' Creates the data list
data <- list(
  S=S,H=H,r0=r0,r1=r1,n0=n0,n1=n1,y=y,m=m
)

#' Points to the txt file where the OpenBUGS model is saved
filein <- here::here("05_ald","EvSynth.txt")

#' Defines the parameters list
params <- c("p0","p1","or","l","c.inf")

# Creates a function to draw random initial values 
inits=function(){
  list(alpha=rnorm(data$S,0,.5),delta=rnorm(data$S,0,.5))
}

#' Sets the number of iterations, burnin and thinning
n.iter <- 10000
n.burnin <- 9500
n.thin <- 2

#' Finally calls JAGS to do the MCMC run and saves results to the object "es"
es=jags(
  data=data,
  parameters.to.save=params,
  inits=inits,n.chains=2,n.iter=25000,n.burnin=5000,n.thin=10,DIC=TRUE,pD=TRUE,
  model.file=filein
)

#' Displays the summary statistics
print(es,digits=3,intervals=c(0.025, 0.975))

#' Convergence check through traceplots (example for node p1)
bmhe::traceplot(es,"p1")

################################################################################
#' NB: If you want/need to run `BUGS` and `R2OpenBUGS`, remember to 
#' 1. Run `library(R2OpenBUGS)` to load the relevant package to connect 
#'    `R` to `OpenBUGS`
#' 2. Run `bugs(...)` instead of `jags(...)`
#' 3. `R2OpenBUGS` stores the MCMC output under the root sub-object 
#'    so there is no `$BUGSoutput` sub-object. For instance, you need to run 
#'    `plot(es$sims.list$p1[1:500],t="l",col="blue",ylab="p1")`
#'    to make a traceplot of the first 500 simulations for p1
#'    But the functions in `bmhe` will know how to do this automatically!
################################################################################


#' Runs economic analysis 
#' Define the variables 'mu.e' and 'mu.c' as empty matrices with 'n.sims' 
#' rows (one per each MCMC simulation) and 2 columns (one per treatment arm)
n.sims=es$BUGSoutput$n.sims
mu.e=mu.c=matrix(NA,n.sims,2)

# Fixed costs
# Cost of GP visit to issue NIs prescription
c.gp=19
# Cost of course of treatment with NIs
c.ni=2.40*6*7*1.175

# Extracts simulations from the object 'flu' (for simplicity in the code)
p0=es$BUGSoutput$sims.list$p0
p1=es$BUGSoutput$sims.list$p1
c.inf=es$BUGSoutput$sims.list$c.inf
l=es$BUGSoutput$sims.list$l

# Population average costs
mu.c[,1]=(1-p0)*c.gp + p0*(c.gp+c.inf)
mu.c[,2]=(1-p1)*(c.gp+c.ni) + p1*(c.gp+c.ni+c.inf)

# Population average effects
mu.e[,1]=-l*p0
mu.e[,2]=-l*p1

# Now uses 'BCEA' to make the economic analysis
library(BCEA)
m=bcea(mu.e,mu.c,interventions=c("Status quo","NIs prophylaxis"),ref=2)
# Economic analysis
summary(m,wtp=1000)
