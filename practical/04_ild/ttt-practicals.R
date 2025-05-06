# Loads the 'tidyverse' library to use 'ggplot'
library(tidyverse)

# Loads the data & discards records with missing values
# Loads the dataset
ttt=readRDS(here::here("04_ild/ttt.rds"))

# Visualise the data
ttt |> ggplot(aes(qaly)) + geom_histogram(fill="grey",col="black") + 
  facet_grid(.~as.factor(arm)) + xlab("QALYs") + ylab("")
ttt |> ggplot(aes(totalcost)) + geom_histogram(fill="grey",col="black") + 
  facet_grid(.~as.factor(arm)) + xlab("Total costs (GBP)") + ylab("")

# Define the relevant variables
e=ttt$qaly
c=ttt$totalcost
Trt=ttt$Trt
u0star=scale(ttt$qol_0,scale=F) |> as.numeric()
N=ttt |> nrow()

# Creates the datalist
data=list(e=e,c=c/1000,Trt=Trt,u0star=u0star,N=N)

# Initialises the object 'model' as a list to store all the results
model=list()

# Runs JAGS in the background
library(R2jags)

# Normal/Normal independence
nn_indep=function(){
  for (i in 1:N) {
    # Model for the effects
    e[i] ~ dnorm(phi.e[i], tau.e[Trt[i]])
    phi.e[i] <- alpha0 + alpha1*(Trt[i]-1) + alpha2*u0star[i]
    # Model for the costs
    c[i] ~ dnorm(phi.c[i], tau.c[Trt[i]])
    phi.c[i] <- beta0 + beta1*(Trt[i]-1)
  }
  # Rescales the main economic parameters
  # NB: Cost must be rescaled (x1000)
  for (t in 1:2) {
    mu.e[t] <- alpha0 + alpha1*(t-1)
    mu.c[t] <- 1000*(beta0 + beta1*(t-1))
  }
  # Minimally informative priors on the regression coefficients
  alpha0 ~ dnorm(0,0.0001)
  alpha1 ~ dnorm(0,0.0001)
  alpha2 ~ dnorm(0,0.0001)
  beta0 ~ dnorm(0,0.0001)
  beta1 ~ dnorm(0,0.0001)
  for (t in 1:2) {
    # PC-prior on the sd with Pr(sigma_e>0.8) \approx 0.01
    sigma.e[t] ~ dexp(5.75)
    tau.e[t] <- pow(sigma.e[t],-2)
    # PC-prior on the sd with Pr(sigma_c>2) \approx 0.5
    sigma.c[t] ~ dexp(0.35)
    tau.c[t] <- pow(sigma.c[t],-2)
  }
}

# Visualises the PC priors for the model coefficients
ggplot() + stat_function(fun=dnorm,args=list(mean=0,sd=100)) + xlim(-500,500) +
  xlab("") + ylab("")
ggplot() + stat_function(fun=dexp,args=list(rate=5.75)) + xlim(0,2) +
  xlab("") + ylab("")
ggplot() + stat_function(fun=dexp,args=list(rate=0.35)) + xlim(0,10) +
  xlab("£ x 1,000") + ylab("")


# Stores the ooutput as an element named 'nn_indep' in the object 'model'
model$nn_indep=jags(
  data=data,
  parameters.to.save=c(
    "mu.e","mu.c","alpha0","alpha1","alpha2","beta0","beta1","sigma.e","sigma.c"
  ),
  inits=NULL,n.chains=2,n.iter=5000,n.burnin=3000,n.thin=1,DIC=TRUE,
  # This specifies the model code as the function 'nn_indep'
  model.file=nn_indep
)
print(model$nn_indep,digits=3,interval=c(0.025,0.5,0.975))


# Normal/Normal MCF
nn_mcf=function(){
  for (i in 1:N) {
    # Marginal model for the effects
    e[i] ~ dnorm(phi.e[i], tau.e[Trt[i]])
    phi.e[i] <- alpha0 + alpha1*(Trt[i]-1) + alpha2*u0star[i]
    # *Conditional* model for the costs
    c[i] ~ dnorm(phi.c[i], tau.c[Trt[i]])
    phi.c[i] <- beta0 + beta1*(Trt[i]-1) + beta2*(e[i]-mu.e[Trt[i]]) 
  }
  # Rescales the main economic parameters
  for (t in 1:2) {
    mu.e[t] <- alpha0 + alpha1*(t-1)
    mu.c[t] <- 1000*(beta0 + beta1*(t-1))
  }
  # Minimally informative priors on the regression coefficients
  alpha0 ~ dnorm(0,0.0001)
  alpha1 ~ dnorm(0,0.0001)
  alpha2 ~ dnorm(0,0.0001)
  beta0 ~ dnorm(0,0.0001)
  beta1 ~ dnorm(0,0.0001)
  beta2 ~ dnorm(0,0.0001)
  for (t in 1:2) {
    # PC-prior on the *marginal* sd with Pr(sigma_e>0.8) \approx 0.01
    sigma.e[t] ~ dexp(5.75)
    tau.e[t] <- pow(sigma.e[t],-2)
    # PC-prior on the *conditional* sd with Pr(lambda_c>2) \approx 0.5
    lambda.c[t] ~ dexp(0.34)
    tau.c[t] <- pow(lambda.c[t],-2)
    # Retrieves the correlation coefficients
    rho[t] <- beta2*sigma.e[t]/sigma.c[t]
    # And the *marginal* standard deviation for the cost
    sigma.c[t] <- sqrt(pow(lambda.c[t],2) + pow(sigma.e[t],2)*pow(beta2,2))
  }
}

model$nn_mcf=jags(
  data=data,
  parameters.to.save=c(
    "mu.e","mu.c","alpha0","alpha1","alpha2","beta0","beta1","beta2",
    "sigma.e","sigma.c","lambda.c","rho"
  ),
  inits=NULL,n.chains=2,n.iter=5000,n.burnin=3000,n.thin=1,DIC=TRUE,
  # This specifies the model code as the function 'nn_mcf'
  model.file=nn_mcf
)

# Shows the summary statistics from the posterior distributions. 
print(model$nn_mcf,digits=3,interval=c(0.025,0.5,0.975))


# Gamma/Gamma MCF
gg_mcf=function(){
  for (i in 1:N) {
    # Marginal model for the *rescaled* effects
    estar[i] ~ dgamma(nu.e[Trt[i]], gamma.e[Trt[i],i])
    gamma.e[Trt[i],i] <- nu.e[Trt[i]]/phi.e[i]
    log(phi.e[i]) <- alpha0 + alpha1*(Trt[i]-1) + alpha2*u0star[i]
    # Conditional model for the costs
    c[i] ~ dgamma(nu.c[Trt[i]], gamma.c[Trt[i],i])
    gamma.c[Trt[i],i] <- nu.c[Trt[i]]/phi.c[i]
    log(phi.c[i]) <- beta0 + beta1*(Trt[i]-1) + beta2*(estar[i]-mustar.e[Trt[i]]) 
  }
  # Rescales the main economic parameters
  for (t in 1:2) {
    mustar.e[t] <- exp(alpha0 + alpha1*(t-1))
    mu.e[t] <- 3 - mustar.e[t]
    mu.c[t] <- 1000*exp(beta0 + beta1*(t-1))
  }
  # Minimally informative priors on the regression coefficients
  alpha0 ~ dnorm(0,0.0001)
  alpha1 ~ dnorm(0,0.0001)
  alpha2 ~ dnorm(0,0.0001)
  beta0 ~ dnorm(0,0.0001)
  beta1 ~ dnorm(0,0.0001)
  beta2 ~ dnorm(0,0.0001)
  # PC prior on the shape parameters 
  # assume that Pr(nu.e>30) = Pr(nu.c>30) \approx 0.01
  for (t in 1:2) {
    nu.e[t] ~ dexp(0.15)
    nu.c[t] ~ dexp(0.15)
  }
}

# Adds the variable 'estar' to the datalist
data$estar=3-data$e

# Runs JAGS in the background and stores the output in the element 'gg_mcf' 
model$gg_mcf=jags(
  data=data,
  parameters.to.save=c(
    "mu.e","mustar.e","mu.c","alpha0","alpha1","alpha2","beta0",
    "beta1","beta2","nu.e","nu.c"
  ),
  inits=NULL,n.chains=2,n.iter=5000, n.burnin=3000,n.thin=1,DIC=TRUE,
  # This specifies the model code as the function 'model.code'
  model.file=gg_mcf
)

# Shows the summary statistics from the posterior distributions. 
print(model$gg_mcf,digits=3,interval=c(0.025,0.5,0.975))


# Cost-effectiveness analysis
library(BCEA)
# Defines a list of labels for the interventions
interventions=c("Standard","Active intervention")
# Sets the reference (t=2=intensive case management)
ref=2

# 1. Normal/Normal independent
# Defines the variables for effects and costs
eff=model$nn_indep$BUGSoutput$sims.list$mu.e
cost=model$nn_indep$BUGSoutput$sims.list$mu.c
# Runs BCEA  
m_nn_indep=bcea(
  eff=eff,cost=cost,interventions=interventions,ref=ref
)

# 2. Normal/Normal MCF
# Defines the variables for effects and costs
eff=model$nn_mcf$BUGSoutput$sims.list$mu.e
cost=model$nn_mcf$BUGSoutput$sims.list$mu.c
# Runs BCEA
m_nn_mcf=bcea(
  eff=eff,cost=cost,interventions=interventions,ref=ref
)

# 3. Gamma/Gamma MCF
# Defines the variables for effects and costs
eff=model$gg_mcf$BUGSoutput$sims.list$mu.e
cost=model$gg_mcf$BUGSoutput$sims.list$mu.c
# Runs BCEA
m_gg_mcf=bcea(
  eff=eff,cost=cost,interventions=interventions,ref=ref
)

# Plots the results
contour2(m_nn_indep,graph="gg")
contour2(m_nn_mcf,graph="gg")
contour2(m_gg_mcf,graph="gg")

p1=ceac.plot(m_nn_indep,graph="gg")
p2=ceac.plot(m_nn_mcf,graph="gg")
p3=ceac.plot(m_gg_mcf,graph="gg")

# Post-processes the graphs to superimpose the three CEACs
# First stack the underlying data (extracted from the three 'ggplot' objects)
# And associate each with a 'model' variable
p1$data |> mutate(model="Normal/Normal indep") |>
  bind_rows(p2$data |> mutate(model="Normal/Normal MCF")) |>
  bind_rows(p3$data |> mutate(model="Gamma/Gamma MCF")) |> 
  # Then use 'ggplot' to make the actual plot using the three original datasets
  # but different colours for each model
  ggplot(aes(k,ceac,color=model)) + geom_line() +
  # And then provide further customisation
  xlab("Willingness to pay") + ylim(0,1) +
  ylab("Probability of cost effectiveness") + 
  labs(title="Cost Effectiveness Acceptability Curves",color="Model") +
  theme(
    legend.position="inside",legend.position.inside=c(.6,.75),
    legend.background=element_blank()
  ) 

