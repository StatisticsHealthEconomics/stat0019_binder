# Loads the relevant packages
library(survHE)

# Loads & visualises the data
load(here::here("07_survival","survival_data.Rdata"))
head(dat)

# Defines the model formula, including the treatment arm as the only covariate
formula=Surv(time,event)~as.factor(arm)

# Defines a vector of strings identifying the models used (Exponential and Weibull)
mods=c("exp","weibull")

# Runs survHE to estimate the model parameters
m1=fit.models(formula=formula,data=dat,distr=mods)

# Explores the model output to check what has been stored
names(m1)
lapply(m1,names)
# Uses the 'print' and 'plot' methods (specialised for survHE) to summarise the output
print(m1)     # First print the summary for the Exponential model (the first in 'mods')
print(m1,2)   # Then the output for the Weibull model (the second element of 'mods')
plot(m1)

# Runs survHE to estimate the two models using INLA
m2=fit.models(formula=formula,data=dat,distr=mods,method="inla")
print(m2)     # First print the summary for the Exponential model (the first in 'mods')
print(m2,2)   # Then the output for the Weibull model (the second element of 'mods')

## Runs survHE to estimate the two models using INLA
## **NB**: The Binder VM breaks if you try to install the 'hmc' module.
##         This is because 'survHEhmc' pre-compiles all the HMC/rstan
##         models, which makes them very quick to run. But in doing so, it
##         creates a very large package (that's why the three modules have been
##         separated into three separate packages: 'survHE', 'survHEinla' and
##         'survHEhmc') and it crashes the Binder memory. 
## If you want to do this on your computer, you can use the following commands
## (simply uncomment them)
#m3=fit.models(formula=formula,data=dat,distr=mods,method="hmc")

## NB: the model has been run on a more powerful (local) machine and the output is
##     saved in the object 'm3_HMC.rds', so you can still use it in the Binder too!
m3=readRDS(here::here("07_survival","m3_HMC.rds"))
print(m3)     # First print the summary for the Exponential model (the first in 'mods')
print(m3,2)   # Then the output for the Weibull model (the second element of 'mods')

# Plots all the results in the same graph
plot(
  # Defines the objects and the relevant models to plot
  m1,m2,m3,mod=c(1,2,3,4,5,6),
  # Selects a vector of colours & labels for each
  # More details on the underlying survHE options can be found at:
  # https://gianluca.statistica.it/software/survhe/tutorial.html
  colour=c("blue","green","red","yellow","magenta","orange"),
  lab.model=c("Exponential (MLE)", "Weibull (MLE)","Exponential (INLA)", 
         "Weibull (INLA)","Exponential (HMC)", "Weibull (HMC)"
  ),
  # Selects the time frame for the plot (extrapolation beyond data)
  t=seq(0,50)
)

# Perform PSA on the survival curves (induced by uncertainty in the model parameters)
psa1=make.surv(m1,mod=2,t=seq(0.01,50),nsim=1000)
psa.plot(psa1,col=c("blue","red"))

psa2=make.surv(m2,mod=2,t=seq(0.01,50),nsim=1000)
psa.plot(psa2,col=c("blue","red"))

psa3=make.surv(m3,mod=2,t=seq(0.01,50),nsim=1000)
psa.plot(psa3,col=c("blue","red"))

