library(BCEA)
library(mgcv)

setwd(here::here("12_evppi"))

## Load the data from our model
psa <- read.csv("psa_model_1.csv")
head(psa)
dim(psa)

## Calculate the net benefit with a willingness to pay of 10000
K <- 10000
nb1 <- psa$effects1 * K - psa$costs1
nb2 <- psa$effects2 * K - psa$costs2
inb <- nb2 - nb1

## Investigate the mean net benefit and incremental net benefit
mean(nb1)
mean(nb2)
mean(inb)

## Calculate EVPI for this model
# You can run ?pmax for more information on the pmax function
EVPI <- mean(pmax(nb1, nb2)) - max(mean(nb1), mean(nb2)) 
EVPI

# equivalently, we can calculate EVPI using the incremental net benefit. 
# This is useful later when we use regression to estimate EVPPI.
EVPI <- mean(pmax(0, inb)) - max(0, mean(inb))
EVPI

## Calculate EVPPI using regression
theta <- psa[,1:31][-c(7, 16)]
model <- gam(inb ~ s(theta[, 5]))

# Check the model fit
plot(fitted(model), residuals(model))
gam.check(model)

## EVPPI for each parameter
# initialise a vector to hold the results. 
evppi_vector <- numeric(29)

for (i in (1:29)) { # theta 7 and 16 are not defined
  print(i)
  parameter_of_interest <- theta[, i]
  model <- gam(inb ~ s(parameter_of_interest))
  fittedValues <- fitted(model)
  evppi_vector[i] <- mean(pmax(0, fittedValues)) - max(0, mean(fittedValues))
}

names(evppi_vector) <- paste("theta", 1:29)
print(round(evppi_vector, 1))

# plot
x.points <- barplot(
  evppi_vector, ylab = "partial EVPI", xaxt = "n", ylim = c(0, 500)
)
axis(
  side = 1, at = x.points, labels = names(evppi_vector), 
  tick = FALSE, hadj = 0.8, las = 3
)

# In BCEA
effects <- cbind(psa$effects1, psa$effects2)
costs <- cbind(psa$costs1, psa$costs2)
bcea.out <- bcea(eff = effects, cost = costs, ref = 2,
                 plot = F)
inputs <- createInputs(theta)
info.rank(bcea.out, inputs, wtp = K)

## Mutli-parameter EVPPI
model <- gam(inb ~ te(theta[, 5], theta[, 14]))
fittedValues <- fitted(model)
mean(pmax(0, fittedValues)) - max(0, mean(fittedValues))

# 4 parameters
system.time(
  model <- gam(inb ~ te(theta[, 5], theta[, 6], theta[, 14], theta[, 15]))
)
fittedValues <- fitted(model)
mean(pmax(0, fittedValues)) - max(0, mean(fittedValues))
gam.check(model)

### you can also use the earth package to speed up computation but this is 
### not available on the binder
# install.packages("earth")
# library(earth)

## EVPPI with multiple decision options
rm(theta, nb1, nb2, inb, model, psa)

# Load data
psa <- read.csv("psa_model_2.csv")
head(psa)
dim(psa)

# Extract the costs, effects and parameters
costs <- psa[, grep("costs", names(psa))] # grep matches text strings
effects <- psa[, grep("effects", names(psa))]
theta <- psa[, grep("theta", names(psa))]

# Willingness to pay
K <- 10000 # set K to 10000 again

# Net benefit
nb <- effects * K - costs
inb <- nb - nb[, 1]

colMeans(nb)
colMeans(inb)

# EVPI
inb <- as.data.frame(inb) # convert to data.frame so do.call works
mean(do.call(pmax, inb)) - max(colMeans(inb))

## Single parameter EVPPI
model_2v1 <- gam(inb[, 2] ~ s(theta[, 5]))
model_3v1 <- gam(inb[, 3] ~ s(theta[, 5]))

fittedValues <- data.frame(
  fittedValue_1v1 = 0, # note that the incremental NB of option 1 vs 1 is zero
  fittedValue_2v1 = fitted(model_2v1), 
  fittedValue_3v1 = fitted(model_3v1)
)

mean(do.call(pmax, fittedValues)) - max(colMeans(fittedValues))

## Try calculating for all the single parameters

## EVPPI for groups
model_2v1 <- gam(inb[, 2] ~ te(theta[, 5], theta[, 6], theta[, 7]))
model_3v1 <- gam(inb[, 3] ~ te(theta[, 5], theta[, 6], theta[, 7]))

fittedValues <- data.frame(
  fittedValue_1v1 = 0, # note that the incremental NB of option 1 vs 1 is zero
  fittedValue_2v1 = fitted(model_2v1), 
  fittedValue_3v1 = fitted(model_3v1)
)

mean(do.call(pmax, fittedValues)) - max(colMeans(fittedValues))

## NOTE earth is not available on the binder.

## Using BCEA
# Create the BCEA object
bceaObject <- bcea(
  e = as.matrix(effects), c = as.matrix(costs), ref = 1, 
  wtp = 10000, plot = FALSE
)

# Calculate EVPPI
evppiObject <- evppi(bceaObject, 5, theta)

# Residual checking
diag.evppi(evppiObject, bceaObject, int = 1)
diag.evppi(evppiObject, bceaObject, int = 2)
# Print EVPPI
evppiObject$evppi


## Multi parameter EVPPI
bceaObject <- bcea(
  e = as.matrix(effects[1:1000, ]), c = as.matrix(costs[1:1000, ]), 
  ref = 1, wtp = 10000, plot = FALSE
)

evppiObject <- evppi(bceaObject, c(5, 6, 7), theta[1:1000, ], method = "GAM")
diag.evppi(evppiObject, bceaObject, int = 1)
diag.evppi(evppiObject, bceaObject, int = 2)
# The residuals here have structure, meaning that the GAM model may not be appropriate 
# for this EVPPI calculation.
evppiObject$evppi

### Launch SAVI
library(remotes) # If you are running this not on the binder, you may need to install this package
install_github('Sheffield-Accelerated-VoI/SAVI-package')
library(SAVI)
SAVI()
