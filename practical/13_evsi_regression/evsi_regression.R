## Load Packages
library(mgcv)

## Set the working directory
setwd(here::here("13_evsi_regression"))
# This needs to be changed to your working directory on your own computer

## Load the data
psa <- read.csv("psa_model_3.csv")
head(psa)
dim(psa)

# Extract costs, effects and the parameters
costs <- psa[, grep("costs", names(psa))]
effects <- psa[, grep("effects", names(psa))]
theta <- psa[, grep("theta", names(psa))]

# Set the willingness to pay
K <- 10000 # set K to 10000 

# Calculate net monetary benefit and incremental net benefit
nb <- effects * K - costs
inb <- nb[, 2] - nb[, 1]
colMeans(nb)
mean(inb)

## EVPI
nb <- as.data.frame(nb) # convert to data.frame so do.call works
mean(do.call(pmax, nb)) - max(colMeans(nb))

## Calculate single parameter EVPPI
# initialise a vector to hold the results. 
evppi_vector <- numeric(19)

for (i in 1:19) {
  print(i)
  parameter_of_interest <- theta[, i]
  model <- gam(inb ~ s(parameter_of_interest))
  fittedValues <- data.frame(0, fitted(model))
  evppi_vector[i] <- mean(do.call(pmax, fittedValues)) - max(colMeans(fittedValues))
}

names(evppi_vector) <- paste("theta", 1:19)
print(round(evppi_vector, 1))

### EVSI calculation
# create an object called theta14 that holds the 14th column of theta
theta14 <- theta[, 14]
theta14[1:10]
summary(theta14)

# Plot theta14
hist(theta14)

## Simulate data
# For a single dataset
set.seed(1)
N <- 100
x <- rbinom(1, N, theta14[1])
x
# summary statistic
summary.statistics <- x / N
summary.statistics

# For multiple datasets
set.seed(1)
N <- 100
x <- rbinom(10, N, theta14[1:10]) # Example for the first 10 values of theta14
summary.statistics <- x / N
summary.statistics

N <- 100
x <- rbinom(length(theta14), N, theta14) # For all the parameters
summary.statistics <- x / N

# Fit the regression model
model <- gam(inb ~ s(summary.statistics))
fittedValues <- data.frame(0, fitted(model))
evsi <- mean(do.call(pmax, fittedValues)) - max(colMeans(fittedValues))
evsi

# Check model
gam.check(model)

## Range of study sizes
set.seed(1)

evsi_values <- c() # initialise empty vector for results
study_sizes <- c(10, 20, 50, 100, 200, 500, 1000)
for (N in study_sizes) {
  print(N)
  x <- rbinom(length(theta14), N, theta14)
  summary.statistics <- x / N
  model <- gam(inb ~ s(summary.statistics))
  fittedValues <- data.frame(0, fitted(model))
  evsi <- mean(do.call(pmax, fittedValues)) - max(colMeans(fittedValues))
  evsi_values <- c(evsi_values, evsi)
}
data.frame(study_sizes, evsi_values)
plot(study_sizes, evsi_values, type = "b")

# Compare to EVPPI
evppi_vector[14]

## EVSI for theta 5
set.seed(1)
theta5 <- theta[, 5]

evsi_values <- NULL
study_sizes <- c(10, 20, 50, 100, 200, 500, 1000)
for (N in study_sizes) {
  print(N)
  x <- rbinom(length(theta5), N, theta5)
  summary.statistics <- x / N
  model <- gam(inb ~ s(summary.statistics))
  fittedValues <- data.frame(0, fitted(model))
  evsi <- mean(do.call(pmax, fittedValues)) - max(colMeans(fittedValues))
  evsi_values <- c(evsi_values, evsi)
}
data.frame(study_sizes, evsi_values)
plot(study_sizes, evsi_values, type = "b")

# Compare to EVPPI
evppi_vector[5]

## Ensuring monotonicity by replicating the data
theta5_rep <- rep(theta5, 100)
inb_rep <- rep(inb, 100)

# Calculate EVSI across sample size
set.seed(1)

evsi_values <- NULL
study_sizes <- c(10, 20, 50, 100, 200, 500, 1000)
for (N in study_sizes) {
  print(N)
  x <- rbinom(length(theta5_rep), N, theta5_rep)
  summary.statistics <- x / N
  model <- gam(inb_rep ~ s(summary.statistics))
  fittedValues <- data.frame(0, fitted(model))
  evsi <- mean(do.call(pmax, fittedValues)) - max(colMeans(fittedValues))
  evsi_values <- c(evsi_values, evsi)
}
data.frame(study_sizes, evsi_values)
plot(study_sizes, evsi_values, type = "b")

## Simulating an RCT
set.seed(1)

N_arm1 <- 100
N_arm2 <- 100

x_arm1 <- rbinom(length(theta5), N_arm1, theta5)
x_arm2 <- rbinom(length(theta14), N_arm2, theta14)

# Summarising as an RCT - only informs relative effect, i.e., log odds ratio
N_arm1_adj <- ifelse(x_arm1 == 0 | x_arm1 == N_arm1, N_arm1 + 1, N_arm1)
N_arm2_adj <- ifelse(x_arm2 == 0 | x_arm2 == N_arm2, N_arm2 + 1, N_arm2)
x_arm1_adj <- ifelse(x_arm1 == 0 | x_arm1 == N_arm1, x_arm1 + 0.5, x_arm1)
x_arm2_adj <- ifelse(x_arm2 == 0 | x_arm2 == N_arm2, x_arm2 + 0.5, x_arm2)

p_arm1_adj <- x_arm1_adj / N_arm1_adj
p_arm2_adj <- x_arm2_adj / N_arm2_adj
odds_ratio <- (p_arm1_adj / (1 - p_arm1_adj) ) / (p_arm2_adj / (1 - p_arm2_adj))
l_OR <- log(odds_ratio)

# Calculate EVSI
model_OR <- gam(inb ~ s(l_OR))
fittedValues <- data.frame(0, fitted(model_OR))
evsi <- mean(do.call(pmax, fittedValues)) - max(colMeans(fittedValues))
evsi
# Model checking
gam.check(model_OR)

# EVSI informing both absolute effects
model_p1_p2 <- gam(inb ~ te(p_arm1_adj, p_arm2_adj))
fittedValues <- data.frame(0, fitted(model_p1_p2))
evsi <- mean(do.call(pmax, fittedValues)) - max(colMeans(fittedValues))
evsi # EVSI is higher

gam.check(model_p1_p2)

# EVSI for absolute effects and log odds ratio
model_p1_p2_OR <- gam(inb ~ te(p_arm1_adj, p_arm2_adj, l_OR))
fittedValues <- data.frame(0, fitted(model_p1_p2_OR))
evsi <- mean(do.call(pmax, fittedValues)) - max(colMeans(fittedValues))
evsi # No extra information

gam.check(model_p1_p2_OR)

## EVSI for time to event data
# Simulate data
theta16 <- theta[, 16] # theta16 is the mean survivial
set.seed(1)

shape <- 0.8 # Variability specified separately
scale <- theta16[1] / gamma(1 + 1 / shape)
duration <- rweibull(10, shape, scale)
duration

# Summarise the data with median survival time
x <- numeric(length(theta16))
sample_size <- 20
set.seed(1)

for (i in 1:length(theta16)) {
  shape <- 0.8
  scale <- theta16[i] / gamma(1 + 1 / shape)
  duration <- rweibull(sample_size, shape, scale)
  x[i] <- median(duration) 
}

# Calculate EVSI
model <- gam(inb ~ s(x))
fittedValues <- data.frame(0, fitted(model))
evsi <- mean(do.call(pmax, fittedValues)) - max(colMeans(fittedValues))
evsi
# Model checking
gam.check(model)
evppi_vector[16]

# Change the sample size
x <- numeric(length(theta16))
sample_size <- 50
set.seed(1)

for (i in 1:length(theta16)) {
  shape <- 0.8
  scale <- theta16[i] / gamma(1 + 1 / shape)
  duration <- rweibull(sample_size, shape, scale)
  x[i] <- median(duration) 
}

model <- gam(inb ~ s(x))
fittedValues <- data.frame(0, fitted(model))
evsi <- mean(do.call(pmax, fittedValues)) - max(colMeans(fittedValues))
evsi

