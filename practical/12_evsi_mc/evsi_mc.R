library(voi)
library(R2OpenBUGS)
library(ggplot2)

## Set working directory to source the model code
setwd(here::here("12_evsi_mc"))
# Change if working directory 
source("model_functions_example.R")

## Run the health economic model for initial PSA
n_psa_size <- 10000
# Generate parameter simulations
params <- generate_psa_parameters(n_psa_size)
# Compute costs and effects for each parameter simulation
m_costs_effects <- array(NA, dim = c(n_psa_size, 2, 2))
for(s in 1:n_psa_size){
  # Set the values of the calculate_costs_effects function
  formals(calculate_costs_effects) <- params[s, ] 
  # Compute costs and effects
  m_costs_effects[s, , ] <- calculate_costs_effects()
}
# Set the column names for the output
dimnames(m_costs_effects)[2:3] <- dimnames(
  calculate_costs_effects()
)

### For EVSI analysis within the voi package
# We need to specific a model object with three elements:
# e: the matrix of effects for each intervention
# c: the matrix of costs for each intervention
# k: a vector of willingness to pay values
model_output <- list(e = m_costs_effects[, "Effects", ],
                     c = m_costs_effects[, "Costs", ],
                     k = seq(0, 50000, length.out = 501))

### EVPI within voi package
evpi_df <- evpi(model_output)
K = 10000 # willingness to pay threshold
evpi_df$evpi[evpi_df$k == K]

### calculate EVPPI within voi package
evppi_df <- evppi(model_output,
                  params,
                  pars = list("params_of_interest" = 
                                c("theta_5", "theta_14", "theta_7", "theta_16")),
                  method = "earth")

### Computing EVSI
# Create function to generate data
data_sim <- function(inputs, n = 250){
  if(length(n) > 1){ # If multiple sample sizes are passed, select the first
    n <- n[1]
  }
  # Extract the key parameters
  theta_5 <- inputs[, "theta_5"]
  theta_14 <- inputs[, "theta_14"]
  
  # Specify individual level standard deviation and correlation
  sd_X_5 <- 0.2
  sd_X_14 <- 0.2
  rho_X <- 0.6
  
  # Generate data as a conditional and marginal
  sd_14 <- sd_X_14 * sqrt(1 - rho_X * rho_X) # Conditional standard deviation
  # Save data as a row for each parameter simulation
  X <- matrix(NA, nrow = length(theta_5), ncol = 2 * n)
  for(i in 1:length(theta_5)){
    # First n values simulated from marginal conditional on theta_5
    X[i, 1:n] <- rnorm(n, mean = theta_5[i], sd = sd_X_5)
    # Define conditional mean for each individual
    mu_14 <- theta_14[i] + sd_X_14 / sd_X_5 * rho_X * (X[i, 1:n] - theta_5[i])
    # Second set of n values simulated from conditional distribution
    X[i, (n + 1):(2 * n)] <- rnorm(n, mean = mu_14, sd = sd_14)
  }
  # Return X as a data.frame
  return(data.frame(X)) 
}

data_analysis_fn <- function(data, args, pars){
  # Function to provide updated parameter values conditional on data
  Model<- "model{
    # Correlated PSA distributions requires multivariate normal distribution
    cor_theta_1[1:4] ~ dmnorm(Mu_1[], Tau_1[,])
    theta_5 <- cor_theta_1[1]
    theta_7 <- cor_theta_1[3]
    theta_14 <- cor_theta_1[2]
    theta_16 <- cor_theta_1[4]

    ## Sampling distribution for the data
    # As decomposed multivariate normal - remember this is in terms of PRECISON
    tau_14 <- tau_X_14 / (1 - rho_X * rho_X)
    for(i in 1:N){
      X_5[i] ~ dnorm(theta_5, tau_X_5)
      mu_14[i] <- theta_14 + sqrt(tau_X_5 / tau_X_14) * rho_X *
        (X_5[i] - theta_5)
      X_14[i] ~ dnorm(mu_14[i], tau_14)
    }
  }"
  filein <- file.path(tempdir(),fileext="datmodel.txt")
  cat(Model, file=filein)
  
  # Data in BUGS needs to be as a numeric vector
  data <- as.vector(as.matrix(data))
  data_load <- list(
    Mu_1 = args$Mu_1, # Prior mean for parameters
    Tau_1 = args$Tau_1, # Prior precision matrix for parameters
    tau_X_14 = args$tau_X_14, # individual-level standard deviations
    tau_X_5 = args$tau_X_5,
    rho_X = args$rho_X,
    N = args$n, # Sample size
    X_5 = data[1:args$n], # Data
    X_14 = data[(args$n + 1):(2 * args$n)]
  )
  
  # Update parameters using BUGS
  bugs.data <- bugs(
    data = data_load,
    parameters.to.save = pars,
    model.file = filein,
    inits = NULL,
    n.chains = 1,
    n.iter = args$n.iter, # number of iteractions can be changed to improve accuracy
    n.thin = 1,
    n.burnin = 250)
  
  # Save the parameters of interest
  theta_5 <- bugs.data$sims.matrix[, "theta_5"]
  theta_14 <- bugs.data$sims.matrix[, "theta_14"]
  theta_7 <- bugs.data$sims.matrix[, "theta_7"]
  theta_16 <- bugs.data$sims.matrix[, "theta_16"]
  # Return parameters as data frame
  return(data.frame(theta_5 = theta_5,
                    theta_14 = theta_14,
                    theta_7 = theta_7,
                    theta_16 = theta_16))
}

## Specify the data that are required in the analysis function
Mu_1 <- c(.7, .8, 3, 3) # Prior mean
# Define prior precision matrix
sd_5 <- 0.1
sd_14 <- 0.1
sd_7 <- 0.5
sd_16 <- 1
rho <- 0.6 # Correlation

# Definie variance matrices
S_1 <- matrix(
  #Vector containing all the values to enter into the matrix
  c(sd_5^2,rho*sd_14*sd_5,rho*sd_5*sd_7,rho*sd_5*sd_16,
    rho*sd_5*sd_14,sd_14^2,rho*sd_14*sd_7,rho*sd_14*sd_16,
    rho*sd_5*sd_7,sd_14*rho*sd_7,sd_7^2,rho*sd_16*sd_7,
    rho*sd_5*sd_16,sd_14*rho*sd_16,rho*sd_7*sd_16,sd_16^2),
  #Specify the size of the matrix
  nrow=4)

# Specify individual level standard devation
sd_X_5 <- 0.2
sd_X_14 <- 0.2
rho_X <- 0.6

args_list <- list(
  Mu_1 = Mu_1,
  Tau_1 = solve(S_1), # precision matrix
  tau_X_14 = 1 / sd_X_14 ^ 2,
  tau_X_5 = 1 / sd_X_5 ^ 2,
  rho_X = rho_X,
  n = 250,
  n.iter = 1000
)

## Calculate EVSI using the evsi function
evsi_single_n <- evsi(outputs = model_output, # load the costs, effects and WTP
             inputs = params, # load the parameters
             # define the parameters that will be updated
             pars = c("theta_5", "theta_14", "theta_7", "theta_16"),
             # define the parameters that will be used to generate the data
             pars_datagen = c("theta_5", "theta_14"),
             n = 250, # sample size to compute the data
             method = "mm", # use the moment matching method
             datagen_fn = data_sim, # load the data generation function
             # load the function to calculate costs and effects
             model_fn = calculate_costs_effects, 
             # load the list of elements needed for the analysis function
             analysis_args = args_list, 
             analysis_fn = data_analysis_fn, # load data analysis function
             par_fn = generate_psa_parameters, # load parameter simulation function
             Q = 50 # number of nested samples
             )

# Calculate multiparameter EVSI - similar syntax
evsi_multi_n <- evsi(outputs = model_output,
     inputs = params,
     pars = c("theta_5", "theta_14", "theta_7", "theta_16"),
     pars_datagen = c("theta_5", "theta_14"),
     n = seq(50,500, by = 10), # specify multiple sample sizes
     method = "mm",
     datagen_fn = data_sim,
     model_fn = calculate_costs_effects,
     analysis_args = args_list,
     analysis_fn = data_analysis_fn,
     par_fn = generate_psa_parameters,
     Q = 50)


### Plotting
source("01_plotting_functions.R")

## Calculates EVPI and EVPPI for plotting
plotting_obj <- evsi.plot.adapt(model_output,
                                params,
                                list("params_of_interest" = 
                                       c("theta_5", "theta_14", "theta_7", "theta_16")),
                                evsi_multi_n,
                                method = "earth")
## Plots EVSI across willingness to pay for each sample size
evsi.wtp.plot(plotting_obj)
## Plots EVSI across sample size for a fixed willingness to pay
evsi.ss.plot(plotting_obj, k = K)
## Plots the probability of a cost effective trial
evsi.prob.plot(plotting_obj, setup = 1e5, pp = 1e3,
               Time = c(0, 10), k = K)
## Plot the expected net benefit of sampling
evsi.enbs.plot(plotting_obj, setup = c(7000,1e5), pp = c(500, 1e3), Time = 10,
               Pop = 2500, k = K)
coss(plotting_obj, setup = 1e5, pp = 1e3, Time = 10, Pop = 2500)
