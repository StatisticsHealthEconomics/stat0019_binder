library(voi)

## Set working directory to source the model code
setwd(here::here("13_evsi_regression"))
source(here::here("12_evsi_mc","model_functions_example.R"))

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
EVPI <- evpi(model_output)
K = 10000 # willingness to pay threshold
EVPI$evpi[EVPI$k == K]

### Computing EVSI
# Create function to generate data
data_sim_summary <- function(inputs, n = 250){
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
  X <- matrix(NA, nrow = length(theta_5), ncol = 2)
  for(i in 1:length(theta_5)){
    # Summarise data as mean
    X[i, 1] <- mean(rnorm(n, mean = theta_5[i], sd = sd_X_5))
    # Define conditional mean for each individual
    mu_14 <- theta_14[i] + sd_X_14 / sd_X_5 * rho_X * (X[i, 1] - theta_5[i])
    # Summarise data as mean
    X[i, 2] <- mean(rnorm(n, mean = mu_14, sd = sd_14))
  }
  # Return X as a data.frame
  return(data.frame(X)) 
}

## Use evsi function
## Calculate EVSI using the evsi function
evsi_single_n <- evsi(outputs = model_output, # load the costs, effects and WTP
                      inputs = params, # load the parameters
                      # define the parameters that will be updated
                      pars = c("theta_5", "theta_14", "theta_7", "theta_16"),
                      n = 250, # sample size to compute the data
                      method = "gam", # use the regression based method
                      datagen_fn = data_sim_summary, # load the data generation function
)

# Calculate multiparameter EVSI - similar syntax
evsi_multi_n <- evsi(outputs = model_output,
                     inputs = params,
                     pars = c("theta_5", "theta_14", "theta_7", "theta_16"),
                     n = seq(50,500, by = 10), # specify multiple sample sizes
                     method = "gam",
                     datagen_fn = data_sim_summary)

### Plotting
source(here::here("14_evsi_mc","01_plotting_functions.R"))

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



