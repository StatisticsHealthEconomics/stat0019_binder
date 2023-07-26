library(voi)
library(R2OpenBUGS)

### Model Specification
generate_psa_parameters <- function(n){
  
  ## Initialise parameter matrix
  theta <- matrix(NA, ncol = 19, nrow = n)
  #PSA distributions for the parameters
  #All parameters are normal distributions
  theta[, 1] <- rnorm(n, 10000, 10)
  theta[, 11] <- rnorm(n, 15000, 10)
  theta[, 2] <- rnorm(n, 0.1, 0.02)
  theta[, 12] <- rnorm(n, .08, 0.02)
  theta[, 3] <- rnorm(n, 5.2, 1)
  theta[, 13] <- rnorm(n, 6.1, 1)
  theta[, 4] <- rnorm(n, 4000, 2000)
  theta[, 8] <- rnorm(n, 0.25, 0.01)
  theta[, 17] <- rnorm(n, 0.20, 0.05)
  theta[, 9] <- rnorm(n, -0.1, 0.02)
  theta[, 18] <- rnorm(n, -0.1, 0.02)
  theta[, 10] <- rnorm(n, 0.5, 0.2)
  theta[, 19] <- rnorm(n, 0.5, 0.2)
  
  ## Correlated parameters PSA
  Mu.1 <- c(.7, .8, 3, 3)
  Mu.2 <- c(.3, .3)
  
  # Correlation Matrices - standard deviations
  sd.5 <- 0.1
  sd.6 <- 0.1
  sd.7 <- 0.5
  sd.14 <- 0.1
  sd.15 <- 0.05
  sd.16 <- 1
  rho <- 0.6 # Correlation
  
  # Definie variance matrices
  S.1 <- matrix(
    #Vector containing all the values to enter into the matrix
    c(sd.5^2,rho*sd.14*sd.5,rho*sd.5*sd.7,rho*sd.5*sd.16,
      rho*sd.5*sd.14,sd.14^2,rho*sd.14*sd.7,rho*sd.14*sd.16,
      rho*sd.5*sd.7,sd.14*rho*sd.7,sd.7^2,rho*sd.16*sd.7,
      rho*sd.5*sd.16,sd.14*rho*sd.16,rho*sd.7*sd.16,sd.16^2),
    #Specify the size of the matrix
    nrow=4)
  
  
  S.2<-matrix(
    #Vector containing all the values to enter into the matrix
    c(sd.6^2,rho*sd.6*sd.15,
      rho*sd.6*sd.15,sd.15^2),
    #Specify the size of the matrix
    nrow=2)
  
  theta[, c(5, 14, 7, 16)] <- MASS::mvrnorm(n, Mu.1, S.1)
  theta[, c(6, 15)] <- MASS::mvrnorm(n, Mu.2, S.2)
  
  colnames(theta) <- paste("theta_", 1:19, sep = "")
  return(data.frame(theta))
}

calculate_costs_effects <- function(theta_1, theta_2, theta_3,
                                    theta_4, theta_5, theta_6,
                                    theta_7, theta_8, theta_9,
                                    theta_10, theta_11, theta_12,
                                    theta_13, theta_14, theta_15,
                                    theta_16, theta_17, theta_18,
                                    theta_19){
  ## Effects
  e_treatment1 <- (theta_5*theta_6*theta_7)+(theta_8*theta_9*theta_10)
  e_treatment2 <- (theta_14*theta_15*theta_16)+(theta_17*theta_18*theta_19)
  e<-c(e_treatment1, e_treatment2)
  ## Costs
  c_treatment1 <- (theta_1+theta_2*theta_3*theta_4)
  c_treatment2 <- (theta_11+theta_12*theta_13*theta_4)
  c<-c(c_treatment1, c_treatment2)
  
  output <- array(NA, dim = c(2, length(e)),
                  dimnames = list(c("Effects", "Costs"),
                                  c("SoC", "Novel")))
  output[1, ] <- e
  output[2, ] <- c
  return(output)
  
}

calculate_net_benefit <- function(
    costs_effects,
    wtp)
{
  if(!is.null(dim(costs_effects))){
    nb <- wtp * costs_effects[, 1, ] - 
      costs_effects[, 2, ]
  }
  
  return(nb)
}

## Run the model for initial PSA
n_psa_size <- 10000
params <- generate_psa_parameters(n_psa_size)
m_costs_effects <- array(NA, dim = c(n_psa_size, 2, 2))
for(s in 1:n_psa_size){
  formals(calculate_costs_effects) <- params[s, ]
  m_costs_effects[s, , ] <- calculate_costs_effects()
}
# Set the column names for the output
dimnames(m_costs_effects)[2:3] <- dimnames(
  calculate_costs_effects()
)

### For EVSI analysis, we need to compute a model object to save the model
# results
model_output <- list(e = m_costs_effects[, "Effects", ],
                     c = m_costs_effects[, "Costs", ],
                     k = seq(0, 50000, length.out = 501))

### EVPI
EVPI <- evpi(model_output)
K = 10000
EVPI$evpi[EVPI$k == K]

### EVSI
# Create function to generate data
data_sim <- function(inputs, n = 250){
  if(length(n) > 1){
    n <- n[1]
  }
  theta_5 <- inputs[, "theta_5"]
  theta_14 <- inputs[, "theta_14"]
  
  # Specify individual level standard devation
  sd.X.5 <- 0.2
  sd.X.14 <- 0.2
  rho.X <- 0.6
    sd.14 <- sd.X.14 * sqrt(1 - rho.X * rho.X)
  X <- matrix(NA, nrow = length(theta_5), ncol = 2 * n)
  for(i in 1:length(theta_5)){
    X[i, 1:n] <- rnorm(n, mean = theta_5[i], sd = sd.X.5)
    mu.14 <- theta_14[i] + sd.X.14 / sd.X.5 * rho.X * (X[i, 1:n] - theta_5[i])
    X[i, (n + 1):(2 * n)] <- rnorm(n, mean = mu.14, sd = sd.14)
  }
  
  return(data.frame(X)) 
}

data_analysis_fn <- function(data, args, pars){
  Model<- "model{
    #Correlated PSA distributions requires multivariate normal distribution
    cor.theta.1[1:4] ~ dmnorm(Mu.1[], Tau.1[,])
    theta_5 <- cor.theta.1[1]
    theta_7 <- cor.theta.1[3]
    theta_14 <- cor.theta.1[2]
    theta_16 <- cor.theta.1[4]

    ##Sampling distribution for the data
    #Decompose multivariate normal into two scalar data points with conditional mean and variance
    #Decomposition found in Baio (2012), Page 156
    tau.14 <- tau.X.14 / (1 - rho.X * rho.X)
    for(i in 1:N){
      X.5[i] ~ dnorm(theta_5, tau.X.5)
      mu.14[i] <- theta_14 + sqrt(tau.X.5 / tau.X.14) * rho.X *
        (X.5[i] - theta_5)
      X.14[i] ~ dnorm(mu.14[i], tau.14)
    }
  }"
  filein <- file.path(tempdir(),fileext="datmodel.txt")
  cat(Model, file=filein)
  
  data <- as.vector(as.matrix(data))
  data_load <- list(
    Mu.1 = args$Mu.1,
    Tau.1 = args$Tau.1,
    tau.X.14 = args$tau.X.14,
    tau.X.5 = args$tau.X.5,
    rho.X = args$rho.X,
    N = args$n,
    X.5 = data[1:args$n],
    X.14 = data[(args$n + 1):(2 * args$n)]
  )
  
  bugs.data <- bugs(
    data = data_load,
    parameters.to.save = pars,
    model.file = filein,
    inits = NULL,
    n.chains = 1,
    n.iter = args$n.iter,
    n.thin = 1,
    n.burnin = 250)
  
  theta_5 <- bugs.data$sims.matrix[, "theta_5"]
  theta_14 <- bugs.data$sims.matrix[, "theta_14"]
  theta_7 <- bugs.data$sims.matrix[, "theta_7"]
  theta_16 <- bugs.data$sims.matrix[, "theta_16"]
  return(data.frame(theta_5 = theta_5,
                    theta_14 = theta_14,
                    theta_7 = theta_7,
                    theta_16 = theta_16))
}

## Data to be loaded into the analysis function
Mu.1 <- c(.7, .8, 3, 3)
# Correlation Matrices - standard deviations
sd.5 <- 0.1
sd.6 <- 0.1
sd.7 <- 0.5
sd.14 <- 0.1
sd.15 <- 0.05
sd.16 <- 1
rho <- 0.6 # Correlation

# Definie variance matrices
S.1 <- matrix(
  #Vector containing all the values to enter into the matrix
  c(sd.5^2,rho*sd.14*sd.5,rho*sd.5*sd.7,rho*sd.5*sd.16,
    rho*sd.5*sd.14,sd.14^2,rho*sd.14*sd.7,rho*sd.14*sd.16,
    rho*sd.5*sd.7,sd.14*rho*sd.7,sd.7^2,rho*sd.16*sd.7,
    rho*sd.5*sd.16,sd.14*rho*sd.16,rho*sd.7*sd.16,sd.16^2),
  #Specify the size of the matrix
  nrow=4)

# Specify individual level standard devation
sd.X.5 <- 0.2
sd.X.14 <- 0.2
rho.X <- 0.6

args.list <- list(
  Mu.1 = Mu.1,
  Tau.1 = solve(S.1), # precision matrix
  tau.X.14 = 1 / sd.X.14 ^ 2,
  tau.X.5 = 1 / sd.X.5 ^ 2,
  rho.X = rho.X,
  n = 250,
  n.iter = 1000
)

evsi_single_n <- evsi(outputs = model_output,
             inputs = params,
             pars = c("theta_5", "theta_14", "theta_7", "theta_16"),
             pars_datagen = c("theta_5", "theta_14"),
             n = 250,
             method = "mm",
             datagen_fn = data_sim,
             model_fn = calculate_costs_effects,
             analysis_args = args.list,
             analysis_fn = data_analysis_fn,
             par_fn = generate_psa_parameters,
             Q = 50)

evsi_multi_n <- evsi(outputs = model_output,
     inputs = params,
     pars = c("theta_5", "theta_14", "theta_7", "theta_16"),
     pars_datagen = c("theta_5", "theta_14"),
     n = seq(30,200, by = 10),
     method = "mm",
     datagen_fn = data_sim,
     model_fn = calculate_costs_effects,
     analysis_args = args.list,
     analysis_fn = data_analysis_fn,
     par_fn = generate_psa_parameters,
     Q = 50)



