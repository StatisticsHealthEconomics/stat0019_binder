## The voi package requires model functions

generate_psa_parameters <- function(n){
  # Function to generate parameter distributions for a probabilistic analysis
  set.seed(123)
  ## Initialise parameter matrix
  theta <- matrix(NA, ncol = 19, nrow = n)
  # PSA distributions for the parameters from normal distributions
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
  
  ## theta 5, 14, 7 and 16 are correlated 
  Mu.1 <- c(.7, .8, 3, 3) # Mean vector for these parameters
  ## theta 6 and 15 are also correlated
  Mu.2 <- c(.3, .3) # Mean vector for these parameters
  
  # Define marginal standard deviations to develop varaince matrix
  sd.5 <- 0.1
  sd.14 <- 0.1
  sd.7 <- 0.5
  sd.16 <- 1
  
  sd.6 <- 0.1
  sd.15 <- 0.05

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
  
  # Generate from multivariate normal distributions
  theta[, c(5, 14, 7, 16)] <- MASS::mvrnorm(n, Mu.1, S.1)
  theta[, c(6, 15)] <- MASS::mvrnorm(n, Mu.2, S.2)
  
  colnames(theta) <- paste("theta_", 1:19, sep = "")
  # Parameters should be a dataframe.
  return(data.frame(theta))
}

calculate_costs_effects <- function(theta_1, theta_2, theta_3,
                                    theta_4, theta_5, theta_6,
                                    theta_7, theta_8, theta_9,
                                    theta_10, theta_11, theta_12,
                                    theta_13, theta_14, theta_15,
                                    theta_16, theta_17, theta_18,
                                    theta_19){
  ## Function to compute costs and effects for a single set of parameters
  ## Effects - from decision tree model
  e_treatment1 <- (theta_5*theta_6*theta_7)+(theta_8*theta_9*theta_10)
  e_treatment2 <- (theta_14*theta_15*theta_16)+(theta_17*theta_18*theta_19)
  e<-c(e_treatment1, e_treatment2)
  ## Costs = from decision tree model
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
  # Calculates net monetary benefit based on costs and effects and
  # a willingness to pay parameter
  if(!is.null(dim(costs_effects))){
    nb <- wtp * costs_effects[, 1, ] - 
      costs_effects[, 2, ]
  }
  
  return(nb)
}
