###################################################################
## Install the relevant packages
##################################################################
install.packages("SimJoint")
install.packages("boot")

###################################################################
## Exercise 1
###################################################################

#### Question 1 ####
S <- 1000 
M <- 100
theta_3 <- runif(S, 0.2, 0.3) # Hypothetical distribution for theta_3


#### Question 2 ####
S <- 1000 
# Correlated joint distribution for theta_4 and theta_5 
# (Column 1: theta_4, Column 2: theta_5)
theta_4_5 <- MASS::mvrnorm(S,
                           c(5,6),
                           matrix(c(0.3, 0.1, 0.1, 0.5), nrow = 2)) 

#### Question 3 ####

###################################################################
## Exercise 2
###################################################################

#### Question 1 ####
S <- 1000
theta_6 <- rbeta(S, 70, 15) # Hypothetical distribution for theta_3

calculate_beta_parameters <- function(mean, sd){
  # Function to estimate beta parameters from mean and standard deviation
  shape1 <- ((1 - mean) / sd ^ 2 - 1 / mean) * mean ^ 2
  shape2 <- shape1 * (1 / mean - 1)
  
  # Return the calculated parameters.
  return(list(shape1 = shape1,
              shape2 = shape2))
}


#### Question 2 ####



