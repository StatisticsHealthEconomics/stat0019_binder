# 0. Sets up environment
setwd(ENTER PATH HERE)
library(BCEA)
library(R2OpenBUGS)

# 1. Loads the results of the Bayesian model performed for the 6 statins. 
# NB this contains variables effect and cost.tot that can be used to perform the  economic analysis using BCEA
load("statins_base.Rdata")

# Now loads the results of a new Bayesian model performed for the 6 statins
# In this case, one of the priors has been modified to consider a more robust alternative
# (in particular that's a Half-Cauchy --- instead of a Normal --- model for the 
# statin-specific effectiveness)
load("statins_HC.Rdata")

# 2. Now uses BCEA to perform the two economic analyses
interventions <- c("Atorvastatin","Fluvastatin","Lovastatin","Pravastatin",
                   "Rosuvastatin","Simvastatin")
m1 <- bcea(statins_base$sims.list$effect,statins_base$sims.list$cost.tot,ref=1,interventions=interventions)
m2 <- bcea(statins_HC$sims.list$effect,statins_HC$sims.list$cost.tot,ref=1,interventions=interventions)

# 3. Uses BCEA to do structural PSA
# Defines suitable lists containing:
models <- list(statins_base,statins_HC)             # the BUGS models
effects <- list(statins_base$sims.list$effect,      # the simulated effectiveness variables
                statins_HC$sims.list$effect)
costs <- list(statins_base$sims.list$cost.tot,      # the simulated cost variables
              statins_HC$sims.list$cost.tot)

# Finally uses BCEA to perform the structural PSA to consider the base and HC models
m3 <- struct.psa(models,effects,costs,ref=1,interventions=interventions)

