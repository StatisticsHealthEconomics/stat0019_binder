# Solution to Exercise 1
# Remove weight from regression from nma()
# Remove weight from add_integration()
# Replace prospect_dat with tlv_dat taken from the exercise table
# Run relative_effects on tlv_dat
# IXE_Q2W remains the treatment with greatest probit difference to placebo.

setwd(here::here("08_mlnmr"))


library(multinma)
library(dplyr)
#library(tidyr)

#  Check how many cores current PC has
options(mc.cores = parallel::detectCores())


# Form the IPD dataset
pso_ipd <- transform(plaque_psoriasis_ipd,
                     # Variable transformations
                     bsa = bsa / 100,
                     weight = weight / 10,
                     durnpso = durnpso / 10,
                     prevsys = as.numeric(prevsys),
                     psa = as.numeric(psa),
                     # Treatment classes
                     trtclass = 
                       ifelse(trtc == "PBO", "Placebo",
                              ifelse(trtc %in% c("IXE_Q2W", "IXE_Q4W", "SEC_150", "SEC_300"),  "IL-17 blocker",
                                     ifelse(trtc == "ETN", "TNFa blocker",
                                            ifelse(trtc == "UST", "IL-12/23 blocker", NA)))),
                     # Check complete cases for covariates of interest
                     is_complete = complete.cases(durnpso, prevsys, bsa, weight, psa)
) 

# Only a very small proportion of incomplete rows; we simply remove these
pso_ipd <- subset(pso_ipd, is_complete)

# AgD studies
pso_agd <- transform(plaque_psoriasis_agd,
                     # Variable transformations
                     bsa_mean = bsa_mean / 100, bsa_sd = bsa_sd / 100,
                     weight_mean = weight_mean / 10, weight_sd = weight_sd / 10,
                     durnpso_mean = durnpso_mean / 10, durnpso_sd = durnpso_sd / 10,
                     prevsys = prevsys / 100,
                     psa = psa / 100,
                     # Treatment classes
                     trtclass = 
                       ifelse(trtc == "PBO", "Placebo",
                              ifelse(trtc %in% c("IXE_Q2W", "IXE_Q4W", "SEC_150", "SEC_300"),  "IL-17 blocker",
                                     ifelse(trtc == "ETN", "TNFa blocker",
                                            ifelse(trtc == "UST", "IL-12/23 blocker", NA))))
)

# Create a combined evidence network
pso_net <- combine_network(
  set_ipd(pso_ipd,
          study = studyc,
          trt = trtc,
          r = pasi75,
          trt_class = trtclass),
  set_agd_arm(pso_agd,
              study = studyc,
              trt = trtc,
              r = pasi75_r,
              n = pasi75_n,
              trt_class = trtclass)
)

# Add integration points
pso_net <- add_integration(pso_net,
                           durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
                           prevsys = distr(qbern, prob = prevsys),
                           bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
                           psa = distr(qbern, prob = psa))

# Run the multinma 
pso_fit_FE <- nma(pso_net, 
                  trt_effects = "fixed",
                  link = "probit", 
                  likelihood = "bernoulli2",
                  # Specify the regression model applied to IPD and AgD
                  regression = ~(durnpso + prevsys + bsa + psa)*.trt,
                  class_interactions = "common",
                  prior_intercept = normal(scale = 10),
                  prior_trt = normal(scale = 10),
                  prior_reg = normal(scale = 10),
                  init_r = 0.1,
                  # Do we apply the QR decomposition to the design matrix?
                  QR = TRUE)



# Basic parameter summaries
print(pso_fit_FE, pars = "d")

# Relative effects in the PROSPECT study population
# Specify the characteristics of this population
tlv_dat <- data.frame(
  studyc = "PROSPECT",
  durnpso = 18.6 / 10, durnpso_sd = 8.5 / 10,
  prevsys = 0.802,
  bsa = 17.5 / 100, bsa_sd = 12.4 / 100,
  psa = 0.181)

# Then estimate relative treatment effects in that population
relative_effects(pso_fit_FE, newdata = tlv_dat, 
                 study = studyc, all_contrasts = FALSE)
