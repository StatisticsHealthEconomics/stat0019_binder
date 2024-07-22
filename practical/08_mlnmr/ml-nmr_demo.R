setwd(here::here("08_mlnmr"))

library(multinma)
library(dplyr)
#library(tidyr)

#  Check how many cores current PC has
options(mc.cores = parallel::detectCores())


# Inspect the datasets
head(plaque_psoriasis_ipd)

head(plaque_psoriasis_agd)

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

# Plot the evidence network
plot(pso_net, weight_nodes = TRUE, weight_edges = TRUE, show_trt_class = TRUE)

# Add integration points
pso_net <- add_integration(pso_net,
                           durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
                           prevsys = distr(qbern, prob = prevsys),
                           bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
                           weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
                           psa = distr(qbern, prob = psa))

# Run the multinma 
pso_fit_FE <- nma(pso_net, 
                  trt_effects = "fixed",
                  link = "probit", 
                  likelihood = "bernoulli2",
                  # Specify the regression model applied to IPD and AgD
                  regression = ~(durnpso + prevsys + bsa + weight + psa)*.trt,
                  class_interactions = "common",
                  prior_intercept = normal(scale = 10),
                  prior_trt = normal(scale = 10),
                  prior_reg = normal(scale = 10),
                  init_r = 0.1,
                  # Do we apply the QR decomposition to the design matrix?
                  QR = TRUE)



# Basic parameter summaries
print(pso_fit_FE, pars = "d")

# Population-average conditional effects
relative_effects(pso_fit_FE, all_contrasts = TRUE)                 

# Population-average marginal 
marginal_effects(pso_fit_FE, mtype = "link")  

# Average event probabilities
predict(pso_fit_FE, type = "response")        

# Relative effects in the PROSPECT study population
# Specify the characteristics of this population
prospect_dat <- data.frame(
  studyc = "PROSPECT",
  durnpso = 19.6 / 10, durnpso_sd = 13.5 / 10,
  prevsys = 0.9095,
  bsa = 18.7 / 100, bsa_sd = 18.4 / 100,
  weight = 87.5 / 10, weight_sd = 20.3 / 10,
  psa = 0.202)

# Then estimate relative treatment effects in that population
relative_effects(pso_fit_FE, newdata = prospect_dat, 
                 study = studyc, all_contrasts = FALSE)

# Can also generate marginal treatment effects in PROSPECT
# But need to specify the baseline distribution
marginal_effects(pso_fit_FE,
                 mtype = "link",
                 newdata = prospect_dat,
                 study = studyc,
                 baseline = distr(qbeta, 1156, 1509-1156),
                 baseline_type = "response",
                 baseline_level = "aggregate",
                 baseline_trt = "SEC_300",
                 all_contrasts = TRUE)
