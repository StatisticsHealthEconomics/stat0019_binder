# Solution to exercise 3
# Change both calls to transform() (one on plaque_psoriasis_ipd
# and the other on plaque_psoriasis_agd) so that "UST" is included
# in the same class as c("IXE_Q2W", "IXE_Q4W", "SEC_150", "SEC_300")
# and relabel this class to  "IL-17 or IL-12/23 blocker"
# Run relative_effects() again.
# Effects of UST vs. IXE_Q2W don't change (both about -0.91) as
# they are in the same class and the effect is 'common'
# Effects of UST vs. ETN change (0.45 in CLEAR and 0.47 in ERASURE) as they are still in different classes.


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
                              ifelse(trtc %in% c("IXE_Q2W", "IXE_Q4W", "SEC_150", "SEC_300", "UST"),  "IL-17 or IL-12/23 blocker",
                                     ifelse(trtc == "ETN", "TNFa blocker", NA))),
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
                              ifelse(trtc %in% c("IXE_Q2W", "IXE_Q4W", "SEC_150", "SEC_300", "UST"),  "IL-17 or IL-12/23 blocker",
                                     ifelse(trtc == "ETN", "TNFa blocker", NA)))
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


# Population-average conditional effects
relative_effects(pso_fit_FE, all_contrasts = TRUE)                 
