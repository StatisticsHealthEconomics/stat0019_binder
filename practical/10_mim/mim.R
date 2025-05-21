#' Sets the working directory
setwd(here::here("10_mim/"))

#' First, we need to load the necessary R packages.
#' If you haven't installed them, you'll need to do so.
#' install.packages(
#' c("dplyr", "tidyr", "boot", "copula", "rstanarm", "remotes"))
#'
#' remotes::install_github("n8thangreen/simcovariates") # For gen_data
#' remotes::install_github("StatisticsHealthEconomics/outstandR")
#' More simply, you can use r-universe.dev to install the packages as
#' install.packages(
#'   "outstandR", 
#'   repos = c(
#'     "https://statisticshealtheconomics.r-universe.dev", 
#'     "https://cloud.r-project.org")
#' )
#' install.packages(
#'   "simcovariates", 
#'   repos = c(
#'     "https://statisticshealtheconomics.r-universe.dev", 
#'     "https://cloud.r-project.org")
#' )
#' NB: In the docker VM, these are already installed so you can just proceed

#' Loads relevant packages
library(outstandR)
library(simcovariates) # For gen_data
library(dplyr)
library(tidyr)

# For reproducibility of simulated data
set.seed(123)

#' ## Part 1: Data Simulation & Preparation - Binary Outcomes
#' ### 1.1 Simulation Parameters
N <- 200             # Sample size per trial

#' Active treatment vs. placebo allocation ratio (2:1 implies ~2/3 on active)
allocation <- 2/3

#' Conditional log-OR for active treatment vs. common comparator C
b_trt <- log(0.17)

#' Conditional log-OR for each unit increase in prognostic variables (X3, X4)
b_X <- -log(0.5)

#' Conditional log-OR for interaction term (treatment * effect modifier) 
#' for X1, X2
b_EM <- -log(0.67)

#' Mean of prognostic factors (X3, X4) in AC trial
meanX_AC <- c(0.45, 0.45)

#' Mean of prognostic factors (X3, X4) in BC trial (DIFFERENT from AC)
meanX_BC <- c(0.6, 0.6)      
meanX_EM_AC <- c(0.45, 0.45) # Mean of effect modifiers (X1, X2) in AC trial

#' Mean of effect modifiers (X1, X2) in BC trial (DIFFERENT from AC)
meanX_EM_BC <- c(0.6, 0.6)  
sdX <- c(0.4, 0.4)           # Standard deviation of prognostic factors
sdX_EM <- c(0.4, 0.4)        # Standard deviation of effect modifiers
corX <- 0.2                  # Covariate correlation coefficient
b_0 <- -0.6                  # Baseline intercept coefficient on logit scale

#' **Effect Modifiers vs. Prognostic Variables:**\
#' - **Prognostic variables** (`X3`, `X4` here) predict the outcome 
#'   regardless of treatment.
#' - **Effect modifiers** (`X1`, `X2` here) change the magnitude or 
#'   direction of the treatment effect.
#' Differences in the distribution of effect modifiers between trials are a 
#' key reason for population adjustment.

#' ### 1.2 Generate IPD for AC Trial (Binary Outcome)
#' We simulate Individual Patient Data (IPD) for a trial comparing 
#' treatments A and C.
ipd_trial_bin <- gen_data(N, 
                          b_trt, 
                          b_X, 
                          b_EM, 
                          b_0 = b_0,
                          meanX_AC, 
                          sdX, 
                          meanX_EM_AC, 
                          sdX_EM, 
                          corX, 
                          allocation,
                          family = binomial("logit"))

#' Treatment 'trt' is 0 or 1
#' We map 0 to 'C' (comparator) and 1 to 'A' (new treatment)
ipd_trial_bin$trt <- factor(ipd_trial_bin$trt, labels = c("C", "A"))

#' Visualise and summarise the simulated data
head(ipd_trial_bin)
summary(ipd_trial_bin)

#' The `ipd_trial_bin` dataframe contains patient-level data: covariates 
#' (`X1`-`X4`), treatment assignment (`trt`), and outcome (`y`).

#' ### 1.3 Generate ALD for BC Trial (Binary Outcome)
#' For the BC trial (comparing B vs C), we only have Aggregate Level Data (ALD).
#' We first simulate IPD for BC and then summarize it.
#' The key here is that `meanX_BC` and `meanX_EM_BC` are different from the AC trial, creating a population imbalance.

#' Simulate IPD for BC trial (using BC trial's covariate means)
BC_IPD_bin <- gen_data(N, 
                       b_trt, 
                       b_X, 
                       b_EM, 
                       b_0,
                       meanX_BC, # Using BC means
                       sdX, 
                       meanX_EM_BC, # Using BC means
                       sdX_EM, 
                       corX, 
                       allocation,
                       family = binomial("logit"))

BC_IPD_bin$trt <- factor(BC_IPD_bin$trt, labels = c("C", "B")) # 0=C, 1=B

#' Now, aggregate BC_IPD_bin to create ald_trial_bin. This mimics having 
#' only published summary statistics.

#' Covariate summaries (mean, sd for X1-X4, assumed same across 
#' arms in BC trial for simplicity)
cov_summary_bin <- BC_IPD_bin %>%
  select(X1, X2, X3, X4) %>% # Select covariate columns
  summarise(across(everything(), list(mean = mean, sd = sd))) %>%
  pivot_longer(everything(), names_to = "stat_var", values_to = "value") %>%
  # 'stat_var' will be like "X1_mean", "X1_sd". We need to separate these.
  separate(stat_var, into = c("variable", "statistic"), sep = "_") %>%
  # Covariate summaries are often reported for the overall trial population
  mutate(trt = NA_character_) 

#' Outcome summaries (number of events 'sum', mean proportion 'mean', 
#' sample size 'N' for y by trt)
outcome_summary_bin <- BC_IPD_bin %>%
  group_by(trt) %>%
  summarise(
    sum_y = sum(y),      # Number of events
    mean_y = mean(y),    # Proportion of events
    N = n()              # Sample size in this arm
  ) %>%
  ungroup() %>%
  pivot_longer(cols = -trt, names_to = "stat_var", values_to = "value") %>%
  # 'stat_var' will be "sum_y", "mean_y", "N". We need to parse this.
  mutate(
    variable = case_when(
      grepl("_y$", stat_var) ~ "y", # If it ends with _y, variable is y
      stat_var == "N" ~ NA_character_, # For N, variable can be NA
      TRUE ~ stat_var # Default
    ),
    statistic = case_when(
      grepl("sum_", stat_var) ~ "sum",
      grepl("mean_", stat_var) ~ "mean",
      stat_var == "N" ~ "N",
      TRUE ~ stat_var # Default
    )
  ) %>%
  select(variable, statistic, value, trt)

#' Combine covariate and outcome summaries for the final ALD structure
ald_trial_bin <- bind_rows(cov_summary_bin, outcome_summary_bin) %>%
  select(variable, statistic, value, trt)

print(as.data.frame(ald_trial_bin))

#' The `ald_trial_bin` is in a 'long' format with columns: `variable` 
#' (e.g., "X1", "y"), `statistic` (e.g., "mean", "sd", "sum", "N"), `value`, 
#' and `trt` (treatment arm, or NA if overall). This is the format 
#' `{outstandR}` expects.

#' ## Part 2: Model Fitting - Binary Outcomes
#' Now we use `{outstandR}` to perform population adjustments.
#' We'll compare treatment A (from AC trial IPD) with treatment B (from BC 
#' trial ALD), using C as the common anchor. The target population for 
#' comparison will be the BC trial population.

#' ### 2.1 Define the Model Formula
#' The model formula specifies the relationship between the outcome (`y`), 
#' prognostic variables (`X3`, `X4`), treatment (`trt`), and effect 
#' modifiers (`X1`, `X2`).
#' For a binary outcome with a logit link, the model is:
#' This translates to the R formula: 
#' `y ~ X3 + X4 + trt + trt:X1 + trt:X2` 
#' (The intercept $\beta_0$ is implicit).
lin_form_bin <- as.formula("y ~ X3 + X4 + trt + trt:X1 + trt:X2")

#' ### 2.2 Matching-Adjusted Indirect Comparison (MAIC)
#' MAIC reweights the IPD from the AC trial so that the mean covariate values 
#' of the effect modifiers match those of the BC trial population.

#' MAIC involves bootstrapping, which can take a moment. The number of 
#' bootstrap replicates can sometimes be controlled in `strategy_maic()` 
#' for speed, e.g. n_boot = 100 for a quick check, but higher (e.g., 1000) 
#' is better for stable results. We'll use the default for now.
out_maic_bin <- outstandR(
  ipd_trial = ipd_trial_bin, 
  ald_trial = ald_trial_bin,
  strategy = strategy_maic(
    formula = lin_form_bin,
    family = binomial(link = "logit") 
    # If your package allows, you might add:
    # , n_boot = 200 # for faster demo
  )
)

#' The MAIC results (default: Log-Odds Ratio scale):
print(out_maic_bin)

#' The output provides `contrasts` (e.g., A vs B) and `absolute_effects` 
#' in the target (BC) population. By default, for `binomial(link="logit")`, 
#' the effect measure is the log-odds ratio.

#' ### 2.3 Changing the Outcome Scale (MAIC Example)
#' Often, we want results on a different scale, like log-relative risk 
#' or risk difference. The `scale` argument in `outstandR()` allows this.
out_maic_bin_lrr <- outstandR(
  ipd_trial = ipd_trial_bin, 
  ald_trial = ald_trial_bin,
  strategy = strategy_maic(
    formula = lin_form_bin,
    family = binomial(link = "logit")
  ),
  scale = "log_relative_risk" # Key change!
)

#' The MAIC results on the log-relative risk scale,
print(out_maic_bin_lrr)

#' **Your Turn!** Try getting MAIC results on the **risk difference** scale.
#' Hint: `scale = "risk_difference"`.
out_maic_bin_rd <- outstandR(
  ipd_trial = ipd_trial_bin, 
  ald_trial = ald_trial_bin,
  strategy = strategy_maic(
    formula = lin_form_bin,
    family = binomial(link = "logit")
  ),
  scale = "risk_difference" # Key change!
)

#' The MAIC results on the risk difference scale,
print(out_maic_bin_rd)

#' ### 2.4 Parametric G-computation with Maximum Likelihood (G-comp ML)
#' G-computation fits an outcome regression model to the IPD (AC trial) and 
#' then uses this model to predict outcomes for each patient *as if* they had 
#' received treatment A and *as if* they had received treatment C, but 
#' standardized to the covariate distribution of the target (BC) population.
out_gcomp_ml_bin <- outstandR(
  ipd_trial = ipd_trial_bin, 
  ald_trial = ald_trial_bin,
  strategy = strategy_gcomp_ml(
    formula = lin_form_bin,
    family = binomial(link = "logit")
  )
)

print(out_gcomp_ml_bin)

#' ## Part 3: Adapting for Continuous Outcomes
#' What if our outcome is continuous, like change in blood pressure or a 
#' quality-of-life score? The principles are similar, but we need to adjust 
#' the data generation and model specification.

#' ### 3.1 Simulate Continuous Data
#' We'll use `family = gaussian("identity")` for the `gen_data` function. 
#' We might also adjust some coefficients to be more sensible for a 
#' continuous scale.

#' Adjust some parameters for a continuous outcome
b_0_cont <- 5      # Intercept on the continuous scale
b_trt_cont <- -1.5 # Mean difference for treatment A vs C
b_X_cont <- 0.5    # Effect of prognostic vars on continuous outcome

#' Effect of effect modifiers on treatment effect (continuous)
b_EM_cont <- 0.3

#' #### 3.1.1 IPD for AC Trial (Continuous)
ipd_trial_cont <- gen_data(N, 
                           b_trt_cont, 
                           b_X_cont, 
                           b_EM_cont, 
                           b_0_cont,
                           meanX_AC, 
                           sdX, 
                           meanX_EM_AC, 
                           sdX_EM, 
                           corX, 
                           allocation,
                           family = gaussian("identity")) # Key change!

ipd_trial_cont$trt <- factor(ipd_trial_cont$trt, labels = c("C", "A"))

head(ipd_trial_cont)
summary(ipd_trial_cont$y)


#' #### 3.1.2 ALD for BC Trial (Continuous)
BC_IPD_cont <- gen_data(N, 
                        b_trt_cont, 
                        b_X_cont, 
                        b_EM_cont, 
                        b_0_cont,
                        meanX_BC, # Using BC means
                        sdX, 
                        meanX_EM_BC, # Using BC means
                        sdX_EM, 
                        corX, 
                        allocation,
                        family = gaussian("identity")) # Key change!

BC_IPD_cont$trt <- factor(BC_IPD_cont$trt, labels = c("C", "B"))

#' Aggregate BC_IPD_cont for ALD
#' Covariate summaries structure remains the same
cov_summary_cont <- BC_IPD_cont %>%
  select(X1, X2, X3, X4) %>%
  summarise(across(everything(), list(mean = mean, sd = sd))) %>%
  pivot_longer(everything(), names_to = "stat_var", values_to = "value") %>%
  separate(stat_var, into = c("variable", "statistic"), sep = "_") %>%
  mutate(trt = NA_character_)

#' Outcome summaries for continuous data: mean, sd, N for y by trt
outcome_summary_cont <- BC_IPD_cont %>%
  group_by(trt) %>%
  summarise(
    mean_y = mean(y),    # Mean outcome
    sd_y = sd(y),        # Standard deviation of outcome
    N = n()              # Sample size
  ) %>%
  ungroup() %>%
  pivot_longer(cols = -trt, names_to = "stat_var", values_to = "value") %>%
  mutate(
    variable = case_when(
      grepl("_y$", stat_var) ~ "y",
      stat_var == "N" ~ NA_character_,
      TRUE ~ stat_var
    ),
    statistic = case_when(
      grepl("mean_", stat_var) ~ "mean",
      grepl("sd_", stat_var) ~ "sd", # Changed from sum to sd
      stat_var == "N" ~ "N",
      TRUE ~ stat_var
    )
  ) %>%
  select(variable, statistic, value, trt)

ald_trial_cont <- bind_rows(cov_summary_cont, outcome_summary_cont) %>%
  select(variable, statistic, value, trt)

print(as.data.frame(ald_trial_cont))

#' ### 3.2 Model Fitting for Continuous Outcomes
#' The model formula structure can remain the same if we assume linear 
#' relationships. The key change is in the `family` argument of the 
#' strategy function.
lin_form_cont <- as.formula("y ~ X3 + X4 + trt + trt:X1 + trt:X2")

#' Let's use G-computation ML as an example.
out_gcomp_ml_cont <- outstandR(
  ipd_trial = ipd_trial_cont, 
  ald_trial = ald_trial_cont,
  strategy = strategy_gcomp_ml(
    formula = lin_form_cont,
    family = gaussian(link = "identity") # Key change!
  )
  # For Gaussian family, the default scale is typically
  # "mean_difference", # which is often what we want.
  # We could explicitly state: scale = "mean_difference"
)

print(out_gcomp_ml_cont)

#' **Your Turn!** Try applying **MAIC** to the continuous outcome data.
#' 1. Use `family = gaussian(link = "identity")` within `strategy_maic()`.
#' 2. What `scale` would be appropriate if not the default? 
#'    (e.g., `"mean_difference"`)

# # Solution for MAIC with continuous data:
# out_maic_cont <- outstandR(
#   ipd_trial = ipd_trial_cont,
#   ald_trial = ald_trial_cont,
#   strategy = strategy_maic(
#     formula = lin_form_cont,
#     family = gaussian(link = "identity")
#   ),
#   scale = "mean_difference"
# )
# print(out_maic_cont)

#' ### 2.5 Other Methods
#' `{outstandR}` supports other methods. Here's how you might call them.
#' These are set to `eval=FALSE` to save time in this practical.
#' - **Simulated Treatment Comparison (STC):** 
#'   A conventional outcome regression approach.
# out_stc_bin <- outstandR(
#   ipd_trial = ipd_trial_bin,
#   ald_trial = ald_trial_bin,
#   strategy = strategy_stc(
#     formula = lin_form_bin,
#     family = binomial(link = "logit")
#   )
# )
# print(out_stc_bin)

#' - **Bayesian G-computation (G-comp Bayes):** 
#'   Similar to G-comp ML but uses Bayesian methods (e.g., MCMC via 
#'   `rstanarm`), which can better propagate uncertainty but is 
#'   computationally more intensive.
# # This would require rstanarm and can be slow.
# out_gcomp_stan_bin <- outstandR(
#   ipd_trial = ipd_trial_bin,
#   ald_trial = ald_trial_bin,
#   strategy = strategy_gcomp_stan(
#     formula = lin_form_bin,
#     family = binomial(link = "logit")
#     # For a faster demo if options are passed through:
#     # stan_args = list(iter = 500, chains = 2, refresh = 0)
#   )
# )
# print(out_gcomp_stan_bin)

#' - **Multiple Imputation Marginalisation (MIM):** 
#'   Another approach for marginalization.
# out_mim_bin <- outstandR(
#   ipd_trial = ipd_trial_bin,
#   ald_trial = ald_trial_bin,
#   strategy = strategy_mim(
#     formula = lin_form_bin,
#     family = binomial(link = "logit")
#   )
# )
# print(out_mim_bin)

#' ## Part 4: Understanding Output & Wrap-up
#' Let's briefly revisit one of the binary outcome results to understand the 
#' structure of the `{outstandR}` output.
str(out_maic_bin)

#' The output object (here `out_maic_bin`) is a list containing:
#' - `$contrasts`: This list provides the estimated treatment effects 
#'    (e.g., mean difference, log-OR), their variances, and confidence 
#'    intervals for each pairwise comparison, adjusted to the target population 
#'    (BC trial).
#' - `$contrasts$means$AB`: The estimated effect of A versus B. This is 
#'    often the primary interest.
#' - `$contrasts$means$AC`: The estimated effect of A versus C.
#' - `$contrasts$means$BC`: The estimated effect of B versus C (usually 
#'    derived directly from the ALD).
#' - `$absolute_effects`: This list provides the estimated mean outcome for 
#'    each treatment (A, B, C) in the target population. This can be useful for 
#'    understanding the baseline and treated outcomes.
#' 
#' For example, to extract the estimated log-odds ratio for A vs. B and 
#' its variance:
log_or_AB <- out_maic_bin$contrasts$means$AB
variance_log_or_AB <- out_maic_bin$contrasts$variances$AB

cat(paste("Estimated Log-OR for A vs. B:", round(log_or_AB, 3), "\n"))
cat(paste("Variance of Log-OR for A vs. B:", round(variance_log_or_AB, 3), "\n"))

#' The vignette for `{outstandR}` (which this practical is based on) shows how 
#' to combine results from multiple methods into tables and forest plots for a 
#' comprehensive comparison. This is highly recommended for actual analyses.

#' ### Key Takeaways
#' - Population adjustment is crucial when comparing treatments indirectly using 
#'   IPD and ALD from trials with different patient characteristics (especially 
#'   different distributions of effect modifiers).
#' - The `{outstandR}` package provides a unified interface (`outstandR()` 
#'   function) to apply various adjustment methods.
#' - You need to:
#'     1. Prepare your IPD (for the "anchor" trial, e.g., AC) and ALD 
#'        (for the "comparator" trial, e.g., BC, which also serves as the 
#'        target population).
#'     2. Define an appropriate model `formula`.
#'     3. Choose a `strategy_*()` function corresponding to the desired 
#'        adjustment method (MAIC, STC, G-comp, etc.).
#'     4. Specify the outcome `family` (e.g., `binomial()`, `gaussian()`) 
#'        within the strategy.
#'     5. Optionally, use the `scale` argument in `outstandR()` to transform 
#'        results to a desired effect measure scale.
#' - The methods can be adapted for different outcome types (binary, 
#'   continuous, count, time-to-event, though we only covered binary 
#'   and continuous here).
