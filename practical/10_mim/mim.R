#' Sets the working directory
setwd(here::here("10_mim/"))

library(outstandR)
library(dplyr)
library(ggplot2)

# Loads  & shows the data
data(AC_IPD_binY_contX) # IPD for treatments A vs C
data(BC_ALD_binY_contX) # ALD for treatments B vs C

head(AC_IPD_binY_contX)
print(BC_ALD_binY_contX)

# 1. Outcome Regression Model (used for G-computation)
lin_form_bin <- as.formula(
  "y ~ PF_cont_1 + PF_cont_2 + trt + trt:EM_cont_1 + trt:EM_cont_2")

# 2. Balance Model (used for MAIC weighting)
bal_form_bin <- as.formula("~ EM_cont_1 + EM_cont_2")

formula_list_bin <- list(outcome_model = lin_form_bin, 
                         balance_model = bal_form_bin)

# MAIC analysis
out_maic <-
  outstandR(
    ipd_trial = AC_IPD_binY_contX,
    ald_trial = BC_ALD_binY_contX,
    strategy = strategy_maic(
      formula = list(outcome_model = formula("y ~ trt"),
                     balance_model =  bal_form_bin),
      family = binomial(link = "logit")),
    seed = 12345)

print(out_maic)

# g-computation + MIM + STC analysis
out_gcomp_ml <- outstandR(
  ipd_trial = AC_IPD_binY_contX, 
  ald_trial = BC_ALD_binY_contX,
  strategy = strategy_gcomp_ml(
    formula = formula_list_bin,
    family = binomial(link = "logit") 
  ),
  seed = 12345
)

out_gcomp_bayes <- outstandR(
  ipd_trial = AC_IPD_binY_contX, 
  ald_trial = BC_ALD_binY_contX,
  strategy = strategy_gcomp_bayes(
    formula = formula_list_bin,
    family = binomial(link = "logit") 
  ),
  seed = 12345
)

out_mim <- outstandR(
  ipd_trial = AC_IPD_binY_contX, 
  ald_trial = BC_ALD_binY_contX,
  strategy = strategy_mim(
    formula = formula_list_bin,
    family = binomial(link = "logit") 
  ),
  seed = 12345
)

out_stc <- outstandR(
  ipd_trial = AC_IPD_binY_contX, 
  ald_trial = BC_ALD_binY_contX,
  strategy = strategy_stc(
    formula = formula_list_bin,
    family = binomial(link = "logit") 
  ),
  seed = 12345
)


# Plots all the results
plot(out_stc, out_gcomp_bayes, out_gcomp_ml, out_maic, out_mim)

# MAIC for risk difference
out_maic_rd <- outstandR(
  ipd_trial = AC_IPD_binY_contX, 
  ald_trial = BC_ALD_binY_contX,
  strategy = strategy_maic(
    formula = list(outcome_model = formula("y ~ trt"),
                   balance_model =  bal_form_bin),
    family = binomial(link = "logit")),
  scale = "risk_difference",
  var_method = "sandwich",
  seed = 12345
)

print(out_maic_rd)

# Other models for risk difference
out_gcomp_ml_rd <- outstandR(
  ipd_trial = AC_IPD_binY_contX, 
  ald_trial = BC_ALD_binY_contX,
  strategy = strategy_gcomp_ml(
    formula = formula_list_bin,
    family = binomial(link = "logit") 
  ),
  scale = "risk_difference",
  seed = 12345
)

print(out_gcomp_ml_rd)

out_gcomp_bayes_rd <- outstandR(
  ipd_trial = AC_IPD_binY_contX, 
  ald_trial = BC_ALD_binY_contX,
  strategy = strategy_gcomp_bayes(
    formula = formula_list_bin,
    family = binomial(link = "logit") 
  ),
  scale = "risk_difference",
  seed = 12345
)

print(out_gcomp_bayes_rd)

out_mim_rd <- outstandR(
  ipd_trial = AC_IPD_binY_contX, 
  ald_trial = BC_ALD_binY_contX,
  strategy = strategy_mim(
    formula = formula_list_bin,
    family = binomial(link = "logit") 
  ),
  scale = "risk_difference",
  seed = 12345
)

print(out_mim_rd)

out_stc_rd <- outstandR(
  ipd_trial = AC_IPD_binY_contX, 
  ald_trial = BC_ALD_binY_contX,
  strategy = strategy_stc(
    formula = formula_list_bin,
    family = binomial(link = "logit") 
  ),
  scale = "risk_difference",
  seed = 12345
)

# Compares all the results
plot(out_stc_rd, out_gcomp_bayes_rd, out_gcomp_ml_rd, out_maic_rd, out_mim_rd)

# 1. Raw probabilities from IPD
prob_A_raw <- mean(AC_IPD_binY_contX$y[AC_IPD_binY_contX$trt == "A"])
prob_C_raw_AC <- mean(AC_IPD_binY_contX$y[AC_IPD_binY_contX$trt == "C"])

# 2. Raw probabilities from ALD
prob_B_raw <- BC_ALD_binY_contX |> 
  filter(variable == "y", statistic == "mean", trt == "B") |> 
  pull(value)
prob_C_raw_BC <- BC_ALD_binY_contX |> 
  filter(variable == "y", statistic == "mean", trt == "C") |> 
  pull(value)

# 3. Naive Anchored Risk Difference
naive_anchored_rd <- (prob_A_raw - prob_C_raw_AC) - (prob_B_raw - prob_C_raw_BC)

# 4. Extract and bind results
adjusted_results <- data.frame(
  Estimate = c(out_maic_rd$results$contrasts$means$AB,
               out_gcomp_ml_rd$results$contrasts$means$AB), 
  lower.0.95 = c(out_maic_rd$results$contrasts$CI$AB[1],
                 out_gcomp_ml_rd$results$contrasts$CI$AB[1]), 
  upper.0.95 = c(out_maic_rd$results$contrasts$CI$AB[2],
                 out_gcomp_ml_rd$results$contrasts$CI$AB[2]), 
  Method = c("MAIC (Anchored, Adjusted)", "G-comp ML (Anchored, Adjusted)"))

naive_result <- data.frame(
  Treatments = "AB",
  Estimate = naive_anchored_rd,
  `Std. Error` = NA, lower.0.95 = NA, upper.0.95 = NA,
  Method = "Naive Indirect (Anchored, Unadjusted)",
  check.names = FALSE
)

all_results_intro <- bind_rows(adjusted_results, naive_result)

# Order the factor for plotting
all_results_intro$Method <- factor(all_results_intro$Method, levels = c(
  "Naive Indirect (Anchored, Unadjusted)",
  "MAIC (Anchored, Adjusted)",
  "G-comp ML (Anchored, Adjusted)"
))

# Builds the forest plot
ggplot(all_results_intro, aes(x = Estimate, y = Method, color = Method)) +
  geom_point(size = 4) +
  geom_errorbarh(aes(xmin = lower.0.95, xmax = upper.0.95), height = 0.2, linewidth = 1, na.rm = TRUE) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Relative Treatment Effect (A vs B)",
    subtitle = "Adjusted vs. Unadjusted Anchored Estimates",
    x = "Risk Difference",
    y = ""
  ) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 11, face = "bold"),
    plot.title = element_text(face = "bold", size = 14)
  )

# Advanced Applications: Mixed Data and Custom Distributions
data(AC_IPD_contY_mixedX)
data(BC_ALD_contY_mixedX)

lin_form_mixed <- as.formula("y ~ X1 + X2 + X3 + trt + trt:(X1 + X2 + X4)")
bal_form_mixed <- as.formula("~ X1 + X2 + X4")

formula_list_mixed <- list(outcome_model = lin_form_mixed, balance_model = bal_form_mixed)

custom_distns <- c(X1 = "gamma", X2 = "binom")

out_gcomp_custom <- outstandR(
  ipd_trial = AC_IPD_contY_mixedX,
  ald_trial = BC_ALD_contY_mixedX,
  strategy = strategy_gcomp_ml(
    formula = formula_list_mixed,
    family = gaussian(link = "identity"),
    marginal_distns = custom_distns,
    N = 1000
  ),
  seed = 12345
)

print(out_gcomp_custom)

# Bayesian G-Computation
out_gcomp_bayes <- outstandR(
  ipd_trial = AC_IPD_contY_mixedX,
  ald_trial = BC_ALD_contY_mixedX,
  strategy = strategy_gcomp_bayes(
    formula = formula_list_mixed,
    family = gaussian(link = "identity")
  ),
  seed = 12345
)

# Multiple Imputation Marginalization (MIM)
out_mim <- outstandR(
  ipd_trial = AC_IPD_contY_mixedX,
  ald_trial = BC_ALD_contY_mixedX,
  strategy = strategy_mim(
    formula = formula_list_mixed,
    family = gaussian(link = "identity")
  ),
  seed = 12345
)

# Model Diagnostics, Misspecification & Variable Selection
cat("MAIC Effective Sample Size:", out_maic$model$ESS, "\n")

ggplot(data.frame(weights = out_maic$model$weights), aes(x = weights)) +
  geom_histogram(binwidth = 0.1, color = "black", fill = "lightgray") +
  labs(title = "Distribution of MAIC Weights", x = "Weight", y = "Frequency") +
  theme_minimal()

bal_form_minimal <- as.formula("~ EM_cont_1")

out_maic_minimal <- outstandR(
  ipd_trial = AC_IPD_binY_contX, 
  ald_trial = BC_ALD_binY_contX,
  strategy = strategy_maic(
    formula = list(
      outcome_model = formula("y ~ trt"),
      balance_model = bal_form_minimal
    ),
    family = binomial(link = "logit") 
  ),
  seed = 12345
)

cat("Minimal Model ESS:", out_maic_minimal$model$ESS, "\n")
cat("Minimal Model Std. Error (A vs B):", out_maic_minimal$results$contrasts$variances$AB, "\n")

# Kitchen sink MAIC
bal_form_maximal <- as.formula("~ PF_cont_1 + PF_cont_2 + EM_cont_1 + EM_cont_2")

out_maic_maximal <- outstandR(
  ipd_trial = AC_IPD_binY_contX, 
  ald_trial = BC_ALD_binY_contX,
  strategy = strategy_maic(
    formula = list(
      outcome_model = formula("y ~ trt"),
      balance_model = bal_form_maximal
    ),
    family = binomial(link = "logit") 
  ),
  seed = 12345
)

cat("Maximal Model ESS:", out_maic_maximal$model$ESS, "\n")
cat("Maximal Model Std. Error (A vs B):", out_maic_maximal$results$contrasts$variances$AB, "\n")

# Unanchored Comparisons and Visualization
full_form <- as.formula("y ~ PF_cont_1 + PF_cont_2 + EM_cont_1 + EM_cont_2 + trt")

out_gcomp_unanchored <- outstandR(
  ipd_trial = AC_IPD_binY_contX, 
  ald_trial = BC_ALD_binY_contX,
  strategy = strategy_gcomp_ml(
    formula = list(outcome_model = full_form),
    family = binomial(link = "logit"),
    N = 1000
  ),
  scale = "risk_difference", # Ensures absolute effects are on the probability scale
  seed = 12345
)

# 1. Extract the marginalized probability for A in the target population
prob_A_unanchored <- out_gcomp_unanchored$results$absolute$means$A 

# 2. Extract raw observed probability for B (from the ALD)
prob_B_raw <- BC_ALD_binY_contX |> 
  filter(variable == "y", statistic == "mean", trt == "B") |> 
  pull(value)

# 3. Calculate the unanchored Risk Difference
unanchored_rd <- prob_A_unanchored - prob_B_raw

cat("Unanchored Risk Difference (A vs B):", round(unanchored_rd, 3), "\n")

# Raw probability from IPD for A
prob_A_raw <- mean(AC_IPD_binY_contX$y[AC_IPD_binY_contX$trt == "A"])

# Naive Unanchored Risk Difference: A - B directly
naive_unanchored_rd <- prob_A_raw - prob_B_raw

adjusted_results <- data.frame(
  Estimate = c(out_maic_rd$results$contrasts$means$AB,
               out_gcomp_ml_rd$results$contrasts$means$AB), 
  lower.0.95 = c(out_maic_rd$results$contrasts$CI$AB[1],
                 out_gcomp_ml_rd$results$contrasts$CI$AB[1]), 
  upper.0.95 = c(out_maic_rd$results$contrasts$CI$AB[2],
                 out_gcomp_ml_rd$results$contrasts$CI$AB[2]), 
  Method = c("MAIC (Anchored, Adjusted)", "G-comp ML (Anchored, Adjusted)"))

manual_results <- data.frame(
  Treatments = rep("AB", 3),
  Estimate = c(unanchored_rd, naive_anchored_rd, naive_unanchored_rd),
  `Std. Error` = NA, 
  lower.0.95 = NA,   
  upper.0.95 = NA,
  Method = c("G-comp ML (Unanchored, Adjusted)", 
             "Naive Indirect (Anchored, Unadjusted)", 
             "Naive Direct (Unanchored, Unadjusted)"),
  check.names = FALSE
)

all_results_final <- bind_rows(adjusted_results, manual_results)

all_results_final$Method <- factor(all_results_final$Method, levels = c(
  "Naive Direct (Unanchored, Unadjusted)",
  "G-comp ML (Unanchored, Adjusted)",
  "Naive Indirect (Anchored, Unadjusted)",
  "MAIC (Anchored, Adjusted)",
  "G-comp ML (Anchored, Adjusted)"
))

ggplot(all_results_final, aes(x = Estimate, y = Method, color = Method)) +
  geom_point(size = 4) +
  geom_errorbarh(aes(xmin = lower.0.95, xmax = upper.0.95), height = 0.2, linewidth = 1, na.rm = TRUE) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Relative Treatment Effect (A vs B)",
    subtitle = "The Impact of Population Adjustment and Anchoring Assumptions",
    x = "Risk Difference",
    y = ""
  ) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 11, face = "bold"),
    plot.title = element_text(face = "bold", size = 14)
  )

