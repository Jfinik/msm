#################################################################
# Code for Lab 2; Question 3 [IPTW/MSM time-varying treatment]  #
# 0. Generate SWIG for the study design                         #
# 1. Estimate IPTW weights                                      #
# 2. Fit weighted MSM outcome model & predict counterfactuals   #
# 3. Compute avg. expected counterfactuals for regimes          # 
# 4. Contrast regimes of interest                               #
# 5. Use bootstrapping to get confidence intervals              #
#################################################################

# Load libraries
library(WeightIt)
library(cobalt)
library(marginaleffects)
library(boot)
library(tidyverse)
library(tinytex)
library(quickdag)

#################################################################

# [0] Create SWIG
# Encode paths
edges <- c("A1 -> {A2 S1 A3 Y}",
           "A2 -> {A3 S2 Y}",
           "A3 -> Y",
           "S1 -> {A2 S2}",
           "S2 -> {A3, Y}",
           "U -> {S1 S2 Y}")
sep.point.size = 15

# DAG Input
dag  <- qd_dag(edges, verbose = FALSE)
dag %>% render_graph()
# Generate SWIG
swig <- dag %>%
  qd_swig(fixed.nodes = c("A1", "A2", "A3"))
swig %>% render_graph()

#################################################################

# Load data
comp <- read.csv("C:/Users/jfinik/Desktop/comp_exam_biostat_question2_final.csv")
head(comp)

# Study Design: RCT; 2 parallel arms [treatment, placebo]
# Y is a continuous outcome [blood pressure]
# S1 and S2 [AE presence] are time varying covariates measured after the first and second treatments
# A1, A2, A3 [treatment] are indicators of time varying treatment adherence; three timepoints [all adhered to treatment at A1]

# Objective: "estimate the average causal effect (and 95%CI) of taking all three doses of the new drug compared to placebo"
# We need to compare counterfactual Ys for treated [A2=1 & A3=1] vs. not [A2=0 & A3 =0] 

# Estimand: the average treatment effect in the treated [ATT = diff btw avg. Ys for the treated &  avg Ys they would have had, had they NOT been treated]
# Effect: marginal treatment effect in the sub-population that received treatment; so we subset our data to arm = 1 for treated
# Effect Measure: our outcome is continuous [blood pressure] so effect will be measured by mean difference [which is collapsible]

# Subset data to arm = 1 for treated
comp = comp[comp$arm==1,] 

# Using the cobalt package to glimpse data imbalance
bal.tab(list(A2 ~ S1,
             A3 ~ S1 + S2 + A2),
        data = comp, stats = c("m", "ks"),
        which.time = .all)

# Use weightitMSM() function to specify IPTW weights
# We will use "glm" and stabilize = TRUE to get stabilized weights
# List = formulas for each time point [syntax: 'time-specific treatment variable'  ~ 'pre-treatment covariates']
# Formulas must be in temporal order, and must contain all covariates for that time point and all preceding timepoints  
# Interactions and functions of covariates are allowed [include interactions to minimize misspecification]

# [1] Weighting [note this function generates weights for each timepoint & multiplies across timepoints]
Wmsm.out <- weightitMSM(list(A2 ~ S1,
                             A3 ~ A2*S1*S2),
                        data = comp, method = "glm",
                        stabilize = TRUE)
                        Wmsm.out
# Check mean of weights
summary(Wmsm.out)

# Attach the weights to the data
comp$msm_weights <- Wmsm.out$weights

# [2] Fit a MSM outcome model using weighted lm since our outcome is continuous
fit <- lm(Y ~ A2 * A3, data = comp, weights = msm_weights)
# Note confounding was addressed via weighting so we don't have to include confounders in this model

# [3] Compute avg. expected counterfactuals from predicted values of our MSM using marginaleffects::avg_predictions() 
(p <- avg_predictions(fit,
                      vcov = "HC3",
                      newdata = datagridcf(A2 = 0:1, A3 = 0:1),
                      by = c("A2", "A3"),
                      type = "response"))
# This outputs the expected counterfactual estimates for each regime [0,0] [0,1] [1,0] [1,1]
# For always treated A2=1 & A3=1; mean Y for [1,1]: M=155 (95%CI: 154, 156) 
# For never treated A2=0 & A3=0; mean Y for [0,0]: M=165 (95%CI: 164, 167)
# Same estimates as the 'by hand' coding from the lab
# Note results output "robust" sandwich SEs; 95% of the time estimates will be in the CI derived this way 

# [4] Directly compare these expected counterfactuals to get the ATT using hypotheses()
# To specify individual regimes, we identify the rows of the avg_predictions() output 
# E.g. we want to compare the regime of always treated [1,1] vs. never treated [0,0] indexed by b
# 4th row [b4] corresponds to always treated [1,1]; 1st row [b1] corresponds to never treated [0,0] 
hypotheses(p, "b4 - b1 = 0")

# We can also output all contrasts of expected counterfactuals between regimes using hypotheses(="pairwise")
# (p <- avg_predictions(fit,
#                      vcov = "HC3",
#                      newdata = datagridcf(A2 = 0:1, A3 = 0:1),
#                      by = c("A2", "A3"),
#                      type = "response",
#                      hypothesis = "pairwise"))

# [5] Bootstrap to get precise confidence intervals
boot_fun <- function(data, i) {
  boot_data <- data[i,]
  
  # Weighting
  Wmsm.out <- weightitMSM(list(A2 ~ S1,
                               A3 ~ A2*S1*S2),
                          data = comp, method = "glm",
                          stabilize = TRUE)

  # Bring weights into the dataset
  boot_data$msm_weights <- Wmsm.out$weights
  
  #Fit outcome model
  fit <- lm(Y ~ A2 * A3, data = boot_data, weights = msm_weights)
  
  # G-computation
  p <- avg_predictions(fit,
                        vcov = "HC3",
                        newdata = subset(boot_data, A2 = 0:1, A3 = 0:1),
                        by = c("A2", "A3"),
                       type = "response")
  p$estimate
}

# Run bootstrap with R=199, would increase this for a real analysis
set.seed(54321)
boot_out <- boot(comp, boot_fun, R = 199)
boot_out

# Compute confidence interval for always treated [1,1]
at_m <- 155
at_se <- 0.4949
at_CI <- c(at_m - 1.96 * se, at_m + 1.96 * se)
at_CI
# Final estimate for always treated: M=155 (95%CI: 154, 156)

# Compute confidence interval for never treated [0,0]
nt_m <- 165
nt_se <- 0.6439
nt_CI <- c(nt_m - 1.96 * se, nt_m + 1.96 * se)
nt_CI
# Final estimate for never treated: M=165 (164, 166)

# CI via bootstrap  nearly identical to 'robust sandwich' SE in this case
# Bootstrap is preferred when using weighting methods 

######################################################
# Method 2: TMLE with cross-fitting via lmtp package #         
######################################################

# TMLE with cross-fitting (lmtp_tmle) using SuperLearner
# lmtp_tmle() uses candidate learners appropriate for outcome type
# To make your own custom stack of learners, make_learner_stack()

# Load libraries
library(lmtp)
library(SuperLearner)
library(earth)
library(caret)
library(future)
library(ranger)
library(twang)

# Recall data has been subsetted to arm=1
head(comp)

# Encode variables in lmtp syntax as follows
  A <- c("A2", "A3")
  L <- list(c("S1"), c("S2"))
  Y <- "Y"

# Always treated
  y1 <- lmtp_tmle(comp, A, Y, time_vary = L, shift = static_binary_on, 
            intervention_type = "static", outcome_type = "continuous",folds = 5)
# Never treated
  y0 <- lmtp_tmle(comp, A, Y, time_vary = L, shift = static_binary_off, 
            intervention_type = "static", outcome_type = "continuous",folds = 5)

# Note folds =  # of cross-validation folds [V] for Superlearner]
  
# Estimate the causal mean difference
  lmtp_contrast(y1, ref=y0, type = "additive")
  
# Same results

#################################################################
# Assumptions for causal inference [time varying settings]      #
# Sequential exchangeability                                    #
# Positivity [enough overlap]                                   #
# Consistency/SUTVA                                             #
# No measurement error                                          #
# No model misspecification                                     #
#################################################################
  
