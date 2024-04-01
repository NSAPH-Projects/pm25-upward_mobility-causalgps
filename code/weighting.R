rm(list = ls())
gc()
library(data.table)
library(tidyr)
library(dplyr)
library(parallel)
library(CausalGPS)
library(MASS)
library(foreign)
library(arrow)
library(WeightIt)
library(lme4)
library(npmlreg)
library(HLMdiag)
library(usmap)
library(ggplot2)

atlas_subset <- read_parquet("atlas_subset.parquet")
covariates <- read_parquet("covariates.parquet")

treat = atlas_subset$pm25_1982

trim = 0.95
atlas_subset = subset(atlas_subset, treat <= quantile(treat,trim) & treat >= quantile(treat,1-trim))
covariates = subset(covariates, treat <= quantile(treat,trim) & treat >= quantile(treat,1-trim))

covariates$med_hhinc1990 = log(covariates$med_hhinc1990)
covariates$popdensity2000 = log(covariates$popdensity2000)




# Causal Inference -- weighting
set.seed(5761)
pseudo_pop <- generate_pseudo_pop(
  atlas_subset$kfr_pooled_pooled_p25,
  atlas_subset$pm25_1982,
  covariates,
  ci_appr = "weighting",
  gps_model = "parametric",
  use_cov_transform = T, # changed it to TRUE
  transformers = list("pow2", "pow3"),
  bin_seq = NULL,
  trim_quantiles = c(0, 1), # c(0.05, 0.95), 
  sl_lib = c("m_xgboost"), # c("SL.xgboost", "SL.earth","SL.gam","SL.ranger")
  # use these two lines to specify tuning parameters for the XGBoost model
  params = list(xgb_nrounds=seq(5,30,1),
                  xgb_eta=seq(0.1,0.3,0.01)),
  nthread = 24,
  covar_bl_method = "absolute",
  covar_bl_trs = 0.1,
  covar_bl_trs_type = "mean",
  # increased the # of attempts to achieve covariate balance (the larger 'max_attempt' is, the slower it gets)
  max_attempt = 4) 

summary(pseudo_pop$pseudo_pop$counter_weight)
plot(pseudo_pop)
pseudo_pop$adjusted_corr_results

cap = .995
qw = quantile(pseudo_pop$pseudo_pop$counter_weight, probs = cap)
q1w = quantile(pseudo_pop$pseudo_pop$counter_weight, probs = 1 - cap)
counter_weight_cap = pseudo_pop$pseudo_pop$counter_weight
counter_weight_cap[pseudo_pop$pseudo_pop$counter_weight >= qw] = qw
counter_weight_cap[pseudo_pop$pseudo_pop$counter_weight <= q1w] = q1w

pseudo_pop$adjusted_corr_results = absolute_weighted_corr_fun(w = data.table(pseudo_pop$pseudo_pop$w),
                           vw = data.table(counter_weight_cap),
                           c = data.table(pseudo_pop$pseudo_pop[,-c(1:5)]))

summary(counter_weight_cap)
(sum(counter_weight_cap)^2) / sum(counter_weight_cap^2)
pseudo_pop$adjusted_corr_results



# analysis stage
# atlas_subset$zkfr_pooled_pooled_p25 = c(scale(atlas_subset$kfr_pooled_pooled_p25))
atlas_subset$zpm25_1982 = as.numeric(scale(atlas_subset$pm25_1982, center = T, scale = F))
atlas_subset$lkfr_pooled_pooled_p25 = log(atlas_subset$kfr_pooled_pooled_p25)
no_counties = length(unique(atlas_subset$county_state))

lm0 = lm(lkfr_pooled_pooled_p25 ~ pm25_1982, data = atlas_subset)
lmw = lm(lkfr_pooled_pooled_p25 ~ zpm25_1982, weights = counter_weight_cap, data = atlas_subset)
rlmw = rlm(lkfr_pooled_pooled_p25 ~ zpm25_1982, weights = counter_weight_cap, data = atlas_subset)

# random effect models
library(lme4)
library(npmlreg)
lmi = lmer(formula = lkfr_pooled_pooled_p25 ~ zpm25_1982 + (1 | county_state), data = atlas_subset, weights = counter_weight_cap)
lms = lmer(formula = lkfr_pooled_pooled_p25 ~ zpm25_1982 + (1 + zpm25_1982 | county_state), data = atlas_subset, weights = counter_weight_cap)

# NPML
K = 1:11 # number of finite mixture components
tmp.npi = lapply(K, function(k) 
  allvc(lkfr_pooled_pooled_p25 ~ zpm25_1982, random = ~1|county_state, data = atlas_subset, random.distribution = "np", k = k, weights = counter_weight_cap, lambda = 0.99, plot.opt = 0))
BIC = sapply(1:length(K), function(k) tmp.npi[[k]]$disparity + 1e250 * any(tmp.npi[[k]]$masses < 0.01)) + log(no_counties) * (1 + K + K - 1 + K)
BIC[which(BIC == 1e250):length(BIC)] = 1e250
ICL = BIC - 2 * sapply(1:length(K), function(k) sum(tmp.npi[[k]]$post.prob * ifelse(tmp.npi[[k]]$post.prob > 0, log(tmp.npi[[k]]$post.prob), 0)))

BIC
which.min(BIC)
npi = tmp.npi[[which.min(BIC)]]
rm(tmp.npi)

tmp.nps = lapply(K, function(k) 
  allvc(lkfr_pooled_pooled_p25 ~ zpm25_1982, random = ~zpm25_1982|county_state, data = atlas_subset, random.distribution = "np", k = k, weights = counter_weight_cap, lambda = 0.99, plot.opt = 0))
BIC = sapply(1:length(K), function(k) tmp.nps[[k]]$disparity + 1e250 * any(tmp.nps[[k]]$masses < 0.01)) + log(no_counties) * (K + K + K - 1 + K)
BIC[which(BIC == 1e250):length(BIC)] = 1e250
ICL = BIC - 2 * sapply(1:length(K), function(k) sum(tmp.nps[[k]]$post.prob * ifelse(tmp.nps[[k]]$post.prob > 0, log(tmp.nps[[k]]$post.prob), 0)))

BIC
which.min(BIC)
nps = tmp.nps[[which.min(BIC)]]
rm(tmp.nps)
