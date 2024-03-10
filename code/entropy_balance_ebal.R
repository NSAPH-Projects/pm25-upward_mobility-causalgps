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

# entropy weights (Vegetabile et al., 2021) + cap at the 995th percentile
mod = weightit(formula = atlas_subset$pm25_1982 ~ covariates, method = "ebal", d.moments = 3, moments = 2, base.weight = NULL)
# mod = weightit(formula = atlas_subset$pm25_1982 ~ covariates, method = "ebal", moments = 2)
cap = .995
mod$weights_cap <- mod$weights
mod$weights_cap[mod$weights >= quantile(mod$weights, cap)] = quantile(mod$weights, cap)
mod$weights_cap[mod$weights <= quantile(mod$weights, 1 - cap)] = quantile(mod$weights, 1 - cap)
mod = weightit(formula = atlas_subset$pm25_1982 ~ covariates, method = "ebal", d.moments = 3, moments = 2, base.weight = mod$weights_cap)
mod$weights_cap <- mod$weights
mod$weights_cap[mod$weights >= quantile(mod$weights, cap)] = quantile(mod$weights, cap)
mod$weights_cap[mod$weights <= quantile(mod$weights, 1 - cap)] = quantile(mod$weights, 1 - cap)
summary(mod$weights)
summary(mod$weights_cap)

##################
# analysis stage
atlas_subset$zpm25_1982 = as.numeric(scale(atlas_subset$pm25_1982, center = T, scale = F))
atlas_subset$lkfr_pooled_pooled_p25 = log(atlas_subset$kfr_pooled_pooled_p25)
no_counties = length(unique(atlas_subset$county_state))

lm0 = lm(lkfr_pooled_pooled_p25 ~ zpm25_1982, data = atlas_subset)
lmw = lm(lkfr_pooled_pooled_p25 ~ zpm25_1982, weights = mod$weights_cap, data = atlas_subset)
rlmw = rlm(lkfr_pooled_pooled_p25 ~ zpm25_1982, weights = mod$weights_cap, data = atlas_subset)

# random effect models
lmi = lmer(formula = lkfr_pooled_pooled_p25 ~ zpm25_1982 + (1 | county_state), data = atlas_subset, weights = mod$weights_cap)
lms = lmer(formula = lkfr_pooled_pooled_p25 ~ zpm25_1982 + (1 + zpm25_1982 | county_state), data = atlas_subset, weights = mod$weights_cap) 
# if you see any warning add the following command:
# control = lmerControl(optCtrl=list(maxfun=2e5), optimizer = "bobyqa")
