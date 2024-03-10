rm(list = ls())
gc()
graphics.off()
# E values
# library(EValue)
library(arrow)
source("utils.R")

atlas_subset <- read_parquet("atlas_subset.parquet")
covariates <- read_parquet("covariates.parquet")

treat = atlas_subset$pm25_1982

trim = 0.95
atlas_subset = subset(atlas_subset, treat <= quantile(treat,trim) & treat >= quantile(treat,1-trim))
covariates = subset(covariates, treat <= quantile(treat,trim) & treat >= quantile(treat,1-trim))
atlas_subset$lkfr_pooled_pooled_p25 = log(atlas_subset$kfr_pooled_pooled_p25)
sdY = sd(atlas_subset$lkfr_pooled_pooled_p25)


# load the results for entropy balance, weighting and matching
load("results_entropy.RData")
res_e = res_models
load("results_weighting.RData")
res_w = res_models
load("results_matching.RData")
res_m = res_models

mod = names(res_e) # models list
delta = 1 # one-unit increase in the exposure
true = 0

evalues = matrix(NA, 3, length(mod))
rownames(evalues) = c("Entropy Balance", "IPTW", "GPS Matching")
colnames(evalues) = mod
for(j in 1:length(mod)) {
  idx = ifelse(j <= 5, 2, 1)
  evalues[1,j] = evalue.OLS(est = res_e[[j]][idx,1], se = res_e[[j]][idx,2], sd = sdY, delta = delta, true = true)[2,1]
  evalues[2,j] = evalue.OLS(est = res_w[[j]][idx,1], se = res_w[[j]][idx,2], sd = sdY, delta = delta, true = true)[2,1]
  evalues[3,j] = evalue.OLS(est = res_m[[j]][idx,1], se = res_m[[j]][idx,2], sd = sdY, delta = delta, true = true)[2,1]
}
evalues

xtable(evalues, digits = 3)


# DO NOT RUN
# GPS adjustment
# evalue.OLS(est = lm_gps$coefficients[2], se = summary(lm_gps)$coefficients['w', 'Std. Error'], sd = sd(lm_gps$residuals), delta = 1, true = 0)

# GPS weighting
# evalue.OLS(est = lm_wgps$coefficients[2], se = summary(lm_wgps)$coefficients['w', 'Std. Error'], sd = sd(lm_wgps$residuals), delta = 1, true = 0)
