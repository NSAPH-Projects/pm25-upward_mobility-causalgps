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

# compare the goodness of fit of the two models (random intercepts vs random intercepts and slope)
anova(lmi, lms)
nc <- detectCores()
cl <- makeCluster(rep("localhost", nc))
pbkrtest::PBmodcomp(largeModel = lms, smallModel = lmi, nsim = 50, cl = cl, seed = 101)

# US Map plot (just an example to see the results)
dta.map = data.frame(fips = rownames(ranef(lms)$county_state), value = ranef(lms)$county_state[,1])
plot_usmap(data = dta.map, values = "value") + 
  scale_fill_continuous(
    low = "blue", high = "red", name = "Random intercepts", label = scales::comma
  ) + 
  theme(legend.title=element_text(size=18), legend.text = element_text(size=15), legend.position = "right")

# US Map plot (random intercepts and slopes quartiles)
qnle = 4
county.slope = percent_change((fixef(lms)[2] + est$zpm25_1982)) # 100 * (fixef(lms)[2] + est$zpm25_1982)
dta.map = data.frame(fips = rownames(ranef(lms)$county_state), value = cut(county.slope, quantile(county.slope, 0:qnle/qnle), dig.lab = 2))
plot_usmap(data = dta.map, values = "value", exclude = c("AK","HI")) + 
  scale_fill_manual(values = c(colorRampPalette(c("blue", "red"))(qnle)), name = "", guide = guide_legend(reverse = TRUE), na.translate = F) +
  theme(legend.title=element_text(size=18), legend.text = element_text(size=30))

# NPML
K = 1:10 # number of finite mixture components
tmp.npi = lapply(K, function(k) 
  allvc(lkfr_pooled_pooled_p25 ~ zpm25_1982, random = ~1|county_state, data = atlas_subset, random.distribution = "np", k = k, weights = mod$weights_cap, lambda = 0.99, plot.opt = 0))
BIC = sapply(1:length(K), function(k) tmp.npi[[k]]$disparity + 1e250 * any(tmp.npi[[k]]$masses < 0.01)) + log(no_counties) * (1 + K + K - 1 + K)
ICL = BIC - 2 * sapply(1:length(K), function(k) sum(tmp.npi[[k]]$post.prob * ifelse(tmp.npi[[k]]$post.prob > 0, log(tmp.npi[[k]]$post.prob), 0)))

BIC
which.min(BIC)
npi = tmp.npi[[which.min(BIC)]]

tmp.nps = lapply(K, function(k) 
  allvc(lkfr_pooled_pooled_p25 ~ zpm25_1982, random = ~zpm25_1982|county_state, data = atlas_subset, random.distribution = "np", k = k, weights = mod$weights_cap, lambda = 0.99, plot.opt = 0))
BIC = sapply(1:length(K), function(k) tmp.nps[[k]]$disparity + 1e250 * any(tmp.nps[[k]]$masses < 0.01)) + log(no_counties) * (K + K + K - 1 + K)
ICL = BIC - 2 * sapply(1:length(K), function(k) sum(tmp.nps[[k]]$post.prob * ifelse(tmp.nps[[k]]$post.prob > 0, log(tmp.nps[[k]]$post.prob), 0)))

BIC
which.min(BIC)
nps = tmp.nps[[which.min(BIC)]] 

# US Map plot (just an example to see the results)
dta.map = aggregate(apply(nps$post.prob, 1, which.max), by = list(nps$data$county_state), mean)
colnames(dta.map) = c("fips", "value")
plot_usmap(data = dta.map, values = "value") + 
  scale_fill_continuous(
    low = "blue", high = "red", name = "Class posterior probability", label = scales::comma
  ) + 
  theme(legend.title=element_text(size=18), legend.text = element_text(size=15), legend.position = "right")
