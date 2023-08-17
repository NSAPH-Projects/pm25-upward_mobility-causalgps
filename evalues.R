# E values
# library(EValue)
source("utils.R")

# GPS adjustment
evalue.OLS(est = lm_gps$coefficients[2], se = summary(lm_gps)$coefficients['w', 'Std. Error'], sd = sd(lm_gps$residuals), delta = 1, true = 0)

# GPS weighting
evalue.OLS(est = lm_wgps$coefficients[2], se = summary(lm_wgps)$coefficients['w', 'Std. Error'], sd = sd(lm_wgps$residuals), delta = 1, true = 0)
