# attributable_fraction.R
# SARS-CoV-2-attributable fraction of pandemic-era TB by a counterfactual that
# sets SARS-CoV-2 activity to zero, with a parametric-bootstrap 95% CI.
# RUN INSIDE THE CIHI SAE. Reproduces Figure 4 / the 9.2% PAF.
library(dlnm); library(splines); library(MASS)

cb <- crossbasis(dat$norm_sars, lag = 12,
                 argvar = list(fun = "lin"), arglag = list(fun = "ns", df = 3))
m  <- glm(n_tb ~ cb + sin1 + cos1 + sin2 + cos2 + series_week +
            offset(log(week_pop)), family = quasipoisson(), data = dat)

pand <- dat$date >= as.Date("2020-03-15")
fitted_all <- predict(m, type = "response")
# counterfactual design: SARS-CoV-2 crossbasis set to zero
cf <- dat; cf$norm_sars <- 0
cb0 <- crossbasis(cf$norm_sars, lag = 12,
                  argvar = list(fun = "lin"), arglag = list(fun = "ns", df = 3))
Xall <- model.matrix(m); Xcf <- Xall
cb_cols <- grep("^cb", colnames(Xall))
Xcf[, cb_cols] <- attr(cb0, "argvar")  # replace CB columns with zero-exposure basis
lin_cf <- Xcf %*% coef(m)
cf_counts <- exp(lin_cf) * 1               # offset already in linear predictor
paf <- sum((fitted_all - cf_counts)[pand]) / sum(fitted_all[pand])

# parametric bootstrap (1000 multivariate-normal draws from coef covariance)
set.seed(1); B <- 1000
draws <- mvrnorm(B, coef(m), vcov(m))
pafs <- apply(draws, 1, function(b){
  fa <- exp(Xall %*% b); fc <- exp(Xcf %*% b)
  sum((fa - fc)[pand]) / sum(fa[pand]) })
cat(sprintf("PAF = %.1f%% (95%% CI %.1f-%.1f)\n",
            100*paf, 100*quantile(pafs,.025), 100*quantile(pafs,.975)))
