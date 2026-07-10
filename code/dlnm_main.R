# dlnm_main.R
# Main and stratified distributed-lag non-linear models for TB and NTM vs
# SARS-CoV-2 activity (cumulative RR per +1 SD, 0-12 week lag). Reproduces the
# combined Table 3. RUN INSIDE THE CIHI SAE.
#
# Expects weekly analytic frame `dat` (see PSEUDOCODE.md for construction) with:
#   n_tb, n_ntm, week_pop, series_week, sin1,cos1,sin2,cos2,
#   norm_sars, age_cat (0 all / 1 <65 / 2 >=65), ich (0/1)
library(dlnm); library(splines)

cumRR <- function(d, outcome) {
  cb <- crossbasis(d$norm_sars, lag = 12,
                   argvar = list(fun = "lin"),
                   arglag = list(fun = "ns", df = 3))
  m  <- glm(as.formula(paste(outcome,
        "~ cb + sin1 + cos1 + sin2 + cos2 + series_week + offset(log(week_pop))")),
        family = quasipoisson(), data = d)
  pr <- crosspred(cb, m, at = 1, cumul = TRUE)
  c(RR = pr$allRRfit["1"], lo = pr$allRRlow["1"], hi = pr$allRRhigh["1"])
}

strata <- list(
  "Overall"   = quote(TRUE),
  "Age <65 y" = quote(age_cat == 1),
  "Age >=65 y"= quote(age_cat == 2),
  "Non-ICH"   = quote(ich == 0),
  "ICH"       = quote(ich == 1))

for (oc in c("n_tb","n_ntm")) {
  cat("\n===", oc, "cumulative RR per SD SARS-CoV-2 (0-12 wk) ===\n")
  for (nm in names(strata)) {
    d <- subset(dat, eval(strata[[nm]]))
    e <- cumRR(d, oc)
    cat(sprintf("%-11s %.3f (%.3f-%.3f)\n", nm, e[1], e[2], e[3]))
  }
}
# Between-stratum heterogeneity: Wald test on the difference of log cumulative
# RRs (age: <65 vs >=65; immune: non-ICH vs ICH) via the pooled interaction model.
