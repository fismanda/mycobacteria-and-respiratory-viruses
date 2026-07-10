# negative_controls_table2.R
# Recompute EVERY cell of Table 2 (influenza A and RSV negative controls vs TB),
# so the values are verified end-to-end rather than transcribed.
# RUN INSIDE THE CIHI SAE. No row-level data leaves the enclave; only the
# aggregate coefficient table below is exported.
#
# Expects a gap-free weekly analytic data frame `dat` with columns:
#   n_tb        weekly incident TB count
#   week_pop    weekly population (offset)
#   series_week integer week index (1..676)
#   date        week start date (for period subsetting)
#   sin1,cos1,sin2,cos2   first/second-order Fourier terms (period = 52.18 wk)
#   norm_flu, norm_rsv, norm_sars   SD-standardized weekly viral activity
library(dlnm); library(splines)

PANDEMIC_START <- as.Date("2020-03-15")
PRE_START      <- as.Date("2013-01-01"); PRE_END <- as.Date("2019-12-31")

## generic fit: cumulative RR per +1 SD of `exp` over a 0..LAG-week window,
## optionally adjusting for a second virus crossbasis (adjust), optional trend.
fit_cell <- function(d, exp, lag = 12, adjust = NULL, trend = TRUE) {
  if (nrow(d) < 30 || sum(d$n_tb) < 40) return(c(RR=NA, lo=NA, hi=NA))
  cb <- crossbasis(d[[exp]], lag = lag,
                   argvar = list(fun = "lin"),
                   arglag = list(fun = "ns", df = if (lag > 12) 4 else 3))
  rhs <- "cb + sin1 + cos1 + sin2 + cos2 + offset(log(week_pop))"
  if (trend) rhs <- paste(rhs, "+ series_week")
  cbA <- NULL
  if (!is.null(adjust)) {
    cbA <- crossbasis(d[[adjust]], lag = lag,
                      argvar = list(fun = "lin"),
                      arglag = list(fun = "ns", df = if (lag > 12) 4 else 3))
    rhs <- paste("cbA +", rhs)
  }
  m  <- glm(as.formula(paste("n_tb ~", rhs)), family = quasipoisson(), data = d)
  pr <- crosspred(cb, m, at = 1, cumul = TRUE)
  c(RR = pr$allRRfit["1"], lo = pr$allRRlow["1"], hi = pr$allRRhigh["1"])
}

pre  <- subset(dat, date >= PRE_START & date <= PRE_END)
pand <- subset(dat, date >= PANDEMIC_START)

rows <- list(
  # label,                                    virus,      subset, lag, adjust,    trend
  c("Influenza | Overall (2013-2024)",         "norm_flu","all", 12,  NA,        TRUE),
  c("Influenza | Pre-pandemic (2013-2019)",    "norm_flu","pre", 12,  NA,        TRUE),
  c("Influenza | Pandemic (2020-2024)",        "norm_flu","pand",12,  NA,        FALSE),
  c("Influenza | Pandemic, SARS-adjusted",     "norm_flu","pand",12,  "norm_sars",FALSE),
  c("Influenza | Overall, 24-week lag",        "norm_flu","all", 24,  NA,        TRUE),
  c("RSV | Overall (2013-2024)",               "norm_rsv","all", 12,  NA,        TRUE),
  c("RSV | Pre-pandemic (2013-2019)",          "norm_rsv","pre", 12,  NA,        TRUE),
  c("RSV | Pandemic (2020-2024)",              "norm_rsv","pand",12,  NA,        FALSE),
  c("RSV | Overall, SARS-adjusted",            "norm_rsv","all", 12,  "norm_sars",TRUE),
  c("RSV | Overall, 24-week lag",              "norm_rsv","all", 24,  NA,        TRUE)
)
pick <- function(s) if (s=="pre") pre else if (s=="pand") pand else dat

out <- do.call(rbind, lapply(rows, function(r){
  d <- pick(r[4-1]); est <- fit_cell(d, r[2], lag=as.numeric(r[4]),
        adjust=if (is.na(r[5])) NULL else r[5], trend=as.logical(r[6]))
  data.frame(cell=r[1], RR=round(est[1],3), lo=round(est[2],3), hi=round(est[3],3),
             row.names=NULL)
}))
out$CI <- sprintf("%.3f-%.3f", out$lo, out$hi)
print(out[,c("cell","RR","CI")], right=FALSE)
write.csv(out, "table2_negative_controls_verified.csv", row.names = FALSE)
# SANITY: influenza and RSV pre-pandemic RRs should NOT be identical.
cat("\nCheck: flu vs RSV pre-pandemic identical? ",
    isTRUE(all.equal(out$RR[2], out$RR[7])), " (should be FALSE)\n")
