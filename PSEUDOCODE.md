# SAE-only data assembly (pseudocode)

Row-level data live inside the CIHI Secure Analytic Environment (SAE) and cannot
leave it, so the steps that build the weekly analytic frame `weekly` are given here
as pseudocode. Only aggregate model coefficients are exported. Ethics: University
of Toronto REB #41690.

## 1. Case ascertainment (Stata/R inside SAE)

- Source: linked CIHI Discharge Abstract Database (DAD, inpatient) + National
  Ambulatory Care Reporting System (NACRS, unscheduled ED).
- TB = ICD-10-CA A15, A17, A18, A19 in any diagnosis field (exclude A16, clinically
  suspected but unconfirmed pulmonary TB). NTM = A31.
- Link on encrypted health-card number; keep the earliest qualifying record per
  person per outcome. A person may contribute once to TB and once to NTM.
  Yields 3,985 first TB and 5,214 first NTM acute-care presentations.
- Immunocompromising host condition (ICH): HIV/AIDS (B20-B24, Z21), solid-organ
  transplant (Z940-Z944), stem-cell/marrow transplant (Z941, Z948), haematologic
  malignancy (C81-C85, C88, C90-C96), primary immunodeficiency (D80-D84).
  Medication-induced immunosuppression is not captured; the ICH-negative stratum
  therefore includes people with unmeasured immune compromise.
- Restrict to the 14 GTHA + central-Ontario public health units.
- Aggregate to weekly counts by (week_date, age_cat, ich): tb_count, ntm_count.

## 2. Population offset

- Weekly population by linear interpolation of annual PHU denominators.
- Enters every model as `offset(log(weekly_pop))`.

## 3. Respiratory viral exposures

- SARS-CoV-2: test-adjusted case counts from CCM (Mar 2020-Aug 2022; correction for
  differential testing by age/sex, Fisman 2021 / Bosco 2025), bridged to RVDSS
  percent positivity (Sep 2022-Mar 2024). Each component is scaled before
  concatenation.
- Influenza A/B and RSV: RVDSS for the full period.
- Five further groups from RVDSS, available Sep 2015 onward: human parainfluenza
  virus, human metapneumovirus, adenovirus, enterovirus/rhinovirus, seasonal
  coronaviruses.
- Divide each series by its SD over the period in which its surveillance metric was
  used, without mean-centring, so that zero continues to denote no detected activity.

## 4. Frame assembly

- Index exposure, outcome and population series on calendar week and merge with
  one-to-one date keys.
- Before modelling, assert: no duplicate week keys; no unmatched weeks; no implicit
  padding of any series to the length of another. Every model must be fit on rows
  where all of outcome, exposure and population are observed.
- Analytic period: 1 April 2011 through 31 March 2024, **676 analysed weeks**.
  The 468 weeks through 14 March 2020 are pre-pandemic; the pandemic era begins
  15 March 2020.
- Twelve weeks preceding 1 April 2011 are retained in the frame **only** to supply
  lagged exposure values, and contribute no outcome observations. Exploratory series
  beginning Sep 2015 comprise 433 analysed weeks.

## 5. Seasonality and trend

- Fourier terms sin1, cos1, sin2, cos2 at period 52.18 weeks, in every model.
- Secular trend, by analysis:
  - DLNM: natural spline in continuous week, `ns(t, 13)` over the 13-year series.
  - Level shift: three alternatives, `1` (none), `t` (linear), `ns(t, 3)`.

## 6. Models

### 6a. Viral activity and mycobacterial presentations (DLNM)

    crossbasis(exposure, lag = 12,
               argvar = list(fun = "lin"),
               arglag = list(fun = "ns", df = 3))

    glm(outcome ~ cb + sin1 + sin2 + cos1 + cos2 + ns(t, 13) +
                  offset(log(weekly_pop)),
        family = quasipoisson)

- Estimand: cumulative RR per 1-SD increase over the full lag window, from
  `crosspred(..., cumul = TRUE)`.
- Repeat within age (<65, >=65) and ICH strata; test between-stratum heterogeneity
  by Wald test on the difference in log cumulative RR.
- All non-SARS-CoV-2 exposures additionally adjusted for concurrent SARS-CoV-2
  activity.
- Sensitivity: repeat for SARS-CoV-2 with `lag = 26` and `lag = 52`. Note that at
  lag 52 the crossbasis competes directly with `ns(t, 13)`, whose knots fall at
  approximately annual intervals; report the minimum detectable effect alongside
  the estimate.

### 6b. Pandemic-era level shift, two families

Family 1 — fit pre-pandemic only, project forward:

    m   <- glm.nb(outcome ~ TREND + sin1 + sin2 + cos1 + cos2 +
                            offset(log(weekly_pop)), data = pre)
    post$expct <- predict(m, newdata = post, type = "response")
    r   <- glm.nb(outcome ~ 1 + offset(log(expct)), data = post)
    RR  <- exp(coef(r)[[1]])

  With the log expectation as an offset, `exp(intercept)` is the observed-to-expected
  rate ratio, and equals `sum(observed)/sum(expected)` exactly. Its interval
  conditions on the projection as known and therefore understates uncertainty.

Family 2 — fit an indicator across the full series:

    glm.nb(outcome ~ pandemic + TREND + sin1 + sin2 + cos1 + cos2 +
                     offset(log(weekly_pop)), data = weekly)

  `pandemic` must be integer 0/1, not logical, so that counterfactual prediction
  with `pandemic <- 0` behaves. Note that a flexible TREND fit across the whole
  series can absorb the level change, leaving indicator and trend competing for the
  same variation; the indicator coefficient is then not separately identified and
  should not be used to derive an implied excess.

Run each family under all three TREND specifications, giving six estimates per
outcome. Compare by AIC and BIC **within** a family only — the two families are fit
to different numbers of observations. Report all six.

### 6c. Accumulation of excess

    excess_by_quarter <- tapply(post$outcome - post$expct, post$qtr, sum)
    cumulative        <- cumsum(excess_by_quarter)

- Repeat for each of the three projections; the envelope across them is the grey band
  of Figure 2B.
- Implied excess under a constant rate ratio is `expct * (RR - 1)`, accumulated the
  same way, with each RR applied to its own projection.

## 7. Outputs

Forest and figure scripts run outside the SAE on exported aggregate estimates.

| Script | Produces |
|---|---|
| `dlnm_main.R` | Supplementary Tables S1-S3, Figure 3 estimates |
| `level_shift.R` | Table 2, Figure 2 projections |
| `figure2.R` | Figure 2, panels A and B |
| `forest.R` | Figure 3 |

**Note.** `attributable_fraction.R` in earlier releases corresponds to an analysis
no longer reported: an attributable fraction is not defined for a null exposure
association, and the figure it produced has been withdrawn.
