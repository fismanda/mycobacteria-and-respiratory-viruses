# SAE-only data assembly (pseudocode)

Row-level data live inside the CIHI Secure Analytic Environment (SAE) and cannot
leave it, so the steps that build the weekly analytic frame `dat` are given here
as pseudocode. Only aggregate model coefficients are exported. Ethics: University
of Toronto REB #41690.

## 1. Case ascertainment (Stata/R inside SAE)
- Source: linked CIHI Discharge Abstract Database (DAD, inpatient) + National
  Ambulatory Care Reporting System (NACRS, ED).
- TB = ICD-10-CA A15, A17, A18, A19 in any diagnosis field (exclude A16).
  NTM = A31. Keep earliest qualifying record per person per outcome (incident).
- Immunocompromising host condition (ICH): HIV/AIDS (B20-B24, Z21), transplant
  (Z940-Z944, Z941, Z948), haematologic malignancy (C81-C85,C88,C90-C96),
  primary immunodeficiency (D80-D84).
- Restrict to the 14 GTHA + central-Ontario public health units.
- Aggregate to weekly counts by (series_week, age_cat, ich): n_tb, n_ntm.

## 2. Population offset
- Weekly population by linear interpolation of annual PHU denominators.

## 3. Respiratory viral exposures
- SARS-CoV-2: test-adjusted case counts from CCM (Mar 2020-Aug 2022; correction
  for differential testing by age/sex, Fisman 2021 / Bosco 2025), bridged to
  RVDSS percent positivity (Sep 2022-Mar 2024).
- Influenza A and RSV: RVDSS for the full period.
- Standardize each series to SD units over the study period (natural zero kept).

## 4. Seasonality and trend
- Fourier terms sin1,cos1,sin2,cos2 at period 52.18 weeks; linear secular trend
  series_week.

## 5. Modelling
- See code/: dlnm_main.R (Table 3), negative_controls_table2.R (Table 2),
  attributable_fraction.R (Figure 4 / PAF). Forest and figure scripts run
  outside the SAE on the exported aggregate estimates.
