# Respiratory viruses and mycobacterial disease — analysis code and interactive scenario tool

Supplementary code and tools for the manuscript *Pandemic-Era Increase in Tuberculosis
Acute-Care Presentations Was Not Associated With SARS-CoV-2 Activity: A Population-Based
Time-Series Study in Ontario, Canada* (Greater Toronto–Hamilton Area Plus, 2011–2024).
University of Toronto REB #41690.

> **Note on earlier releases.** Releases before v2.0.0 accompany a superseded version of
> this analysis. A length mismatch between the mycobacterial outcome series and the viral
> exposure series inflated estimated SARS-CoV-2 effects. All estimates in the current
> release derive from the corrected, week-for-week aligned series (676 analysed weeks).
> Earlier releases are retained for provenance and should not be cited as current.

## What the study found

Acute-care presentations for tuberculosis ran modestly above their pre-pandemic
trajectory, by an amount that depends on how that trajectory is specified — rate ratios
from 1.031 to 1.217 across six specifications, with both BIC-selected models giving 1.095.
NTM presentations did not increase. Weekly SARS-CoV-2 activity was **not** followed by
increased presentations for either outcome, at lags through 12 weeks (TB 0.909,
95% CI 0.822–1.006; NTM 0.953, 0.868–1.047) or 26 weeks. RSV activity was associated with
both outcomes, establishing that the design detects viral-activity signals where they
exist.

## Contents

- `index.html` — interactive **scenario tool** for two mechanisms that require no immune
  dysfunction: disruption of latent-TB preventive treatment, and growth in the size and
  source-country composition of the arriving population. The reader selects which of the
  six counterfactual specifications to treat as the excess to be explained, and sets every
  parameter. Self-contained; open locally or via GitHub Pages at
  `https://fismanda.github.io/mycobacteria-and-respiratory-viruses/`.

  **This tool evaluates order-of-magnitude plausibility under user-selected assumptions.
  It does not estimate causal attribution.** Every quantity entering the calculation is a
  control in the interface; the only fixed constants are catchment population (11,200,000),
  pre-pandemic TB rate (2.96 per 100,000/year) and window length (4 years). Both formulas
  are printed in the interface with the user's values substituted.

- `code/` — analysis scripts:
  - `dlnm_main.R` — distributed lag non-linear models, cumulative rate ratios per 1-SD of
    viral activity, overall and by age and immunocompromise stratum, for all eight
    exposures; extended-lag sensitivity analyses at 26 and 52 weeks. *(SAE)*
  - `comparator_exposures.R` — influenza, RSV and the five 2015-onward exposures, each
    adjusted for concurrent SARS-CoV-2 activity. *(SAE)*
  - `level_shift.R` — pandemic-era level shift under two estimation families
    (pre-pandemic projection, and full-series indicator) crossed with three trend
    specifications, with AIC and BIC. *(SAE)*
  - `cumulative_excess.R` — quarterly accumulation of observed-minus-expected
    presentations, and the excess implied by a constant rate ratio. *(SAE)*
  - `forest_tb_ntm.R`, `figures_tb_ntm.R` — forest plot and figures. *(runs outside SAE)*

- `PSEUDOCODE.md` — data-assembly steps that run only inside the CIHI SAE, and the model
  specifications in full.
- `data/README.md` — data availability. No patient-level data is shared.

*(SAE)* scripts expect an aggregate weekly analytic frame and run inside the CIHI Secure
Analytic Environment; only aggregate coefficients are exported.

### Withdrawn

`attributable_fraction.R` computed a SARS-CoV-2 attributable fraction. An attributable
fraction is not defined for a null exposure association, and both the script and the figure
it produced have been withdrawn. The former "negative controls" framing has also been
dropped: influenza and RSV were pre-specified comparator exposures, and RSV proved to be
positively associated with both outcomes, so neither functions as a negative control.

## Reproducibility

R >= 4.2 with `dlnm`, `splines`, `MASS`. The interactive tool needs only a browser.

## Citation

See `CITATION.cff`. Please cite the release matching the version of the analysis you are
referring to.

## License

MIT (see `LICENSE`).
