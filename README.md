# SARS-CoV-2 and mycobacterial disease — analysis code and interactive model

Supplementary code and tools for the manuscript *SARS-CoV-2 and Mycobacterial
Disease: Population-Level Evidence for Acquired Immune Dysfunction* (Ontario,
2011–2024). University of Toronto REB #41690.

## Contents
- `index.html` — interactive **LTBI treatment-disruption model** (Supplementary
  Appendix 1). Self-contained; open locally or via GitHub Pages. Deployed at
  `https://fismanda.github.io/mycobacteria-and-respiratory-viruses/`.
- `code/` — analysis scripts:
  - `negative_controls_table2.R` — influenza/RSV negative controls, every cell of
    Table 2 (period × lag), with a built-in flu≠RSV sanity check. *(SAE)*
  - `dlnm_main.R` — main and stratified DLNM cumulative RRs, TB and NTM (Table 3). *(SAE)*
  - `attributable_fraction.R` — counterfactual SARS-CoV-2 attributable fraction
    with parametric-bootstrap CI (Figure 4). *(SAE)*
  - `forest_tb_ntm.R`, `figures_tb_ntm.R` — forest plot and figures. *(runs outside SAE)*
- `PSEUDOCODE.md` — data-assembly steps that run only inside the CIHI SAE.
- `data/README.md` — data availability (no patient-level data is shared).

*(SAE)* scripts expect an aggregate weekly analytic frame and run inside the CIHI
Secure Analytic Environment; only aggregate coefficients are exported.

## Reproducibility
R ≥ 4.2 with `dlnm`, `splines`, `MASS`. The interactive model needs only a browser.

## License
MIT (see `LICENSE`).
