# Quinone reduction-energy representation dependence — source release candidate (V7)

**STATUS: V7 AUTHOR-APPROVED RELEASE CANDIDATE — not yet public; license and DOI pending. Numerical source tables and M6 raw records are unchanged from V6.**

## Scope and provenance
This artifact contains (a) source-linked tabular results from the original M1-M6 freeze, (b) all 24 converged M6 raw input/output/geometry triplets, (c) plot artwork and scripts, and (d) a pure-Python numerical audit. Its creation has NOT changed the historical calculations or corrected the M1 raw-output exceptions. It is NOT an experimental potential dataset.

## Verify without ORCA or pip-installed packages
From this directory run:

```bash
python3 scripts/verify_v5_source_data.py  # preserved numerical verifier from V5; results unchanged
```

The verifier checks all M1 entity and parent-relative shift arithmetic, M2R historical and corrected-vintage median statistics, M5 component cancellation, all M6 energy-derived shifts, and presence, normal termination, and XYZ SHA-256 matches for all 24 M6 raw jobs. Passing this script does **not** prove experimental validity, reproduce M1-M5 upstream DFT or validate unreconciled M1 source outputs.

## Optional regenerated figure proofs

```bash
python3 scripts/render_figures.py
```

Requires `rdkit`, `numpy`, `matplotlib` and `cairosvg`. The script reads CSVs in `data/` and writes proofs to `figures/` (may overwrite supplied local copies). Molecular structures use source SMILES; figure captions/interpretations must be read in the manuscript.

## Important method and validity boundaries
- Nominal neutral singlet / radical-anion doublet electronic descriptor: ΔEred = (E[-1]−E[0]) × 27.211386245988 eV/Eh.
- A positive parent-relative shift means less favorable *calculated electronic attachment* by the defined sign convention; not a voltammetric potential.
- Historical M2R ten-site median and M3R four corrected cases are different data vintages; a mixed-vintage 10-site sensitivity is illustrative only.
- M5 geometry and side-chain contribution labels are defined by a chosen counterfactual and algebraic subtraction; not directly isolated chemical mechanisms.
- M6 uses source-fixed geometries; the two targeted sign reversals are functional-robust within that limited panel, not an independent experimental test.
- Four M1 production single-point outputs are not reconciled: see `data/M1_FOUR_RAW_OUTPUT_EXCEPTIONS.csv`.
- The menadione/thiodione historical comparator does not represent a matched experimental one-electron benchmark; no experimental error metric is claimed.
- The 0.05 eV rule belongs to the frozen decision protocol and has not been retuned.

## Upstream archived source data
The historical M1-M5 archive remains separately preserved at `QUINONE_M1_M5_LEAN_RAW.tar.gz` (not included in this small repository candidate because it is approximately 131 MB). Its SHA-256 is **6ed5778dae6b7e7edc79461836bb0fc598c4208208c5729a735b4d8c328175ee**; obtain the underlying archive from the author prior to claiming upstream quantum reproducibility. Four M1 raw gaps remain even in that archive.

## Publication requirements still open
Author name, affiliation, computing-resource statement, funding, competing interests and AI-assisted workflow have been confirmed in V7. Public release has been authorized in principle, but a new, dedicated repository does not yet exist, reuse-rights/license choice remains outstanding, and no DOI has been assigned. No public release or journal submission has occurred. Do not invent a public DOI.

## Author and citation metadata
**Author:** Aris Tsai (ORCID 0009-0004-8641-7647). **Affiliation:** Stanford Online High School. **Correspondence:** aristsai09@gmail.com. See `CITATION.cff` and `AI_ASSISTANCE_AND_FIGURE_PROVENANCE.md`. An institutional email supplied to manuscript editors is not a coauthor address and is not reproduced here.
