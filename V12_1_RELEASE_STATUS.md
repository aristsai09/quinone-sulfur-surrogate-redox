# V12.1 submission data update

This additive update preserves all historical M1–M6 files. It provides the retrospective 34-site and eight-site source-corroborated tables, excluded-record list, descriptive sensitivity tables, and a deterministic figure-rebuilding script.

Run `python3 scripts/render_v12_1_corrob_figures.py` from the repository root (requires matplotlib). Figures are written under `figure_assets/` as publication-resolution PNG/PDF files. The manuscript's original V12.1 figure files are separately provided in the author-reviewed submission bundle; historical plots remain in the original `figures/` directory and the manuscript SI.

The frozen historical ten-site median of 0.175899559 eV and the post hoc source-corroborated eight-site median of 0.162916518 eV are **distinct**. The separate four-case M3R corrections must not be silently merged into original M2R values. The mixed-vintage eight-site summary (0.111650561 eV) is **post hoc** and not a homogeneous confirmatory result.

The V12.1 machine-readable Supporting Information ZIP also includes forensic failed-run records and audits; those large ZIPs are not part of this small GitHub patch. Four historical M1 final-output gaps remain unresolved. M6 is a 24-job fixed-geometry functional-sensitivity analysis, not an experimental electrochemical validation. Further ORCA jobs were not accepted as replacement energy sources.

The author supplied Zenodo DOI https://doi.org/10.5281/zenodo.23459391 for the original research archive. The precise Zenodo deposited file inventory has not been independently verified in this update. Do not imply that the Zenodo version already contains the V12.1 files unless independently confirmed.
