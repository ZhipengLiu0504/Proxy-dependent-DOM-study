# Proxy-dependent DOM study

This repository contains R scripts and supporting data to reproduce the analyses behind main Figures 1-5
and supplementary Figures S3-S7. It does not modify the original scripts,
manuscript, supplementary information, or source data.

## Run

Use R 4.5.2 or a compatible installation. The tested mixOmics version is 6.34.0.
Required packages: readr, readxl, dplyr, tidyr, purrr, stringr, tibble, ggplot2,
scales, patchwork, and mixOmics. The optional preparation script additionally
requires writexl. Dependencies are checked, never installed automatically.
Detailed tested versions are recorded in each output folder's `sessionInfo.txt`.

Set the working directory to this release folder, then run:

```r
source("run_all.R")
source("tests/validate_outputs.R")
```

Or from a terminal in this directory:

```text
Rscript run_all.R --check-only
Rscript run_all.R
Rscript run_all.R --stage=03
Rscript run_all.R --stage=01,02,05
Rscript tests/validate_outputs.R
```

Only this release's generated outputs are replaced when rerunning. Keep any
hand-edited figures elsewhere. Do not run individual module files directly:
the driver provides their input paths and isolated output directories.

## Structure and figure mapping

| Module | Manuscript figures |
|---|---|
| `scripts/01_ms_composition.R` | 1a-d, S3, S4 |
| `scripts/02_community_pathways.R` | 2a, 2c, S5 |
| `scripts/03_proxy_diablo.R` | 3a-c, 4a-d, 5a-b, S7 |
| `scripts/04_formula_diablo.R` | S6 |
| `scripts/05_taxa_function.R` | 2b, 2d |

`figure_map.csv` maps all 22 panels to their plotting objects and export sizes.
Running the workflow creates matched PNG (300 dpi) and vector PDF exports in `figures/`. These are
individual reproducibility panels, not newly assembled journal submission plates.
Recheck typography and dimensions when assembling the final plates.

Running the workflow creates `outputs/` with numeric CSV results, complete intermediate tables in
`analysis_tables.rds`, fitted DIABLO models, diagnostic PDFs, and session details.
For example, use `readRDS("outputs/03_proxy_diablo/analysis_tables.rds")` to inspect
the matrices and data frames underlying Figures 3-5 and S7.

`data/` contains unchanged copies of the six main inputs. The two pathway files
ending in `.xls` are tab-delimited text and are intentionally read as text.
`data/vendor_assignments/` contains the ten vendor formula-assignment CSVs.
These are NOT instrument-raw FT-ICR-MS spectra.

The original SI2, SI3 and SI4 reference-table files are not distributed.
The shared inputs and analysis scripts can regenerate their underlying results.
Public validation checks test internal consistency without requiring these
reference files. SI dataset 5 is a separate pathway statistics table; this
repository does not reconstruct its external statistical workflow.

## Optional upstream MS preparation

```r
source("optional/prepare_ms_assignments.R")
```

This retains the original H + 1 convention, mass/elemental-ratio filters,
AI-mod calculation, and per-mass formula selection. It writes a separate
`outputs/00_optional_ms_preparation/regenerated_ms_assignments.xlsx`.
It never replaces `data/ms_formula_assignments.xlsx`, which remains the verified
manuscript input. The ten regenerated sheets were compared with the manuscript
workbook; shared columns and row counts agreed in the local test.

## Figure contract and scientific interpretation

The purpose is faithful reproduction, not a new analysis or a redesign.
The quantitative panels describe habitat differences in DOM composition,
microbial communities and predicted functions, followed by multiblock
associations and within-reservoir comparisons. Preserve the original data,
panel membership, feature selection, pooling hierarchy and statistical tests.
R is the only plotting/export backend used in this release.

- Fig. 1/S3/S4: mutually exclusive sequential molecular-pool rules; CRAM is
  classified first and condensed aromatics then use AI-mod >= 0.67. RDOC and
  BDOC are compositional proxies, not direct persistence/lability measurements.
- Fig. 2a/2c/S5 use supplied habitat-level summaries. Fig. 2b/2d use the ten
  samples matched to the MS metadata. The full metadata table has 24 entries;
  authors should verify how upstream habitat summaries were constructed.
- S5 retains `Others` in normalization; Fig. 2a excludes it before normalization.
  These original denominators are deliberately not harmonized silently.
- Fig. 2c reports descriptive pathway log2 ratios with a 1e-6 pseudocount on
  percentage abundance. "Enriched" in its original labels does not imply a
  significant test. PICRUSt2 predicts functional potential, not measured activity.
- Fig. 3 uses proxy-level DIABLO; S6 uses formula-level DIABLO. Their feature
  prefilters and genus transformations differ and are preserved. Their scores
  are in-sample and are not independent validation of predictive accuracy.
- Fig. 3b/c show component-1 **loadings**, not fold changes or p values. Signs
  are model-orientation dependent. Validation checks the retained sign-to-habitat
  coloring in the tested fit.
- Fig. 4a/b: five samples per habitat, mean +/- SD, two-sided unpaired Welch tests.
  Fig. 4c/d: three matched pairs in reservoir 1, mean +/- SD, two-sided paired
  t tests. Lines connect paired observations. Stars show unadjusted p values.
- Fig. 5: ten matched samples, Spearman tests with `exact = FALSE`; retain
  |rho| > 0.6 and raw p < 0.05. BH-adjusted p values remain in exported tables
  but are not the plotted inclusion criterion. These are exploratory associations.
- S7: summed assigned peak signal, not absolute DOC concentration; five samples
  per habitat and a two-sided Welch test. The boxplot uses ggplot2 defaults.

The author should decide separately whether to change inferential methods,
multiple-testing policy, manuscript wording, or the use of habitat summaries.
None of those scientific choices is silently changed by this code cleanup.

## Changes and verification

- Retained the analysis blocks needed by the current main and supplementary figures.
- Excluded obsolete/debug alternatives, unused maps, extra absolute-signal plots,
  unused formula-level networks, and exploratory plots not in the selected figures.
- Replaced Chinese tutorial/debug comments with concise English explanations;
  translated executable messages and removed machine-specific input paths.
- Added isolated stages, dependency/input checks, deterministic plotting seed,
  figure-numbered exports, saved models/tables, and numeric regression checks.
- Fixed the habitat-summary Cartesian-completion bug: the old block added zero
  observations for samples in the wrong habitat. Normalized stacked proportions
  were unchanged, but unnormalized means/SD/SE were misleading. The corrected
  block reuses only the ten observed sample-habitat pairs.
- Removed an undefined interactive-session-only plotting object at the old MS
  script's end. Removed duplicate/unused plots, not underlying necessary analyses.
- Adjusted title line breaks and long-label spacing; used locale-safe English
  source and PDF fonts. Statistical values and selection thresholds are unchanged.

After running the workflow and public tests, see `outputs/validation_checks.csv`.
These validate internal consistency; they do not certify every assumption of
the study design or establish compliance with all journal submission requirements.

## Citation and licensing

Repository: https://github.com/ZhipengLiu0504/Proxy-dependent-DOM-study

Article citation and archive DOI will be added when available. No code or data
license has been specified yet. No third-party package source is bundled.
