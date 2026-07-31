# Code audit

## Scope reviewed

The final package was reviewed across all nine R source files for:

- cross-file function and parameter compatibility;
- fixed and random heterogeneity posterior calculations;
- posterior-predictive simulation and evidence updating;
- reuse of simulations across equity weights and decisions;
- deterministic reuse across recruitment-cost scenarios;
- equity-weight normalization;
- common and subgroup-specific decision rules;
- trial-participant consequences;
- participation burdens of 0.010 QALYs for Group A and 0.005 QALYs for Group B;
- staircase recruitment-cost arithmetic;
- wENBS component identity;
- exact grid maximization;
- main and appendix row counts;
- 9-panel figure construction;
- stale-output protection and output signatures;
- logging and repeated-run comparisons.

## Changes made in this audited version

- Rebuilt the workflow using base R only; no `dplyr`, `readr`, or `ggplot2` dependency remains.
- Removed obsolete near-optimal and distinguishability-rule code.
- Retained only the final recruitment-cost scenarios: $0, $3,000, and $6,000.
- Added the base-case participation burden for both groups to every main and appendix wENBS calculation.
- Ensured that each allocation's evidence is simulated before the loops over equity weights and decision structures.
- Applied recruitment costs after wEVSI simulation, so cost scenarios cannot receive different random draws.
- Added simulation fingerprints and explicit reuse-validation files.
- Added analysis signatures to result files and stale-output checks to figure scripts.
- Added clean child R sessions and complete stdout/stderr logs.
- Added function-level and full-output validation tests.
- Added same-code/same-runtime repeated-run comparison of the six key CSV outputs.

## Static checks completed

- Balanced parentheses, brackets, braces, quoted strings, and comments in every R file.
- No native `|>` pipe or contributed-package pipe remains.
- No `library()` or `require()` calls remain.
- No obsolete cost levels ($2,000, $5,000, or $10,000) remain as analysis scenarios.
- No obsolete near-optimal threshold variables remain.
- Every file listed in the analysis signature is included in the package.
- The archive and individual source-file checksums were generated after the audit.

## Numerical validation

The same equations were previously executed independently in Python with 20,000 simulations for the complete main and appendix analyses. That implementation reproduced all component identities, simulation-reuse requirements, cost monotonicity, burden calculations, and same-seed repeatability.

## Changes made after the initial audit

- Fixed a fatal path bug: `Rscript` encodes spaces in the `--file=` argument as
  `~+~`, so every child script resolved its root directory incorrectly and the
  workflow could not run from a path containing spaces. All six scripts that
  parse `--file=` now decode it.
- Set the base-case participation burden to 0.010 QALYs for both subgroups. The
  earlier asymmetric burden (0.005 for Group B) assumed a per-person difference
  that the case study has no evidence for.
- Restricted the appendix evidence-profile analyses to the $0 recruitment-cost
  increment; the recruitment-cost gradient is explored in the main analysis.
- Rebuilt all figures: larger text, the equity-weight legend moved to a reserved
  strip below the panels, rotated row strips, and a padded allocation axis so
  maxima on the first or last grid point are not clipped by the panel border.
- Combined the two decision structures into a single appendix figure per
  evidence profile.

## Execution status

The complete workflow has been executed under R 4.4.1 on macOS (arm64). All 56
validation tests pass, and two consecutive runs produced byte-identical output
for all six key CSV files.

## Limitation

The Monte Carlo standard error of wENBS is approximately 20-50 QALYs per
allocation. Under near-complete borrowing the wENBS curve varies by less than
this across the entire 5%-95% allocation range, so the maximizing allocation is
not identified in that column and the plotted maxima there should not be read as
preferred designs. The curve, not the point, should be used to judge how sharply
an allocation is identified.
