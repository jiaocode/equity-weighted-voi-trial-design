# Final audited R workflow: inclusive clinical-trial design

This package contains the complete main and appendix analyses for the two-group hierarchical wEVSI/wENBS case study.

## Final specifications

- Total trial sample size: 600
- Group A population share: 17%; Group B: 83%
- Equity weights for Group A: 1, 1.5, and 2; Group B weight: 1
- Recruitment-cost increments for Group A: **$0, $3,000, and $6,000** in the
  main analysis; the appendix evidence profiles use **$0 only**
- Trial-participation burden: **0.010 QALYs lost per participant in both
  groups** (the base case does not assume a differential burden)
- Simulation size: 20,000 draws and 50 Monte Carlo batches
- Allocation grid: 5%–95% Group A in one-percentage-point increments
- Borrowing models: near-complete, moderate, and near-no borrowing
- Main decisions: common population-level and subgroup-specific
- Appendix evidence profiles: reversed, near-threshold, and away-from-threshold

Participation burdens are stored in `00_analysis_parameters.R` as positive loss magnitudes and are subtracted from wENBS. The same normalized equity weights are applied to health benefits and participation burdens.

## Figures

- Figure 1: wEVSI under a common population-level decision
- Figure 2: wEVSI under subgroup-specific decisions
- Figure 3: 9-panel wENBS under a common decision
- Figure 4: 9-panel wENBS under subgroup-specific decisions

Figures 3 and 4 use:

- rows: $0, $3,000, and $6,000 recruitment-cost increments;
- columns: near-complete, moderate, and near-no borrowing;
- lines: equity weights 1, 1.5, and 2;
- filled points: exact numerical maxima of the displayed curves;
- dashed vertical line: Group A's 17% population share;
- dotted horizontal line: wENBS = 0.

The appendix produces a three-panel wENBS figure (one column per borrowing
model) for each alternative evidence profile and both decision structures, at
the $0 recruitment-cost increment only. There is no arbitrary near-optimal or
statistical-distinguishability rule.

All figures place the equity-weight legend in a reserved strip below the
panels, so it cannot overlap a curve. Text sizing is controlled centrally by
the `PLOT` settings in `00_analysis_parameters.R`; every figure script resets
`par(cex = 1)` after `mfrow`, because base R otherwise shrinks all text to 0.66
of nominal size in layouts with three or more rows or columns.

## Reproducibility design

For each evidence profile, borrowing model, and allocation:

1. The future evidence is simulated once.
2. That evidence is reused across all equity weights.
3. It is also reused across both decision structures.
4. Recruitment costs are applied deterministically without resimulation.

The code fixes the R random-number algorithms, seeds, simulation size, batching, allocation grid, and numerical threading. Output files contain an analysis signature based on all R source files. Figure scripts refuse to use results created by a different code signature.

## Requirements

- R 4.0.0 or later
- No contributed R packages are required

## Run the complete analysis

From a terminal in the extracted folder:

```bash
./run_final_analysis.sh
```

Or from any operating system with Rscript available:

```bash
Rscript --vanilla 07_run_all.R
```

The driver deletes old `results/` and `figures/` folders before running. Each substantive script runs in a separate clean R session.

## Error logs

If a child script fails, the driver prints its complete standard output and error logs. The files are also retained under:

```text
results/logs/
```

For example:

```text
results/logs/03_run_main_analysis_stdout.log
results/logs/03_run_main_analysis_stderr.log
```

## Repeated-run check

When an earlier run has the same code signature and runtime signature, the workflow preserves its six key CSV files before cleaning and compares them with the new run. Results are written to:

```text
results/reproducibility_check.csv
```

A directly comparable run stops with an error if any key CSV differs byte-for-byte.

## File guide

- `00_analysis_parameters.R`: all fixed inputs, seeds, costs, burdens, scenarios, and plotting settings
- `01_hierarchical_model_functions.R`: hierarchical posterior and predictive calculations
- `02_value_functions.R`: equity weights, wEVSI, wENBS, costs, burdens, and exact maximization
- `03_run_main_analysis.R`: complete main analysis and internal validations
- `04_make_main_figures.R`: Figures 1–4
- `05_run_appendix_analysis.R`: complete appendix analyses and internal validations
- `06_make_appendix_figures.R`: six appendix 9-panel figures
- `07_run_all.R`: clean workflow driver, logging, signatures, and repeated-run comparison
- `08_run_validation_tests.R`: deterministic function tests and validation of all outputs and figures
- `run_final_analysis.sh`: terminal wrapper
- `CODE_AUDIT.md`: audit scope and limitations

## Important interpretation note

The filled point is the exact grid maximum. On a nearly flat wENBS curve, small Monte Carlo differences can move the point even when the curve and substantive conclusion are unchanged. The curve—not only the point—should be used to judge how sharply the preferred allocation is identified.
