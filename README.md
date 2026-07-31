# Equity-weighted value of information for inclusive clinical-trial design

R code reproducing the main and appendix analyses for the two-group
hierarchical wEVSI/wENBS case study: how many trial participants to enrol from
each population subgroup when the decision-maker places extra weight on health
gains in a historically marginalized group.

## Requirements

- R 4.0.0 or later
- No contributed R packages

## Run the analysis

From a terminal in this folder:

```bash
./run_analysis.sh
```

Or, on any system with `Rscript` available:

```bash
Rscript --vanilla 07_run_all.R
```

The driver deletes any existing `results/` and `figures/` folders, then runs
each substantive script in a separate `--vanilla` R session. A complete run
takes roughly three minutes on a current laptop.

## Case-study inputs

| Input | Value |
| --- | --- |
| Population share | Group A 17%, Group B 83% |
| Current estimated INB (standard error) | Group A 0.41 (0.57), Group B −0.15 (0.18) |
| Total trial sample size | 600, randomized 1:1 within each subgroup |
| Outcome standard deviation | 2 in each subgroup |
| Allocation grid | 5%–95% Group A, in one-percentage-point steps |
| Borrowing models | near-complete (ω = 0.001), moderate (ω ~ Half-Normal(0, 0.30)), near-no (ω = 3.0) |
| Decision structures | common population-level, subgroup-specific |
| Raw equity weights | Group A 1, 1.5, or 2; Group B 1 |
| Annual eligible population | 15,000 over a 10-cohort horizon, 3% discount rate |
| Time to reporting | 3 years, for every allocation |
| Fixed trial cost | $19 million |
| Base variable cost | $53,000 per participant |
| Group A recruitment-cost increment | $0, $3,000, or $6,000 per participant, applied in successive blocks of 50 after the first 50 Group A participants |
| Cost-effectiveness threshold | $150,000 per QALY |
| Participation burden | 0.010 QALYs per participant in both subgroups |
| Simulation size | 20,000 draws in 50 batches |

The appendix evidence profiles (reversed, near-threshold, and
away-from-threshold) use the $0 recruitment-cost increment only; the
recruitment-cost gradient is explored in the main analysis.

Participation burdens are stored as positive loss magnitudes and subtracted
from wENBS. The same normalized equity weights are applied to health benefits
and to participation burdens.

## Figures

- Figure 1: per-person wEVSI under a common population-level decision
- Figure 2: per-person wEVSI under subgroup-specific decisions
- Figure 3: wENBS under a common population-level decision
- Figure 4: wENBS under subgroup-specific decisions

Figures 3 and 4 have recruitment-cost rows ($0, $3,000, $6,000) and
borrowing-model columns. The appendix produces one figure per evidence profile,
with a common decision in the upper row and subgroup-specific decisions in the
lower row. In every figure, lines are equity weights, filled points are the
grid maxima of the plotted curves, the dashed vertical line is Group A's 17%
population share, and the dotted horizontal line is zero.

## Reproducibility

For each evidence profile, borrowing model, and allocation, the future evidence
is simulated once and reused across all equity weights, both decision
structures, and every recruitment-cost level, which are applied
deterministically without resimulation. The code fixes the random-number
algorithms, seeds, simulation size, batching, allocation grid, and numerical
threading.

Output files carry a signature derived from the contents of all R source files.
The figure scripts refuse to plot results produced by a different code version.
When a previous run has the same code and runtime signature, the workflow
compares its six key CSV files with the new run and stops if any differ
byte-for-byte; the comparison is written to `results/reproducibility_check.csv`.

`08_run_validation_tests.R` runs as part of the workflow and checks the wENBS
component identity, equity-weight normalization, the participant-burden and
staircase-cost arithmetic, reuse of simulation draws, row counts, and agreement
between the reported optima and the plotted curves. The run aborts if any check
fails, and the results are written to `results/validation_summary.csv`.

## Interpreting the maxima

The filled point in each panel is the grid maximum. The Monte Carlo standard
error of wENBS is roughly 20–50 QALYs per allocation. Under near-complete
borrowing the wENBS curve varies by less than this across the entire allocation
range, so the maximizing allocation is not identified in that column and the
plotted point should not be read as a preferred design. The shape of the curve,
not the location of the point, indicates how sharply an allocation is
identified.

## Error logs

If a script fails, the driver prints its standard output and error logs and
retains them under `results/logs/`, for example
`results/logs/03_run_main_analysis_stderr.log`.

## Files

- `00_analysis_parameters.R`: inputs, seeds, costs, burdens, scenarios, and plotting settings
- `01_hierarchical_model_functions.R`: hierarchical posterior and predictive calculations
- `02_value_functions.R`: equity weights, wEVSI, wENBS, costs, burdens, and grid maximization
- `03_run_main_analysis.R`: main analysis
- `04_make_main_figures.R`: Figures 1–4
- `05_run_appendix_analysis.R`: appendix evidence-profile analyses
- `06_make_appendix_figures.R`: appendix figures
- `07_run_all.R`: workflow driver, logging, signatures, and repeated-run comparison
- `08_run_validation_tests.R`: validation of functions, outputs, and figures
- `run_analysis.sh`: terminal wrapper

## Licence

MIT. See `LICENSE`.
