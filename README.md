# MAIHDA analysis

Intersectional Multilevel Analysis of Individual Heterogeneity and
Discriminatory Accuracy (MAIHDA) for the binary outcomes in `analysis_master`,
following Evans, Leckie, Subramanian, Bell and Merlo (2024), *SSM - Population
Health*, 26, 101664.

## Layout

| Path | Purpose |
| --- | --- |
| `R/run_maihda_analysis.R` | The analysis. Runs end to end and writes every table and figure. |
| `tools/make_simulated_analysis_master.R` | Builds a simulated `analysis_master` with the same schema, for development and testing without patient data. |

## Running it

```bash
Rscript R/run_maihda_analysis.R
```

The script searches the usual locations for `analysis_master.rds` or
`simulated_analysis_master.rds`. To be explicit:

```bash
MAIHDA_DATA_PATH=/full/path/analysis_master.rds \
  MAIHDA_OUTPUT_ROOT=/full/path/outputs \
  Rscript R/run_maihda_analysis.R
```

A single-outcome validation run:

```bash
MAIHDA_TEST_MODE=1 Rscript R/run_maihda_analysis.R
```

To generate test data first:

```bash
Rscript tools/make_simulated_analysis_master.R work/simulated_analysis_master.rds 60000
```

## Intersectional strata

Strata are defined by five axes:

```
age_band x sex x ethnicity x imd_quintile x efi_category
```

Opioid strength is deliberately **not** a stratum axis. It is a treatment
characteristic rather than a social position, and folding an exposure into the
definition of the groups whose inequality is being measured would confuse the
two.

## Models

For each outcome:

- **Model A** — `outcome ~ 1 + (1 | stratum)`. Total between-stratum variation.
- **Model B** — `outcome ~ additive main effects + (1 | stratum)`. The
  remaining stratum random effects are the intersectional interaction
  residuals.

Outcomes are collapsed to one binomial record per stratum before fitting,
which retains the exact binomial likelihood.

## Minimum cell size

Strata below `MAIHDA_MIN_STRATUM_N` individuals are excluded before fitting.
The default is **10**.

The MAIHDA literature does not itself require trimming — partial pooling means
a stratum of four people contributes little to its own estimate and is shrunk
towards the additive prediction. The floor exists to satisfy statistical
disclosure control, since this analysis writes stratum-level tables and names
individual strata on its figures. Ten is the lowest threshold that does so.

Zero-event strata are **retained** by default. A stratum with 40 people and no
events is an observation of low risk, not missing information; discarding
those would bias the predicted-risk range upwards.

Every exclusion is accounted for in `stratum_retention.csv` and
`excluded_strata.csv`, and a sensitivity sweep refits the models across a
range of thresholds so the choice can be shown to be immaterial rather than
merely asserted to be.

## Singular fits

`boundary (singular) fit` on **Model B** is usually a substantive result, not a
numerical failure: the additive main effects account for essentially all
between-stratum variation, so no residual intersectional interaction remains
(VPC B ≈ 0, PCV ≈ 100%). The script suppresses lme4's bare warning and instead
confirms the boundary by refitting with alternative optimisers, records it in
the metrics, and states it explicitly on the affected figure panels.

A singular **Model A** is always flagged as a data problem, because it implies
no detectable between-stratum variation at all.

## Multiplicity

Several hundred strata are compared simultaneously within each analysis, so an
uncorrected interval excluding zero is not on its own evidence of an
interaction. Benjamini-Hochberg correction is applied across the strata within
each analysis. Panel E and Table 9 present only strata significant at
`MAIHDA_FDR_LEVEL` (default 0.05); `all_stratum_predictions.csv` retains every
stratum with both the corrected and uncorrected flags.

## Settings

All settings are environment variables.

| Variable | Default | Purpose |
| --- | --- | --- |
| `MAIHDA_DATA_PATH` | auto-discovered | Input `.rds` |
| `MAIHDA_OUTPUT_ROOT` | `maihda_outputs` | Output directory |
| `MAIHDA_TEST_MODE` | `0` | Single-outcome validation run |
| `MAIHDA_MIN_STRATUM_N` | `10` | Minimum stratum size |
| `MAIHDA_MIN_STRATUM_EVENTS` | `0` | Minimum stratum events |
| `MAIHDA_MIN_RETAINED_STRATA` | `20` | Abandon an analysis below this many strata |
| `MAIHDA_THRESHOLD_SENSITIVITY` | `1` | Run the minimum-cell sweep |
| `MAIHDA_FDR_LEVEL` | `0.05` | Benjamini-Hochberg level |
| `MAIHDA_NAGQ` | `1` | glmer `nAGQ` |
| `MAIHDA_SAVE_MODELS` | `1` | Save fitted model objects |
| `MAIHDA_RESUME` | `0` | Skip analyses whose outputs already exist |
| `MAIHDA_WRITE_SVG` | `1` | Write SVG alongside PNG |
| `MAIHDA_LABEL_SIZE_CALLOUT` | `3.4` | Stratum callout label size (mm) |
| `MAIHDA_LABEL_SIZE_AXIS` | `9.5` | Panel E y-axis label size (pt) |
| `MAIHDA_N_EXTREME` | `6` | Extreme strata retained per end |

## Outputs

Everything is written beneath `MAIHDA_OUTPUT_ROOT`, indexed by
`OUTPUT_MANIFEST.csv`.

- `tables/` — descriptives, model metrics, fixed effects, stratum predictions,
  retention and exclusion accounting, fit diagnostics
- `tables/word/` — APA-formatted Word tables, individually and as one book
- `tables/per_analysis/` — per outcome, for inspecting a single model
- `figures/individual/`, `figures/multipanel/` — PNG and SVG
- `models/` — fitted `glmer` objects
- `logs/` — run log, completion summary, `sessionInfo()`

## Requirements

Required: `lme4`, `dplyr`, `tidyr`, `ggplot2`, `patchwork`, `scales`,
`ggrepel`, `officer`, `flextable`.

Optional: `ggprism` (typography; a matched fallback theme is used when it is
absent), `svglite` (vector output; PNG is always written).

The script adapts to the installed `ggplot2` and `flextable` versions rather
than assuming one, and degrades rather than failing when an optional
dependency or a preferred typeface is missing.
