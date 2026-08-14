# MAIHDA analysis

Intersectional Multilevel Analysis of Individual Heterogeneity and
Discriminatory Accuracy (MAIHDA) for the binary outcomes in `analysis_master`,
following Evans, Leckie, Subramanian, Bell and Merlo (2024), *SSM - Population
Health*, 26, 101664.

## Layout

| Path | Purpose |
| --- | --- |
| `R/run_maihda_analysis.R` | The analysis. Runs end to end and writes every table and figure. |
| `R/generate_manuscript.R` | Reads the analysis output and writes a complete Methods and Results section with the actual values in place. |
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

The standard MAIHDA sequence, for each outcome:

- **Model A** (simple intersectional) — `outcome ~ 1 + (1 | stratum)`. Total
  between-stratum variation.
- **Model 2 family** (partially adjusted) — one model per axis, each adding a
  single axis to Model A. The PCV from each isolates that axis's contribution
  to the between-stratum variance. This is what identifies *which* dimension
  drives the variation; Model A and Model B alone cannot. Contributions
  overlap and do not sum, because the axes are correlated.
- **Model B** (intersectional interaction) — `outcome ~ additive main effects
  + (1 | stratum)`. The remaining stratum random effects are the
  intersectional interaction residuals.

Outcomes are collapsed to one binomial record per stratum before fitting,
which retains the exact binomial likelihood.

Per stratum we report the **absolute risk (AR)** and the **absolute risk due
to interaction (ARI)** — total predicted risk minus additive-only predicted
risk — using the terminology of the applied MAIHDA literature.

VPC, MOR and the variance components carry **profile-likelihood 95%
intervals**, the maximum-likelihood counterpart of the credible intervals that
Bayesian MAIHDA analyses report.

### Estimation

Much of the MAIHDA literature uses Bayesian MCMC (brms/Stan, or MLwiN via
`runmlwin`). This uses maximum likelihood with the Laplace approximation:
deterministic, no priors or convergence diagnostics, and therefore easier to
audit and reproduce inside a TRE. The minimum cell size rule also removes the
smallest strata, which is where the two approaches would most likely diverge.
The generated Methods states and justifies this explicitly.

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

**Note that this is more conservative than the field norm.** The MAIHDA
literature generally argues that precision-weighted shrinkage *is* the
multiplicity control, and is more efficient than Bonferroni-type corrections
because it does not sacrifice power. Applying FDR on top of shrinkage risks
understating the number of interactions rather than overstating it. The
generated Methods says so explicitly. Raise `MAIHDA_FDR_LEVEL`, or filter on
`interaction_distinguishable` instead of `interaction_distinguishable_fdr`, to
move closer to the conventional treatment.

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
| `MAIHDA_AXIS_DECOMPOSITION` | `1` | Fit the Model 2 family (per-axis PCV) |
| `MAIHDA_VARIANCE_INTERVALS` | `1` | Profile-likelihood CIs for VPC/MOR |
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

## Generating the write-up

After the analysis has run:

```bash
MAIHDA_OUTPUT_ROOT=/path/to/outputs Rscript R/generate_manuscript.R
```

This reads the analysis CSVs and writes `Manuscript.docx`, `.md` and `.txt`
containing a complete Methods and Results section with the actual values
substituted in. Nothing is recomputed — every figure is read from the analysis
tables, so the text and the results tables cannot drift apart.

It adapts to the run it is given: the axes, categories, reference levels,
thresholds, outcomes and software versions all come from the output, so
changing the analysis needs no change here. It also handles the cases that
matter — a Model B variance at zero is described as a substantive finding, and
a run with no interactions surviving FDR correction produces a properly worded
null result rather than an empty section.

### Study metadata

Some Methods content cannot come from the analysis: data source, setting,
inclusion criteria, ethics approval. Set these in section 1 of the script or
via `MANUSCRIPT_*` environment variables. Anything left unset appears in the
text as a conspicuous `[PLACEHOLDER: ...]` marker **and** is listed in
`placeholders_to_complete.txt`, so nothing can be missed on a read-through.

### Disclosure control

The generated text names individual strata and quotes counts within them,
which is what a TRE output checker will scrutinise. `MANUSCRIPT_DISCLOSURE_MODE`
controls this:

| Mode | Behaviour |
| --- | --- |
| `flag` (default) | Report true values, and list every figure derived from a cell below the threshold in `disclosure_check.csv`. Use while drafting. |
| `apply` | Suppress or round small counts in the text itself. Use for the version that leaves the environment. |
| `off` | No disclosure processing. Only appropriate outside a TRE. |

`MANUSCRIPT_DISCLOSURE_THRESHOLD` defaults to the analysis minimum cell size.
Set `MANUSCRIPT_DISCLOSURE_ROUNDING=5` for the common round-to-nearest-5 rule;
in that case *all* counts are rounded, not only small ones, since a mix of
exact and rounded figures lets an exact value be recovered by difference.

Every run also writes `manuscript_values.csv` — each headline figure quoted in
the text with its source table and column — so any number can be traced back
without re-reading the prose, and an output checker has a single list to work
from.

## Requirements

Required: `lme4`, `dplyr`, `tidyr`, `ggplot2`, `patchwork`, `scales`,
`ggrepel`, `officer`, `flextable`.

Optional: `ggprism` (typography; a matched fallback theme is used when it is
absent), `svglite` (vector output; PNG is always written).

The script adapts to the installed `ggplot2` and `flextable` versions rather
than assuming one, and degrades rather than failing when an optional
dependency or a preferred typeface is missing.
