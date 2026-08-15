# =============================================================================
# MAIHDA analysis of analysis_master
# =============================================================================
#
# PURPOSE
# -------
# This script implements an Intersectional Multilevel Analysis of Individual
# Heterogeneity and Discriminatory Accuracy (MAIHDA) for the binary outcomes in
# analysis_master. It follows the model sequence and presentation logic in:
#
# Evans CR, Leckie G, Subramanian SV, Bell A, Merlo J. (2024).
# A tutorial for conducting intersectional multilevel analysis of individual
# heterogeneity and discriminatory accuracy (MAIHDA).
# SSM - Population Health, 26, 101664.
# https://doi.org/10.1016/j.ssmph.2024.101664
#
# INTERSECTIONAL STRATA
# ---------------------
# Strata are defined by five axes:
#
#   age_band x sex x ethnicity x imd_quintile x efi_category
#
# Opioid strength is deliberately excluded. It is a treatment characteristic
# rather than a social position, and folding an exposure into the definition of
# the groups whose inequality is being measured would confuse the two. If it is
# ever wanted as an adjustment covariate rather than a stratum axis, that is a
# different model and should be added explicitly.
#
# For each specification and outcome, two logistic MAIHDA models are fitted:
#
#   Model A: outcome ~ 1 + (1 | stratum)
#            This null model quantifies the total between-stratum variation.
#
#   Model B: outcome ~ additive main effects + (1 | stratum)
#            The remaining stratum random effects represent departures from
#            the additive prediction. Within the MAIHDA framework these are
#            interpreted as intersectional interaction residuals.
#
# Binary outcomes are collapsed to one binomial record per stratum before
# fitting. This is the computationally efficient approach described in section
# 2.5.C of the tutorial and retains the exact binomial likelihood.
#
# DESCRIPTIVE STATISTICS
# ----------------------
# Section 3B produces the descriptive material that should accompany any
# MAIHDA result: cohort characteristics by axis, crude event rates by category
# with Wilson intervals, a complete stratum-level listing including the strata
# that will later be trimmed, and the stratum size distribution. These are
# computed on the unfiltered cohort so that the analytic sample can always be
# reconciled against the cohort it came from.
#
# MINIMUM CELL SIZE
# -----------------
# A five-way stratification produces a long tail of very small cells. Strata
# below MIN_STRATUM_N individuals, or below MIN_STRATUM_EVENTS events, are
# excluded before fitting. The default threshold is n >= 10, which is the
# conventional statistical disclosure control floor for published cell-level
# figures; see the extended note in section 1 for why the MAIHDA literature
# itself does not require trimming. Every exclusion is accounted for: the
# script reports how many strata, individuals and events were removed, what
# share of the cohort was retained, and how the headline variance measures
# respond to the threshold chosen. See section 5 and the stratum_retention
# tables.
#
# The rule is applied AFTER the binomial collapse and BEFORE model fitting, so
# the reported model quantities all refer to the retained analytic sample. The
# unfiltered stratum size distribution is still reported, so the effect of the
# rule is always visible.
#
# SINGULAR FITS
# -------------
# lme4 reports "boundary (singular) fit" when a random-intercept variance is
# estimated at zero. For Model B this is frequently the correct substantive
# answer rather than a numerical failure: it means the additive main effects
# account for essentially all between-stratum variation and no residual
# intersectional interaction remains. The script therefore suppresses lme4's
# bare warning and instead:
#
#   * confirms the boundary by refitting with alternative optimisers, so a
#     genuine zero variance is distinguished from an optimiser that stalled;
#   * records the outcome in the model metrics and run log;
#   * states the result explicitly on the affected figure panels.
#
# A singular Model A, by contrast, is always flagged as a data problem, because
# it implies no detectable between-stratum variation at all.
#
# OUTPUTS
# -------
# The script writes model objects, CSV tables, APA-formatted Word tables,
# individual figures, combined multi-panel figures and cross-outcome summary
# figures beneath OUTPUT_ROOT. An OUTPUT_MANIFEST.csv lists every file written.
# Nothing in analysis_master is overwritten.
#
# REQUIRED PACKAGES
# -----------------
# Required: lme4, dplyr, tidyr, ggplot2, patchwork, scales, ggrepel, officer,
# flextable.
# Optional: ggprism (typography; a matched fallback theme is used when it is
# unavailable), svglite (vector output; PNG is always written).
#
# RUNNING THE SCRIPT
# ------------------
# The simplest use is:
#
#   Rscript run_maihda_analysis.R
#
# By default, the script searches common locations for analysis_master.rds and
# simulated_analysis_master.rds. A path can be supplied explicitly:
#
#   MAIHDA_DATA_PATH=/full/path/analysis_master.rds \
#     Rscript run_maihda_analysis.R
#
# A smaller validation run, using the 365-day any-event outcome while retaining
# both stratum specifications, can be requested with:
#
#   MAIHDA_TEST_MODE=1 Rscript run_maihda_analysis.R
#
# All tunable settings are environment variables; see section 1.
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260814)

RUN_STARTED_AT <- Sys.time()

# =============================================================================
# 1. USER-EDITABLE SETTINGS
# =============================================================================

# Small helpers so every setting can be supplied from the environment without
# editing the file. This matters on a computing cluster, where the script is
# usually read-only and parameters arrive through the job submission.
env_flag <- function(name, default = FALSE) {
  raw <- tolower(Sys.getenv(name, unset = if (default) "1" else "0"))
  raw %in% c("1", "true", "yes", "y", "on")
}
env_number <- function(name, default) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) return(default)
  parsed <- suppressWarnings(as.numeric(raw))
  if (is.na(parsed)) {
    stop("Environment variable ", name, " must be numeric, received: ", raw)
  }
  parsed
}

# Leave DATA_PATH empty to use the automatic file search below. An environment
# variable is convenient when running from Terminal or a computing cluster.
DATA_PATH <- Sys.getenv("MAIHDA_DATA_PATH", unset = "")

# All result files are placed here. A relative path is interpreted from the
# directory in which R is started.
OUTPUT_ROOT <- Sys.getenv("MAIHDA_OUTPUT_ROOT", unset = "maihda_outputs")

# --- Estimation engine ------------------------------------------------------
# "bayesian" (default) fits the models with brms/Stan by Hamiltonian Monte
# Carlo, following the Bayesian companion code to the Evans et al. tutorial.
# "mle" fits them with lme4 by maximum likelihood.
#
# Bayesian is the default because it is what the tutorial's own R code uses for
# the headline analysis, and because it resolves the one real weakness of the
# maximum-likelihood route for this model class: under MLE there is no
# straightforward way to obtain intervals for the derived quantities that
# matter most here -- the VPC, the PCV, the total predicted risk per stratum,
# and the absolute risk due to interaction on the probability scale. Those all
# require propagating uncertainty through a nonlinear transform of both the
# fixed and the random parts. The MLE route needs a zero-covariance assumption
# and an approximation to get there; the posterior gives them directly and
# exactly, as quantiles of draws.
#
# The MLE engine is retained deliberately. It is far faster, which makes it the
# right choice while iterating, and it needs no C++ toolchain, which matters in
# a locked-down environment where Stan may not be installable. Both engines
# produce identical output columns, so every table, figure and generated
# manuscript sentence works unchanged whichever is used. The engine used is
# recorded in the metrics, the manifest and the generated Methods.
ESTIMATION_ENGINE <- tolower(Sys.getenv("MAIHDA_ENGINE", unset = "bayesian"))
if (!ESTIMATION_ENGINE %in% c("bayesian", "mle")) {
  stop("MAIHDA_ENGINE must be \"bayesian\" or \"mle\", received: ",
       ESTIMATION_ENGINE)
}

# MCMC settings. The defaults are those used in the tutorial's Bayesian code:
# four chains of 2000 iterations with 1000 warmup.
MCMC_CHAINS <- as.integer(env_number("MAIHDA_MCMC_CHAINS", 4))
MCMC_ITERATIONS <- as.integer(env_number("MAIHDA_MCMC_ITER", 2000))
MCMC_WARMUP <- as.integer(env_number("MAIHDA_MCMC_WARMUP", 1000))
MCMC_SEED <- as.integer(env_number("MAIHDA_MCMC_SEED", 1))
MCMC_CORES <- as.integer(env_number("MAIHDA_MCMC_CORES",
                                    min(MCMC_CHAINS, parallel::detectCores())))
# Raised from the Stan default of 0.8 because the funnel geometry of a
# random-intercept model with many small groups is a classic source of
# divergent transitions.
MCMC_ADAPT_DELTA <- env_number("MAIHDA_MCMC_ADAPT_DELTA", 0.95)
MCMC_MAX_TREEDEPTH <- as.integer(env_number("MAIHDA_MCMC_MAX_TREEDEPTH", 12))

# Diagnostic thresholds. The tutorial asks for effective sample sizes above
# 400 and, conventionally, R-hat below 1.01.
MIN_EFFECTIVE_SAMPLE_SIZE <- env_number("MAIHDA_MIN_ESS", 400)
MAX_RHAT <- env_number("MAIHDA_MAX_RHAT", 1.01)

# The tutorial sets a normal(0, 1) prior on the fixed-effect coefficients of
# the logistic additive model, because brms's default flat prior is a poor
# choice on the logit scale: it places most of its mass near probabilities of
# zero and one. Intercept and random-effect priors are left at the brms
# defaults, which are weakly informative and sensible.
PRIOR_FIXED_EFFECTS <- Sys.getenv("MAIHDA_PRIOR_FIXED", unset = "normal(0, 1)")

# The tutorial uses posterior medians rather than means as point estimates.
POSTERIOR_ROBUST <- env_flag("MAIHDA_POSTERIOR_ROBUST", default = TRUE)

# Applies only to the MLE engine. nAGQ = 1 requests the standard Laplace
# approximation in glmer; 0 is faster but less accurate.
N_AGQ <- as.integer(env_number("MAIHDA_NAGQ", 1))

# Set through MAIHDA_TEST_MODE=1 for a one-outcome validation run.
TEST_MODE <- env_flag("MAIHDA_TEST_MODE")

# --- Minimum cell size ------------------------------------------------------
# Strata smaller than MIN_STRATUM_N, or with fewer than MIN_STRATUM_EVENTS
# events, are excluded before fitting.
#
# On what the best practice actually is, because the two relevant literatures
# pull in opposite directions and it is worth being explicit about the choice:
#
#   * The MAIHDA methodological literature does NOT recommend trimming. Partial
#     pooling is the whole point of the approach: a stratum of four people
#     contributes very little to its own estimate and is shrunk towards the
#     additive prediction, so it does no harm. Simulation work on MAIHDA
#     (Bell, Holman and Jones) finds the variance estimates hold up well even
#     with many small strata. On purely statistical grounds the defensible
#     answer is to keep everything.
#
#   * Statistical disclosure control in UK health data does require a floor.
#     A threshold of 10 is the conventional minimum for anything that will be
#     published at cell level, and this analysis writes stratum-level tables
#     and names individual strata on its figures.
#
# The default of 10 is therefore chosen to satisfy disclosure control at the
# lowest threshold that does so, rather than to improve the model. It is the
# defensible middle: a stricter floor such as 30 discards a large share of a
# five-way grid for no statistical gain. The sensitivity sweep in section 9
# refits the models across a range of thresholds so the choice can be shown to
# be immaterial to the conclusions rather than merely asserted to be.
#
# MIN_STRATUM_EVENTS defaults to 0, meaning zero-event strata are RETAINED.
# This is intentional. A stratum with 40 people and no events is a genuine
# observation of low risk, and partial pooling handles it correctly. Requiring
# at least one event would systematically discard the low-risk end of the
# distribution and bias the predicted-risk range upwards. Raise it only if
# disclosure control requires it.
MIN_STRATUM_N <- as.integer(env_number("MAIHDA_MIN_STRATUM_N", 10))
MIN_STRATUM_EVENTS <- as.integer(env_number("MAIHDA_MIN_STRATUM_EVENTS", 0))

# An analysis is abandoned, rather than fitted on an uninformative remnant, if
# fewer than this many strata survive the rule.
MIN_RETAINED_STRATA <- as.integer(env_number("MAIHDA_MIN_RETAINED_STRATA", 20))

# A threshold sweep is run for one outcome so that the choice above can be
# justified rather than asserted. It refits Models A and B at each candidate
# threshold and reports how the variance measures respond.
RUN_THRESHOLD_SENSITIVITY <- env_flag("MAIHDA_THRESHOLD_SENSITIVITY",
                                      default = TRUE)
SENSITIVITY_THRESHOLDS <- c(0, 5, 10, 20, 30, 50)

# The sweep exists to show that the trimming threshold does not drive the
# conclusions, which is a question about the variance components rather than
# about posterior uncertainty. It therefore always runs under maximum
# likelihood, even when the main analysis is Bayesian: running it under MCMC
# would multiply an overnight run by another six fits per specification for no
# gain in what the sweep is actually establishing. This is stated in the
# generated Methods.
SENSITIVITY_ENGINE <- "mle"

# --- Partially adjusted single-axis models ----------------------------------
# An OPTIONAL EXTENSION, not part of the Evans et al. tutorial sequence.
#
# A note on naming, because two MAIHDA lineages number their models
# incompatibly. In the Evans et al. tutorial the logistic models are Model 2A
# (null) and Model 2B (additive main effects) -- what this script calls Model A
# and Model B. In the Merlo/Persmark lineage the sequence is Model 1 (null),
# Model 2 (partially adjusted, one axis at a time) and Model 3 (fully
# adjusted). "Model 2" therefore means completely different things in the two
# traditions, so this script avoids the label entirely and describes these as
# partially adjusted single-axis models.
#
# What they do: one model per axis, each adding a single axis to the null
# model, so the PCV isolates that axis's contribution. This says something
# Model A and Model B cannot -- a PCV of 97% could be almost entirely frailty
# or evenly spread across five axes.
#
# IMPORTANT CAVEAT. The Evans et al. tutorial (section 2.4.C.i) explicitly
# warns that attending to individual axis contributions "seems to result in
# reversion to single-axis thinking about inequity (e.g., asking whether the
# effect of race(ism) is more important than income inequality) ... which is
# counter to the stated purpose of intersectional comparisons". This output is
# therefore secondary and descriptive. It should not displace the collective
# additive effect and the total stratum predictions, which are the point of the
# analysis. The caveat is carried into the generated manuscript text.
#
# One model is fitted per axis per outcome, so this multiplies the fitting work
# by roughly the number of axes. It is the first thing to switch off when
# iterating.
# Default OFF under Bayesian estimation. It multiplies the number of models by
# the number of axes, which under MCMC on a slow processor is hours of extra
# work for a secondary output that the tutorial explicitly cautions against
# over-reading. Under MLE it is cheap and stays on.
RUN_AXIS_DECOMPOSITION <- env_flag("MAIHDA_AXIS_DECOMPOSITION",
                                   default = ESTIMATION_ENGINE != "bayesian")

# --- Uncertainty intervals for the variance components ----------------------
# Published MAIHDA analyses report the VPC with an interval, which Bayesian
# estimation supplies directly from the posterior. Under maximum likelihood the
# equivalent is a profile-likelihood interval on the random-effect standard
# deviation, transformed to the variance, VPC and MOR scales. All three
# transformations are monotonic, so the interval endpoints carry across
# directly.
#
# Profiling costs a further optimisation per model and can be slow on the
# largest specifications; it fails gracefully to NA rather than stopping a run.
RUN_VARIANCE_INTERVALS <- env_flag("MAIHDA_VARIANCE_INTERVALS", default = TRUE)

# --- Numerical and output settings ------------------------------------------
# Tolerance used consistently for every singularity test, so that the metrics
# table, the log and the figure captions cannot disagree with one another.
SINGULAR_TOLERANCE <- 1e-5

# Published VPC values for health outcomes, used only to contextualise the
# results. Range and median across the 21 applied MAIHDA studies reviewed by
# Keller and colleagues (2023), Educational Psychology Review 35:31.
VPC_BENCHMARK_LOW <- 0.5
VPC_BENCHMARK_HIGH <- 41.9
VPC_BENCHMARK_MEDIAN <- 5.5

# Multiplicity control for the stratum interaction residuals. Twelve outcomes
# times two specifications times several hundred strata generates a great many
# simultaneous comparisons, and the uncorrected flag alone will produce false
# positives at a predictable rate. Both flags are reported.
FDR_LEVEL <- env_number("MAIHDA_FDR_LEVEL", 0.05)

# Vector output is useful for journal submission but doubles the figure writing
# time and needs the svglite package. It degrades to PNG-only automatically.
WRITE_SVG <- env_flag("MAIHDA_WRITE_SVG", default = TRUE)

# Fitted model objects are large, and a brms fit carrying several thousand
# posterior draws for every stratum is very much larger than a glmer fit. The
# default is therefore OFF under Bayesian estimation: everything needed for
# tables, figures and the manuscript is extracted to CSV before the object is
# discarded, so keeping the fits is a convenience rather than a requirement.
SAVE_MODELS <- env_flag("MAIHDA_SAVE_MODELS",
                        default = ESTIMATION_ENGINE != "bayesian")

# Skip an analysis whose metrics file already exists.
#
# The default is ON. A Bayesian run over many outcomes takes hours, and this
# makes an interrupted run resumable simply by starting it again: completed
# outcomes are skipped and the run picks up where it stopped. On a fresh output
# directory it does nothing at all. Set MAIHDA_RESUME=0 to force everything to
# be refitted.
RESUME <- env_flag("MAIHDA_RESUME", default = TRUE)

# --- Resource management for long unattended runs ---------------------------
# Fitted models are discarded and the garbage collector is run after every
# analysis, so peak memory is set by the single largest model rather than by
# the number of outcomes. Trimming additionally drops the parts of a brms fit
# that are not needed once the summaries have been extracted.
TRIM_FITS <- env_flag("MAIHDA_TRIM_FITS", default = TRUE)

# Stan compilation takes appreciable time on a slow processor and is pure
# overhead when the same model structure is refitted for each outcome. When
# enabled, the compiled Stan program from the first fit of each structure is
# reused for every subsequent outcome, which removes one compilation per model
# per outcome. It falls back to a fresh compile automatically if the design
# matrix changes.
REUSE_COMPILED_MODELS <- env_flag("MAIHDA_REUSE_COMPILED", default = TRUE)

# Logged after each analysis so a run that is drifting towards trouble is
# visible in the log rather than discovered when it stops.
REPORT_MEMORY <- env_flag("MAIHDA_REPORT_MEMORY", default = TRUE)

# Number of extreme strata retained at each end of the predicted-risk ranking.
N_EXTREME <- as.integer(env_number("MAIHDA_N_EXTREME", 6))

# Number of strata labelled at each end of the dense caterpillar and
# interaction panels. The significant-only interaction panel names every
# retained stratum on its y-axis, so this limit does not apply there.
N_PLOT_LABEL <- as.integer(env_number("MAIHDA_N_PLOT_LABEL", 1))

# Figure dimensions and resolution.
FIGURE_WIDTH <- env_number("MAIHDA_FIGURE_WIDTH", 9)
FIGURE_HEIGHT <- env_number("MAIHDA_FIGURE_HEIGHT", 6.5)
FIGURE_DPI <- env_number("MAIHDA_FIGURE_DPI", 320)

# Size of the stratum condition labels. These name the actual intersectional
# combinations and are the part of the figure a reader most needs to be able to
# read, so they are set well above the ggplot2 annotation default and are
# exposed here rather than buried in each geom.
#
# LABEL_SIZE_CALLOUT applies to the repelled callouts on the dense caterpillar
# and interaction panels; LABEL_SIZE_AXIS applies to the named y-axis of the
# significant-interaction forest panel. Sizes are in millimetres, ggplot2's
# convention for text geoms, where 1 mm is roughly 2.85 points.
LABEL_SIZE_CALLOUT <- env_number("MAIHDA_LABEL_SIZE_CALLOUT", 3.4)
LABEL_SIZE_AXIS <- env_number("MAIHDA_LABEL_SIZE_AXIS", 9.5)

# Preferred typeface, in order. The first family actually present on the system
# is used; if none are found the device default is used rather than allowing
# ggplot2 to emit a font warning for every single panel.
PREFERRED_FONTS <- c("Arial", "Helvetica", "Liberation Sans", "DejaVu Sans")

# All binary outcomes supported by analysis_master. Composite outcomes are
# included alongside their individual components because they answer distinct
# substantive questions.
OUTCOME_LABELS <- c(
  event_90d_any = "Any event, 90 days",
  event_365d_any = "Any event, 365 days",
  event_90d_nonfatal = "Non-fatal event, 90 days",
  event_365d_nonfatal = "Non-fatal event, 365 days",
  event_90d_fall = "Fall, 90 days",
  event_365d_fall = "Fall, 365 days",
  event_90d_fracture = "Fracture, 90 days",
  event_365d_fracture = "Fracture, 365 days",
  event_90d_delirium = "Delirium, 90 days",
  event_365d_delirium = "Delirium, 365 days",
  event_90d_death = "Death, 90 days",
  event_365d_death = "Death, 365 days"
)

if (TEST_MODE) {
  OUTCOME_LABELS <- OUTCOME_LABELS["event_365d_any"]
}

# Compact but fully interpretable axis names for plots and tables. Axis names
# are retained alongside their values because a string such as
# "85+ | Female | Q5" is much easier to misread when several variables have
# overlapping category labels. Defined here, with the other settings, because
# the descriptive tables need them before any figure code is reached.
axis_display_names <- c(
  age_band = "Age", age_3cat = "Age", sex = "Sex",
  ethnicity = "Ethnicity", ethnicity_4cat = "Ethnicity",
  imd_quintile = "IMD", imd_3cat = "IMD",
  opioid_strength = "Opioid", efi_category = "Frailty",
  efi_category_3 = "Frailty"
)
display_name_for <- function(axis) {
  known <- unname(axis_display_names[axis])
  ifelse(is.na(known), axis, known)
}

# The outcome used for the minimum-cell threshold sweep.
SENSITIVITY_OUTCOME <- Sys.getenv(
  "MAIHDA_SENSITIVITY_OUTCOME", unset = names(OUTCOME_LABELS)[1L]
)

# The variable names and their intended reference categories are defined once
# here so that the model specification, output labels and stratum identifiers
# remain consistent throughout the script.
#
# The intersectional strata are defined by five social and clinical axes: age,
# sex, ethnicity, deprivation and frailty. Opioid strength is deliberately NOT
# a stratum axis. It is a treatment characteristic rather than a social
# position, and including it would mix an exposure into the definition of the
# groups whose inequality is being measured.
#
# The structure remains a named list so that further specifications can be
# added without touching any other part of the script; every loop, table and
# figure iterates over whatever is defined here.
SPECIFICATIONS <- list(
  detailed = list(
    label = "Detailed categories",
    axes = c("age_band", "sex", "ethnicity", "imd_quintile", "efi_category"),
    references = c(
      age_band = "65-69",
      sex = "Male",
      ethnicity = "White",
      imd_quintile = "1 (Least deprived)",
      efi_category = "Fit"
    )
  ),
  # The collapsed specification trades resolution for precision. Coarser
  # categories produce fewer, larger strata, which tightens every stratum
  # estimate and reduces the share of the cohort lost to the minimum cell rule.
  # Running both and comparing is the honest way to show that the conclusions
  # are not an artefact of how finely the axes were cut.
  #
  # Frailty is left uncollapsed here so that exactly one thing differs per
  # axis. To collapse it too, change efi_category to efi_category_3 and its
  # reference to the corresponding level.
  collapsed = list(
    label = "Collapsed categories",
    axes = c("age_3cat", "sex", "ethnicity_4cat", "imd_3cat", "efi_category"),
    references = c(
      age_3cat = "65-74",
      sex = "Male",
      ethnicity_4cat = "White",
      imd_3cat = "Q1-Q2",
      efi_category = "Fit"
    )
  )
)

# --- Outcome-specific observability ------------------------------------------
# Each outcome has an observability flag naming the records for which that
# outcome could actually have been ascertained. Analysing an outcome on records
# where it was unobservable would treat "not observed because the person was
# not under observation" as "did not happen", which biases every risk downwards
# and does so unevenly across strata.
#
# The flag name is derived from the outcome name by the rule
#   event_<window>d_<type>  ->  observable_<window>_<type>
# so event_365d_fall uses observable_365_fall. Any outcome whose flag does not
# follow that pattern can be named explicitly in OUTCOME_ELIGIBILITY_OVERRIDES.
derive_eligibility_column <- function(outcome) {
  sub("^event_([0-9]+)d_(.+)$", "observable_\\1_\\2", outcome)
}
OUTCOME_ELIGIBILITY_OVERRIDES <- c(
  # example: event_90d_any = "observable_90_any_custom"
)

# Set to FALSE only to reproduce an earlier unfiltered run; it is not a
# defensible analysis choice.
APPLY_OBSERVABILITY_FILTER <- env_flag("MAIHDA_OBSERVABILITY_FILTER",
                                       default = TRUE)

# =============================================================================
# 2. PACKAGE AND FILE SET-UP
# =============================================================================

required_packages <- c(
  "lme4", "dplyr", "tidyr", "ggplot2", "patchwork", "scales",
  "ggrepel", "officer", "flextable"
)

if (ESTIMATION_ENGINE == "bayesian") {
  required_packages <- c(required_packages, "brms", "rstan")
}

# Missing packages are installed automatically unless this is switched off.
# In a trusted research environment there is often no route to CRAN, so the
# failure has to be informative rather than a bare error: the script reports
# precisely what is missing and how to obtain it.
AUTO_INSTALL <- env_flag("MAIHDA_AUTO_INSTALL", default = TRUE)
CRAN_MIRROR <- Sys.getenv("MAIHDA_CRAN_MIRROR",
                          unset = "https://cloud.r-project.org")

find_missing_packages <- function(packages) {
  packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
}

install_missing_packages <- function(packages, context) {
  if (length(packages) == 0L) return(character(0))
  message("Installing missing ", context, " packages: ",
          paste(packages, collapse = ", "))
  # A writable user library is created if the default library is read-only,
  # which is the usual arrangement on a managed analysis machine.
  library_path <- .libPaths()[1L]
  if (file.access(library_path, mode = 2L) != 0L) {
    library_path <- Sys.getenv("R_LIBS_USER",
                               unset = file.path("~", "R", "library"))
    library_path <- path.expand(library_path)
    dir.create(library_path, recursive = TRUE, showWarnings = FALSE)
    .libPaths(c(library_path, .libPaths()))
    message("Default library is not writable; installing into ", library_path)
  }
  for (package in packages) {
    try(
      utils::install.packages(package, lib = library_path,
                              repos = CRAN_MIRROR, quiet = TRUE),
      silent = TRUE
    )
  }
  find_missing_packages(packages)
}

missing_packages <- find_missing_packages(required_packages)
if (length(missing_packages) > 0L && AUTO_INSTALL) {
  missing_packages <- install_missing_packages(missing_packages, "required")
}

if (length(missing_packages) > 0L) {
  stop(
    "The following required packages are missing and could not be installed: ",
    paste(missing_packages, collapse = ", "),
    "\n\nIf this machine has no route to CRAN, install them from a local ",
    "repository or transfer the sources, then run again. To install ",
    "manually:\n  install.packages(c(",
    paste(sprintf("\"%s\"", missing_packages), collapse = ", "), "))",
    if ("brms" %in% missing_packages || "rstan" %in% missing_packages) {
      paste0(
        "\n\nbrms and rstan additionally need a working C++ toolchain and the ",
        "BH, RcppEigen and StanHeaders packages. If Stan cannot be installed ",
        "here, run with MAIHDA_ENGINE=mle to use maximum likelihood instead; ",
        "all outputs are produced either way."
      )
    } else "",
    "\n\nSet MAIHDA_AUTO_INSTALL=0 to disable automatic installation."
  )
}

# Optional packages improve the output but never block a run.
optional_packages <- c("ggprism", "svglite", "systemfonts")
missing_optional <- find_missing_packages(optional_packages)
if (length(missing_optional) > 0L && AUTO_INSTALL) {
  invisible(install_missing_packages(missing_optional, "optional"))
}

# Optional packages are detected rather than required, so that the absence of a
# purely cosmetic dependency cannot stop an analysis run.
HAS_GGPRISM <- requireNamespace("ggprism", quietly = TRUE)
HAS_SVGLITE <- requireNamespace("svglite", quietly = TRUE)
HAS_SYSTEMFONTS <- requireNamespace("systemfonts", quietly = TRUE)
if (WRITE_SVG && !HAS_SVGLITE) WRITE_SVG <- FALSE

suppressPackageStartupMessages({
  library(lme4)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggrepel)
  library(officer)
  library(flextable)
})

# Determine the script location so that automatic input discovery also works
# when this file is launched from RStudio or from a different working directory.
command_arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", command_arguments, value = TRUE)
script_path <- if (length(file_argument) == 1L) {
  normalizePath(sub("^--file=", "", file_argument), mustWork = FALSE)
} else {
  normalizePath(getwd(), mustWork = FALSE)
}
script_directory <- if (length(file_argument) == 1L) {
  dirname(script_path)
} else {
  getwd()
}
project_directory <- normalizePath(file.path(script_directory, ".."),
                                   mustWork = FALSE)

find_data_file <- function(explicit_path = "") {
  if (nzchar(explicit_path) && !file.exists(explicit_path)) {
    stop("MAIHDA_DATA_PATH was set to a file that does not exist: ",
         explicit_path)
  }
  search_directories <- unique(c(
    getwd(), script_directory, project_directory,
    file.path(project_directory, "outputs"),
    file.path(project_directory, "data"),
    file.path(project_directory, "work"),
    file.path(project_directory, "work", "regeneration_check")
  ))
  candidates <- unique(c(
    explicit_path,
    as.vector(t(outer(
      search_directories,
      c("analysis_master.rds", "simulated_analysis_master.rds"),
      file.path
    )))
  ))
  candidates <- candidates[nzchar(candidates)]
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) {
    stop(
      "Could not find analysis_master.rds. Searched:\n  ",
      paste(candidates, collapse = "\n  "),
      "\nSet MAIHDA_DATA_PATH=/full/path/analysis_master.rds to be explicit."
    )
  }
  normalizePath(existing[1L])
}

DATA_PATH <- find_data_file(DATA_PATH)
OUTPUT_ROOT <- normalizePath(OUTPUT_ROOT, mustWork = FALSE)

directories <- c(
  OUTPUT_ROOT,
  file.path(OUTPUT_ROOT, "tables"),
  file.path(OUTPUT_ROOT, "tables", "word"),
  file.path(OUTPUT_ROOT, "tables", "per_analysis"),
  file.path(OUTPUT_ROOT, "figures"),
  file.path(OUTPUT_ROOT, "figures", "individual"),
  file.path(OUTPUT_ROOT, "figures", "multipanel"),
  file.path(OUTPUT_ROOT, "models"),
  file.path(OUTPUT_ROOT, "logs")
)
for (directory in directories) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(OUTPUT_ROOT)) {
  stop("Could not create the output directory: ", OUTPUT_ROOT)
}
# A write test now is far preferable to discovering a permissions problem after
# an hour of model fitting.
write_probe <- file.path(OUTPUT_ROOT, ".write_probe")
probe_status <- try(
  {
    writeLines("ok", write_probe)
    unlink(write_probe)
  },
  silent = TRUE
)
if (inherits(probe_status, "try-error")) {
  stop("The output directory is not writable: ", OUTPUT_ROOT)
}

# --- Logging ----------------------------------------------------------------
# Everything reported to the console is also written to a durable log, because
# the console output of an unattended cluster job is routinely lost.
LOG_PATH <- file.path(OUTPUT_ROOT, "logs", "run_log.txt")
cat("", file = LOG_PATH)

log_message <- function(..., level = "INFO") {
  line <- paste0(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC"), " UTC  ",
    formatC(level, width = -5), " ", paste0(..., collapse = "")
  )
  message(line)
  cat(line, "\n", file = LOG_PATH, sep = "", append = TRUE)
  invisible(line)
}

# Every file the script produces is recorded, so the run can be audited and so
# an incomplete run is immediately obvious.
output_registry <- new.env(parent = emptyenv())
output_registry$paths <- character(0)
register_output <- function(path) {
  output_registry$paths <- c(output_registry$paths, path)
  invisible(path)
}
write_csv_output <- function(data, path) {
  write.csv(data, path, row.names = FALSE, na = "")
  register_output(path)
}

log_message("MAIHDA analysis starting.")
log_message("Input file      : ", DATA_PATH)
log_message("Output directory: ", OUTPUT_ROOT)
log_message("Outcomes        : ", length(OUTCOME_LABELS),
            if (TEST_MODE) " (test mode)" else "")
log_message("Minimum cell    : n >= ", MIN_STRATUM_N,
            ", events >= ", MIN_STRATUM_EVENTS)
if (!HAS_GGPRISM) {
  log_message("ggprism is not installed; using the matched fallback theme.",
              level = "NOTE")
}
if (!WRITE_SVG) {
  log_message("SVG output disabled (svglite unavailable or turned off).",
              level = "NOTE")
}

log_message("Loading data ...")
analysis_master <- readRDS(DATA_PATH)

if (!is.data.frame(analysis_master)) {
  stop("The RDS object must be a data.frame or tibble, received: ",
       paste(class(analysis_master), collapse = "/"))
}
if (nrow(analysis_master) == 0L) {
  stop("The input data contains no rows.")
}
# tibbles and data.tables both subset differently from data.frame in ways that
# would silently change behaviour further down. Normalising once removes the
# whole class of problem.
analysis_master <- as.data.frame(analysis_master)

# =============================================================================
# 3. INPUT VALIDATION AND AUDIT TABLES
# =============================================================================

all_axes <- unique(unname(unlist(lapply(SPECIFICATIONS, `[[`, "axes"))))
all_required_variables <- unique(c(all_axes, names(OUTCOME_LABELS)))
missing_variables <- setdiff(all_required_variables, names(analysis_master))
if (length(missing_variables) > 0L) {
  stop("Required variables are absent from the input data: ",
       paste(missing_variables, collapse = ", "))
}

# Binary outcomes must resolve to 0, 1 or NA. Logical columns and factor or
# character columns holding "0"/"1" are accepted and coerced, because both
# appear routinely in derived datasets; anything else is a genuine error and is
# caught here rather than halfway through model fitting.
coerce_binary_outcome <- function(x, outcome) {
  if (is.logical(x)) return(as.integer(x))
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    trimmed <- trimws(x)
    recognised <- trimmed %in% c("0", "1") | is.na(trimmed)
    if (!all(recognised)) {
      stop(outcome, " is a character or factor column containing values other ",
           "than 0/1: ",
           paste(utils::head(setdiff(unique(trimmed), c("0", "1", NA)), 5),
                 collapse = ", "))
    }
    return(as.integer(trimmed))
  }
  if (!is.numeric(x)) {
    stop(outcome, " has an unsupported type: ", paste(class(x), collapse = "/"))
  }
  observed <- unique(x[!is.na(x)])
  if (!all(observed %in% c(0, 1))) {
    stop(outcome, " is not coded as a binary 0/1 variable. Observed values: ",
         paste(utils::head(sort(observed), 8), collapse = ", "))
  }
  as.integer(x)
}

for (outcome in names(OUTCOME_LABELS)) {
  analysis_master[[outcome]] <- coerce_binary_outcome(
    analysis_master[[outcome]], outcome
  )
  complete_values <- analysis_master[[outcome]][!is.na(analysis_master[[outcome]])]
  if (length(complete_values) == 0L) {
    stop(outcome, " is missing for every record.")
  }
  # A constant outcome cannot support a random-intercept logistic model, and
  # the failure mode is an obscure lme4 error rather than a clear one.
  if (length(unique(complete_values)) < 2L) {
    stop(outcome, " is constant (all ", unique(complete_values),
         "); a MAIHDA model cannot be fitted to it.")
  }
}

# Axis variables must be usable as categorical predictors. A continuous
# variable reaching this point is nearly always a data-linkage mistake, and
# stratifying on it would silently create one stratum per distinct value.
for (axis in all_axes) {
  distinct_n <- length(unique(analysis_master[[axis]][!is.na(analysis_master[[axis]])]))
  if (distinct_n < 2L) {
    stop("Axis ", axis, " has fewer than two observed categories.")
  }
  if (distinct_n > 50L) {
    stop("Axis ", axis, " has ", distinct_n, " distinct values, which is too ",
         "many for an intersectional stratum. Check that it is categorical.")
  }
}

# --- Resolve and validate the observability flags ---------------------------
OUTCOME_ELIGIBILITY <- vapply(names(OUTCOME_LABELS), function(outcome) {
  if (outcome %in% names(OUTCOME_ELIGIBILITY_OVERRIDES)) {
    unname(OUTCOME_ELIGIBILITY_OVERRIDES[outcome])
  } else {
    derive_eligibility_column(outcome)
  }
}, character(1))

if (APPLY_OBSERVABILITY_FILTER) {
  absent_flags <- OUTCOME_ELIGIBILITY[!OUTCOME_ELIGIBILITY %in% names(analysis_master)]
  if (length(absent_flags) > 0L) {
    stop(
      "Observability flags are missing from the data: ",
      paste(unique(absent_flags), collapse = ", "),
      "\nThese are derived from the outcome names by the rule ",
      "event_<window>d_<type> -> observable_<window>_<type>. If your columns ",
      "are named differently, set them in OUTCOME_ELIGIBILITY_OVERRIDES near ",
      "the top of the script, or set MAIHDA_OBSERVABILITY_FILTER=0 to disable ",
      "filtering entirely."
    )
  }
  # Coerced through the same routine as the outcomes, so logical, 0/1 numeric
  # and "0"/"1" character columns are all accepted and anything else is a
  # clear error rather than a silent misreading.
  for (flag in unique(OUTCOME_ELIGIBILITY)) {
    analysis_master[[flag]] <- coerce_binary_outcome(analysis_master[[flag]],
                                                     flag)
    if (all(analysis_master[[flag]] == 0, na.rm = TRUE)) {
      stop("Observability flag ", flag, " is FALSE for every record.")
    }
  }

  eligibility_audit <- bind_rows(lapply(names(OUTCOME_LABELS), function(outcome) {
    flag <- unname(OUTCOME_ELIGIBILITY[outcome])
    values <- analysis_master[[flag]]
    eligible <- !is.na(values) & values == 1
    outcome_values <- analysis_master[[outcome]]
    data.frame(
      outcome = outcome,
      outcome_label = unname(OUTCOME_LABELS[outcome]),
      eligibility_column = flag,
      total_n = length(values),
      observable_n = sum(eligible),
      unobservable_n = sum(!eligible),
      observable_percent = 100 * mean(eligible),
      # Events occurring among records flagged unobservable indicate the flag
      # and the outcome disagree, which is worth knowing before modelling.
      events_among_unobservable = sum(!eligible & !is.na(outcome_values) &
                                        outcome_values == 1)
    )
  }))
  write_csv_output(eligibility_audit,
                   file.path(OUTPUT_ROOT, "tables",
                             "outcome_observability_audit.csv"))

  for (row_index in seq_len(nrow(eligibility_audit))) {
    row <- eligibility_audit[row_index, ]
    log_message("Observability: ", row$outcome_label, " uses ",
                row$eligibility_column, "; ", format(row$observable_n,
                                                     big.mark = ","),
                " of ", format(row$total_n, big.mark = ","), " records (",
                formatC(row$observable_percent, format = "f", digits = 1),
                "%).")
    if (row$events_among_unobservable > 0L) {
      log_message("  ", row$events_among_unobservable, " events occur among ",
                  "records flagged unobservable for this outcome. The flag ",
                  "and the outcome disagree; check the derivation.",
                  level = "WARN")
    }
  }
} else {
  log_message("Observability filtering is DISABLED. Every outcome will use ",
              "the full cohort.", level = "WARN")
  OUTCOME_ELIGIBILITY <- setNames(rep(NA_character_, length(OUTCOME_LABELS)),
                                  names(OUTCOME_LABELS))
}

input_audit <- data.frame(
  item = c("input_file", "rows", "columns", "test_mode", "n_outcomes",
           "min_stratum_n", "min_stratum_events", "n_agq", "fdr_level",
           "r_version", "lme4_version", "analysis_timestamp_utc"),
  value = c(DATA_PATH, nrow(analysis_master), ncol(analysis_master),
            TEST_MODE, length(OUTCOME_LABELS), MIN_STRATUM_N,
            MIN_STRATUM_EVENTS, N_AGQ, FDR_LEVEL,
            R.version.string,
            as.character(utils::packageVersion("lme4")),
            format(Sys.time(), tz = "UTC"))
)
write_csv_output(input_audit, file.path(OUTPUT_ROOT, "tables",
                                        "input_audit.csv"))

variable_audit <- data.frame(
  variable = all_required_variables,
  class = vapply(analysis_master[all_required_variables],
                 function(x) class(x)[1L], character(1)),
  missing_n = vapply(analysis_master[all_required_variables],
                     function(x) sum(is.na(x)), numeric(1)),
  missing_percent = vapply(analysis_master[all_required_variables],
                           function(x) 100 * mean(is.na(x)), numeric(1)),
  distinct_n = vapply(analysis_master[all_required_variables],
                      function(x) length(unique(x[!is.na(x)])), numeric(1)),
  row.names = NULL
)
write_csv_output(variable_audit,
                 file.path(OUTPUT_ROOT, "tables", "variable_audit.csv"))

# A full category listing removes any ambiguity about which level was used as
# the reference and how each level is spelled in the source data.
category_audit <- bind_rows(lapply(all_axes, function(axis) {
  values <- analysis_master[[axis]]
  counts <- sort(table(values, useNA = "no"), decreasing = TRUE)
  reference_for <- vapply(SPECIFICATIONS, function(specification) {
    if (axis %in% specification$axes) {
      unname(specification$references[axis])
    } else {
      NA_character_
    }
  }, character(1))
  reference_level <- unique(reference_for[!is.na(reference_for)])
  data.frame(
    variable = axis,
    level = names(counts),
    n = as.numeric(counts),
    percent = 100 * as.numeric(counts) / sum(counts),
    is_reference = names(counts) %in% reference_level
  )
}))
write_csv_output(category_audit,
                 file.path(OUTPUT_ROOT, "tables", "category_audit.csv"))

outcome_counts <- bind_rows(lapply(names(OUTCOME_LABELS), function(outcome) {
  values <- analysis_master[[outcome]]
  complete <- !is.na(values)
  data.frame(
    outcome = outcome,
    outcome_label = unname(OUTCOME_LABELS[outcome]),
    total_n = length(values),
    complete_n = sum(complete),
    missing_n = sum(!complete),
    event_n = sum(values[complete] == 1),
    non_event_n = sum(values[complete] == 0),
    event_percent = 100 * mean(values[complete] == 1)
  )
}))
write_csv_output(outcome_counts,
                 file.path(OUTPUT_ROOT, "tables", "outcome_counts.csv"))

log_message("Validation passed: ", nrow(analysis_master), " rows, ",
            length(OUTCOME_LABELS), " outcomes, ",
            length(SPECIFICATIONS), " specifications.")

# =============================================================================
# 3B. DESCRIPTIVE STATISTICS
# =============================================================================
# The MAIHDA results are only interpretable against a clear description of the
# cohort they came from. This section produces the conventional descriptive
# material a reader or reviewer will expect before any modelling:
#
#   1. Cohort characteristics, one row per category of each stratum axis.
#   2. Crude event rates by category, with Wilson confidence intervals, for
#      every outcome and every axis category.
#   3. A complete stratum-level descriptive table, including the strata that
#      the minimum cell size rule will later remove, so nothing disappears
#      without being counted somewhere.
#   4. The stratum size distribution, which is the single most useful summary
#      of how sparse a given intersectional grid is.
#
# These are all computed on the unfiltered cohort. Post-filter descriptions are
# produced separately in the retention tables, so the two are never confused.

# The Wilson score interval is used in preference to the Wald interval because
# many categories here are small or have low event counts, exactly the
# situation in which the Wald interval misbehaves and can run outside [0, 1].
wilson_interval <- function(events, n, level = 0.95) {
  z <- qnorm(1 - (1 - level) / 2)
  proportion <- ifelse(n > 0, events / n, NA_real_)
  denominator <- 1 + z^2 / n
  centre <- (proportion + z^2 / (2 * n)) / denominator
  half_width <- z * sqrt(proportion * (1 - proportion) / n +
                           z^2 / (4 * n^2)) / denominator
  list(
    low = ifelse(n > 0, pmax(0, centre - half_width), NA_real_),
    high = ifelse(n > 0, pmin(1, centre + half_width), NA_real_)
  )
}

descriptive_axes <- SPECIFICATIONS[[1L]]$axes
descriptive_references <- SPECIFICATIONS[[1L]]$references

# --- 1. Cohort characteristics ----------------------------------------------
cohort_characteristics <- bind_rows(lapply(descriptive_axes, function(axis) {
  values <- analysis_master[[axis]]
  complete <- !is.na(values)
  counts <- table(values[complete])
  # Levels are reported in their natural factor order where one exists, so the
  # table reads in a sensible sequence rather than alphabetically.
  level_order <- if (is.factor(values)) levels(droplevels(values[complete])) else
    names(counts)
  data.frame(
    variable = axis,
    variable_label = unname(axis_display_names[axis]),
    level = level_order,
    n = as.numeric(counts[level_order]),
    percent = 100 * as.numeric(counts[level_order]) / sum(counts),
    is_reference = level_order == unname(descriptive_references[axis]),
    missing_n = sum(!complete),
    missing_percent = 100 * mean(!complete)
  )
}))
write_csv_output(cohort_characteristics,
                 file.path(OUTPUT_ROOT, "tables",
                           "descriptive_cohort_characteristics.csv"))

# --- 2. Crude event rates by category ---------------------------------------
descriptive_event_rates <- bind_rows(lapply(names(OUTCOME_LABELS),
                                            function(outcome) {
  bind_rows(lapply(descriptive_axes, function(axis) {
    frame <- data.frame(
      level = analysis_master[[axis]],
      outcome_value = analysis_master[[outcome]]
    )
    frame <- frame[!is.na(frame$level) & !is.na(frame$outcome_value), ]
    summarised <- frame %>%
      group_by(level) %>%
      summarise(n = n(), events = sum(outcome_value == 1), .groups = "drop")
    interval <- wilson_interval(summarised$events, summarised$n)
    data.frame(
      outcome = outcome,
      outcome_label = unname(OUTCOME_LABELS[outcome]),
      variable = axis,
      variable_label = unname(axis_display_names[axis]),
      level = as.character(summarised$level),
      n = summarised$n,
      events = summarised$events,
      event_percent = 100 * summarised$events / summarised$n,
      ci_low_percent = 100 * interval$low,
      ci_high_percent = 100 * interval$high,
      is_reference = as.character(summarised$level) ==
        unname(descriptive_references[axis])
    )
  }))
}))
write_csv_output(descriptive_event_rates,
                 file.path(OUTPUT_ROOT, "tables",
                           "descriptive_event_rates_by_category.csv"))

# --- 3. Stratum-level description -------------------------------------------
# Built once per specification on complete cases across the axes only. Outcome
# specific counts are added for every outcome, and the flag showing whether the
# stratum will survive the minimum cell size rule is attached here so that the
# full grid and the analytic grid can always be reconciled.
descriptive_strata <- bind_rows(lapply(names(SPECIFICATIONS),
                                       function(specification_name) {
  specification <- SPECIFICATIONS[[specification_name]]
  axes <- specification$axes
  complete <- complete.cases(analysis_master[, axes, drop = FALSE])
  axis_data <- droplevels(analysis_master[complete, axes, drop = FALSE])
  axis_data$stratum <- interaction(axis_data[, axes, drop = FALSE],
                                   drop = TRUE, sep = " | ")

  base <- axis_data %>%
    group_by(across(all_of(c(axes, "stratum")))) %>%
    summarise(n = n(), .groups = "drop")

  outcome_columns <- lapply(names(OUTCOME_LABELS), function(outcome) {
    values <- analysis_master[[outcome]][complete]
    aggregated <- tapply(values, axis_data$stratum,
                         function(x) sum(x == 1, na.rm = TRUE))
    aggregated[match(as.character(base$stratum), names(aggregated))]
  })
  names(outcome_columns) <- paste0("events__", names(OUTCOME_LABELS))

  bind_cols(
    data.frame(
      specification = specification_name,
      specification_label = specification$label
    ),
    base,
    as.data.frame(outcome_columns)
  ) %>%
    mutate(
      percent_of_cohort = 100 * n / sum(n),
      meets_minimum_cell = n >= MIN_STRATUM_N
    ) %>%
    arrange(desc(n))
}))
write_csv_output(descriptive_strata,
                 file.path(OUTPUT_ROOT, "tables",
                           "descriptive_strata.csv"))

# --- 4. Stratum size distribution -------------------------------------------
descriptive_stratum_distribution <- descriptive_strata %>%
  group_by(specification, specification_label) %>%
  summarise(
    observed_strata = n(),
    possible_strata = prod(vapply(
      SPECIFICATIONS[[specification[1L]]]$axes,
      function(axis) length(unique(analysis_master[[axis]][
        !is.na(analysis_master[[axis]])
      ])),
      numeric(1)
    )),
    individuals = sum(n),
    minimum_n = min(n),
    q1_n = unname(quantile(n, 0.25)),
    median_n = median(n),
    mean_n = mean(n),
    q3_n = unname(quantile(n, 0.75)),
    maximum_n = max(n),
    strata_below_10 = sum(n < 10),
    strata_10_to_29 = sum(n >= 10 & n < 30),
    strata_30_to_99 = sum(n >= 30 & n < 100),
    strata_100_plus = sum(n >= 100),
    individuals_in_strata_below_minimum = sum(n[n < MIN_STRATUM_N]),
    percent_individuals_below_minimum =
      100 * sum(n[n < MIN_STRATUM_N]) / sum(n),
    .groups = "drop"
  ) %>%
  mutate(empty_strata = possible_strata - observed_strata)
write_csv_output(descriptive_stratum_distribution,
                 file.path(OUTPUT_ROOT, "tables",
                           "descriptive_stratum_distribution.csv"))

for (row_index in seq_len(nrow(descriptive_stratum_distribution))) {
  row <- descriptive_stratum_distribution[row_index, ]
  log_message("Descriptives (", row$specification_label, "): ",
              row$observed_strata, " of ", row$possible_strata,
              " possible strata observed; median n = ", row$median_n,
              "; ", formatC(row$percent_individuals_below_minimum,
                            format = "f", digits = 1),
              "% of individuals sit in strata below n = ", MIN_STRATUM_N, ".")
}

# =============================================================================
# 4. REUSABLE STATISTICAL FUNCTIONS
# =============================================================================

inv_logit <- plogis

# Weighted AUC for collapsed binomial data. Each row represents `events`
# positive observations and `non_events` negative observations at one score.
# Scores are first grouped so ties receive the standard half-credit exactly.
weighted_auc <- function(events, non_events, score) {
  keep <- is.finite(score) & !is.na(events) & !is.na(non_events)
  if (!any(keep)) return(NA_real_)
  auc_data <- data.frame(
    score = score[keep],
    events = as.numeric(events[keep]),
    non_events = as.numeric(non_events[keep])
  ) %>%
    group_by(score) %>%
    summarise(
      events = sum(events),
      non_events = sum(non_events),
      .groups = "drop"
    ) %>%
    arrange(score)

  # Totals are held in double precision throughout. With a cohort of this size
  # the number of possible case-control pairs can exceed R's 32-bit integer
  # limit even though the resulting AUC remains well behaved.
  total_events <- sum(auc_data$events)
  total_non_events <- sum(auc_data$non_events)
  if (total_events == 0 || total_non_events == 0) return(NA_real_)

  non_events_below <- c(0, head(cumsum(auc_data$non_events), -1L))
  concordant_pairs <- sum(
    auc_data$events * (non_events_below + 0.5 * auc_data$non_events)
  )
  concordant_pairs / (total_events * total_non_events)
}

# Extract the random-intercept variance. A variance can legitimately be zero in
# Model B if additive effects explain all between-stratum variation. The result
# is kept rather than replaced with an arbitrary positive value.
random_intercept_variance <- function(model) {
  variance_table <- as.data.frame(VarCorr(model))
  matched <- variance_table$vcov[variance_table$grp == "stratum"]
  if (length(matched) == 0L) return(NA_real_)
  matched[1L]
}

# Latent-response VPC for logistic models. The individual-level variance is
# fixed to pi^2/3, the variance of the standard logistic distribution.
logistic_vpc <- function(variance) {
  variance / (variance + (pi^2 / 3))
}

# Median odds ratio translates the random-intercept variance to the odds-ratio
# scale. It is the median contrast between two otherwise identical individuals
# drawn from higher- and lower-risk strata.
median_odds_ratio <- function(variance) {
  exp(qnorm(0.75) * sqrt(2 * variance))
}

# Profile-likelihood interval for the random-intercept variance, and for the
# VPC and MOR derived from it. lme4 profiles the standard deviation on the
# theta scale; because variance, VPC and MOR are all monotonic increasing
# functions of that standard deviation, the interval endpoints transform
# directly without needing a delta-method approximation.
#
# A variance at the zero boundary has a lower limit of zero by construction,
# which is correct and is reported as such rather than suppressed.
variance_interval <- function(model, level = 0.95) {
  empty <- list(
    variance_low = NA_real_, variance_high = NA_real_,
    vpc_low = NA_real_, vpc_high = NA_real_,
    mor_low = NA_real_, mor_high = NA_real_,
    interval_method = "not computed"
  )
  if (!RUN_VARIANCE_INTERVALS) return(empty)

  profiled <- tryCatch(
    suppressWarnings(suppressMessages(
      confint(model, parm = "theta_", method = "profile", level = level,
              oldNames = FALSE)
    )),
    error = function(condition) NULL
  )
  if (is.null(profiled) || !is.matrix(profiled) || nrow(profiled) == 0L) {
    empty$interval_method <- "profile failed"
    return(empty)
  }

  standard_deviations <- as.numeric(profiled[1L, ])
  if (anyNA(standard_deviations)) {
    empty$interval_method <- "profile failed"
    return(empty)
  }
  variances <- pmax(0, standard_deviations^2)
  list(
    variance_low = variances[1L], variance_high = variances[2L],
    vpc_low = 100 * logistic_vpc(variances[1L]),
    vpc_high = 100 * logistic_vpc(variances[2L]),
    mor_low = median_odds_ratio(variances[1L]),
    mor_high = median_odds_ratio(variances[2L]),
    interval_method = "profile likelihood"
  )
}

extract_convergence_message <- function(model) {
  messages <- model@optinfo$conv$lme4$messages
  if (is.null(messages) || length(messages) == 0L) "OK" else
    paste(messages, collapse = " | ")
}

# --- Robust model fitting ---------------------------------------------------
# A single optimiser occasionally stalls on a flat likelihood and reports a
# zero variance that a different optimiser would not. Because a zero variance
# is a substantive conclusion in MAIHDA, it is worth the extra fits to
# establish whether the boundary is genuine. The ladder is only walked when the
# first fit is singular or raises a convergence warning, so a healthy run pays
# nothing for it.
OPTIMIZER_LADDER <- list(
  list(name = "bobyqa",
       control = list(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))),
  list(name = "Nelder_Mead",
       control = list(optimizer = "Nelder_Mead",
                      optCtrl = list(maxfun = 200000))),
  list(name = "nlminbwrap",
       control = list(optimizer = "nlminbwrap", optCtrl = list()))
)

build_control <- function(entry) {
  glmerControl(
    optimizer = entry$control$optimizer,
    optCtrl = entry$control$optCtrl,
    # lme4's own singularity warning is suppressed here so that the script can
    # report the result with its interpretation attached instead of emitting a
    # bare warning that reads like a failure. Singularity is still detected,
    # recorded in the metrics and written to the log.
    check.conv.singular = .makeCC(action = "ignore", tol = SINGULAR_TOLERANCE)
  )
}

fit_glmer_robust <- function(model_formula, data, model_label) {
  attempts <- list()
  best <- NULL

  for (index in seq_along(OPTIMIZER_LADDER)) {
    entry <- OPTIMIZER_LADDER[[index]]
    captured_warnings <- character(0)
    fitted <- withCallingHandlers(
      tryCatch(
        glmer(model_formula, data = data, family = binomial(link = "logit"),
              nAGQ = N_AGQ, control = build_control(entry)),
        error = function(condition) condition
      ),
      warning = function(condition) {
        captured_warnings <<- c(captured_warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }
    )

    if (inherits(fitted, "condition")) {
      attempts[[length(attempts) + 1L]] <- data.frame(
        optimizer = entry$name, status = "error",
        log_likelihood = NA_real_, variance = NA_real_, singular = NA,
        detail = conditionMessage(fitted)
      )
      next
    }

    singular <- isSingular(fitted, tol = SINGULAR_TOLERANCE)
    attempts[[length(attempts) + 1L]] <- data.frame(
      optimizer = entry$name, status = "fitted",
      log_likelihood = as.numeric(logLik(fitted)),
      variance = random_intercept_variance(fitted),
      singular = singular,
      detail = paste(c(extract_convergence_message(fitted), captured_warnings),
                     collapse = " | ")
    )

    if (is.null(best) ||
        as.numeric(logLik(fitted)) > as.numeric(logLik(best$model)) + 1e-8) {
      best <- list(model = fitted, optimizer = entry$name, singular = singular,
                   warnings = captured_warnings)
    }

    # A clean, non-singular fit needs no further work.
    converged_cleanly <- !singular &&
      length(captured_warnings) == 0L &&
      identical(extract_convergence_message(fitted), "OK")
    if (converged_cleanly) break
  }

  if (is.null(best)) {
    stop("Model ", model_label, " could not be fitted with any optimiser: ",
         paste(vapply(attempts, function(a) paste0(a$optimizer, ": ", a$detail),
                      character(1)), collapse = " || "))
  }

  attempt_table <- bind_rows(attempts)
  successful <- attempt_table[attempt_table$status == "fitted", , drop = FALSE]
  # The boundary is treated as confirmed when every optimiser that produced a
  # fit agreed on it. A disagreement is recorded so it can be inspected.
  boundary_confirmed <- nrow(successful) > 0L && all(successful$singular)

  list(
    model = best$model,
    optimizer = best$optimizer,
    singular = best$singular,
    boundary_confirmed = boundary_confirmed,
    n_optimizers_tried = nrow(attempt_table),
    attempts = attempt_table
  )
}

# =============================================================================
# 4B. BAYESIAN ESTIMATION ENGINE
# =============================================================================
# Follows the Bayesian companion code to the Evans et al. tutorial (brms/Stan,
# adapted by Webb and Bell). Four things differ from the maximum-likelihood
# route, and all four are improvements for this model class:
#
#   1. The VPC and PCV are computed from the posterior draws of the
#      random-intercept standard deviation, so they carry credible intervals
#      directly rather than needing a profile-likelihood approximation.
#   2. The absolute risk and the absolute risk due to interaction are computed
#      per draw on the probability scale and then summarised, so uncertainty is
#      propagated exactly through the inverse-logit transform. The MLE route
#      has to assume zero covariance between the fixed and random parts and
#      approximate; here no such assumption is needed.
#   3. Point estimates are posterior medians, as in the tutorial code.
#   4. Convergence is a property to be checked rather than assumed, so R-hat,
#      effective sample sizes and divergent transitions are recorded for every
#      model and surfaced in the outputs.

stan_environment_ready <- function() {
  if (!requireNamespace("brms", quietly = TRUE)) {
    return(list(ready = FALSE, reason = "the brms package is not installed"))
  }
  if (!requireNamespace("rstan", quietly = TRUE)) {
    return(list(ready = FALSE, reason = "the rstan package is not installed"))
  }
  boost_include <- system.file("include", package = "BH")
  if (!nzchar(boost_include)) {
    return(list(
      ready = FALSE,
      reason = paste(
        "the BH package provides no headers, so Stan cannot compile.",
        "Install the CRAN BH package (not a distribution shim)"
      )
    ))
  }
  list(ready = TRUE, reason = "")
}

if (ESTIMATION_ENGINE == "bayesian") {
  stan_status <- stan_environment_ready()
  if (!stan_status$ready) {
    stop(
      "Bayesian estimation was requested but the Stan toolchain is not usable: ",
      stan_status$reason, ".\n",
      "Either install a working brms/rstan toolchain, or run with ",
      "MAIHDA_ENGINE=mle to use maximum likelihood instead. The MLE engine ",
      "produces the same outputs; see the script header for what differs."
    )
  }
  suppressPackageStartupMessages(library(brms))
  options(mc.cores = MCMC_CORES)
}

# Compiled Stan programs are cached by model structure and reused across
# outcomes. Compilation is a fixed cost of roughly a minute per structure, and
# on a slow processor paying it once per outcome rather than once per structure
# is a substantial and entirely avoidable waste.
compiled_model_cache <- new.env(parent = emptyenv())

# The tutorial leaves the intercept and random-effect priors at the brms
# defaults, which are weakly informative, and replaces the default flat prior
# on the fixed-effect coefficients with normal(0, 1). A flat prior on the logit
# scale is a poor choice: transformed to the probability scale it concentrates
# its mass near zero and one.
build_model_priors <- function(model_formula, data) {
  default_priors <- brms::get_prior(model_formula, data = data)
  has_fixed_effects <- any(default_priors$class == "b")
  if (!has_fixed_effects || !nzchar(PRIOR_FIXED_EFFECTS)) return(default_priors)
  default_priors$prior[default_priors$class == "b"] <- PRIOR_FIXED_EFFECTS
  default_priors
}

# Divergent transitions indicate the sampler could not explore the posterior
# geometry properly and that the estimates may be biased. They are counted
# rather than ignored.
count_divergent_transitions <- function(model) {
  sampler_parameters <- tryCatch(
    rstan::get_sampler_params(model$fit, inc_warmup = FALSE),
    error = function(condition) NULL
  )
  if (is.null(sampler_parameters)) return(NA_integer_)
  sum(vapply(sampler_parameters,
             function(chain) sum(chain[, "divergent__"]), numeric(1)))
}

brms_diagnostics <- function(model, label) {
  summary_table <- tryCatch(
    suppressWarnings(brms::posterior_summary(model)),
    error = function(condition) NULL
  )
  rhat_values <- tryCatch(suppressWarnings(brms::rhat(model)),
                          error = function(condition) NA_real_)
  ess_bulk <- tryCatch(
    suppressWarnings(summary(model)$fixed[, "Bulk_ESS"]),
    error = function(condition) NA_real_
  )
  ess_random <- tryCatch(
    suppressWarnings(summary(model)$random$stratum[, "Bulk_ESS"]),
    error = function(condition) NA_real_
  )
  list(
    max_rhat = suppressWarnings(max(rhat_values, na.rm = TRUE)),
    min_ess = suppressWarnings(min(c(ess_bulk, ess_random), na.rm = TRUE)),
    divergent_transitions = count_divergent_transitions(model),
    label = label
  )
}

# Fits a brms model, reusing a previously compiled Stan program for the same
# model structure where possible. If the design matrix has changed -- which
# happens when trimming removes a whole category for one outcome but not
# another -- the reuse fails and a fresh compile is done automatically.
fit_brms_model <- function(model_formula, data, cache_key, label) {
  cached <- if (REUSE_COMPILED_MODELS && exists(cache_key,
                                                envir = compiled_model_cache)) {
    get(cache_key, envir = compiled_model_cache)
  } else {
    NULL
  }

  fitted_model <- NULL
  if (!is.null(cached)) {
    fitted_model <- tryCatch(
      suppressMessages(suppressWarnings(
        stats::update(cached, newdata = data, recompile = FALSE,
                      chains = MCMC_CHAINS, iter = MCMC_ITERATIONS,
                      warmup = MCMC_WARMUP, cores = MCMC_CORES,
                      seed = MCMC_SEED, refresh = 0,
                      control = list(adapt_delta = MCMC_ADAPT_DELTA,
                                     max_treedepth = MCMC_MAX_TREEDEPTH))
      )),
      error = function(condition) {
        log_message("    Compiled model could not be reused for ", label,
                    " (", conditionMessage(condition),
                    "); recompiling.", level = "NOTE")
        NULL
      }
    )
  }

  if (is.null(fitted_model)) {
    priors <- build_model_priors(model_formula, data)
    fitted_model <- suppressMessages(suppressWarnings(
      brms::brm(
        formula = model_formula, data = data, prior = priors,
        chains = MCMC_CHAINS, iter = MCMC_ITERATIONS, warmup = MCMC_WARMUP,
        cores = MCMC_CORES, seed = MCMC_SEED, refresh = 0,
        control = list(adapt_delta = MCMC_ADAPT_DELTA,
                       max_treedepth = MCMC_MAX_TREEDEPTH)
      )
    ))
    if (REUSE_COMPILED_MODELS) {
      assign(cache_key, fitted_model, envir = compiled_model_cache)
    }
  }
  fitted_model
}

# Posterior draws of the between-stratum standard deviation. Everything the
# variance summaries need is derived from this one vector, which keeps the
# memory footprint small: the alternative, holding the full draws data frame,
# is orders of magnitude larger and is not needed.
stratum_sd_draws <- function(model) {
  draws <- tryCatch(
    as.matrix(model, variable = "sd_stratum__Intercept"),
    error = function(condition) NULL
  )
  if (is.null(draws) || ncol(draws) == 0L) return(numeric(0))
  as.numeric(draws[, 1L])
}

posterior_point <- function(x) {
  if (POSTERIOR_ROBUST) stats::median(x) else mean(x)
}

# VPC, MOR and the variance itself, each summarised from the same draws so the
# point estimates and intervals are mutually consistent.
bayesian_variance_summary <- function(sd_draws) {
  if (length(sd_draws) == 0L) {
    return(list(variance = NA_real_, variance_low = NA_real_,
                variance_high = NA_real_, vpc = NA_real_, vpc_low = NA_real_,
                vpc_high = NA_real_, mor = NA_real_, mor_low = NA_real_,
                mor_high = NA_real_))
  }
  variance_draws <- sd_draws^2
  vpc_draws <- 100 * logistic_vpc(variance_draws)
  mor_draws <- median_odds_ratio(variance_draws)
  quantiles <- function(x) unname(stats::quantile(x, c(0.025, 0.975)))
  list(
    variance = posterior_point(variance_draws),
    variance_low = quantiles(variance_draws)[1L],
    variance_high = quantiles(variance_draws)[2L],
    vpc = posterior_point(vpc_draws),
    vpc_low = quantiles(vpc_draws)[1L],
    vpc_high = quantiles(vpc_draws)[2L],
    mor = posterior_point(mor_draws),
    mor_low = quantiles(mor_draws)[1L],
    mor_high = quantiles(mor_draws)[2L]
  )
}

# The PCV is formed draw by draw from the two models, so its credible interval
# reflects uncertainty in both variances rather than being a point calculation
# on two summaries. Chains are matched by draw index, exactly as in the
# tutorial code.
bayesian_pcv_summary <- function(sd_draws_a, sd_draws_b) {
  n_draws <- min(length(sd_draws_a), length(sd_draws_b))
  if (n_draws == 0L) {
    return(list(pcv = NA_real_, pcv_low = NA_real_, pcv_high = NA_real_))
  }
  variance_a <- sd_draws_a[seq_len(n_draws)]^2
  variance_b <- sd_draws_b[seq_len(n_draws)]^2
  pcv_draws <- 100 * (variance_a - variance_b) / variance_a
  pcv_draws <- pcv_draws[is.finite(pcv_draws)]
  if (length(pcv_draws) == 0L) {
    return(list(pcv = NA_real_, pcv_low = NA_real_, pcv_high = NA_real_))
  }
  quantiles <- unname(stats::quantile(pcv_draws, c(0.025, 0.975)))
  list(pcv = posterior_point(pcv_draws), pcv_low = quantiles[1L],
       pcv_high = quantiles[2L])
}

# --- Fixed effects ----------------------------------------------------------
# Model terms are decomposed into variable and level, and the omitted reference
# categories are reinstated as explicit rows. A fixed-effects table that simply
# drops the reference level forces the reader to reconstruct the comparison
# being made, which is exactly the sort of avoidable friction that generates
# reviewer queries.
build_term_lookup <- function(strata, axes) {
  rows <- lapply(axes, function(axis) {
    levels_present <- levels(strata[[axis]])
    data.frame(
      term = paste0(axis, levels_present),
      variable = axis,
      variable_label = display_name_for(axis),
      level = levels_present,
      is_reference = seq_along(levels_present) == 1L
    )
  })
  bind_rows(
    data.frame(term = "(Intercept)", variable = "(Intercept)",
               variable_label = "Intercept", level = "",
               is_reference = FALSE),
    bind_rows(rows)
  )
}

extract_fixed_effects <- function(model, strata, axes, specification,
                                  specification_label, outcome,
                                  outcome_label) {
  coefficient_table <- as.data.frame(coef(summary(model)))
  coefficient_table$term <- rownames(coefficient_table)
  rownames(coefficient_table) <- NULL
  names(coefficient_table)[1:4] <- c("estimate_log_odds", "standard_error",
                                     "z_value", "p_value")

  lookup <- build_term_lookup(strata, axes)
  estimated <- coefficient_table %>%
    left_join(lookup, by = "term") %>%
    mutate(
      odds_ratio = exp(estimate_log_odds),
      confidence_low = exp(estimate_log_odds - 1.96 * standard_error),
      confidence_high = exp(estimate_log_odds + 1.96 * standard_error),
      is_reference = ifelse(is.na(is_reference), FALSE, is_reference)
    )

  # Reference rows carry an odds ratio of exactly one by construction.
  reference_rows <- lookup %>%
    filter(is_reference, !term %in% estimated$term) %>%
    mutate(
      estimate_log_odds = 0, standard_error = NA_real_, z_value = NA_real_,
      p_value = NA_real_, odds_ratio = 1, confidence_low = NA_real_,
      confidence_high = NA_real_
    )

  bind_rows(estimated, reference_rows) %>%
    mutate(
      specification = specification,
      specification_label = specification_label,
      outcome = outcome,
      outcome_label = outcome_label,
      # Terms that lme4 dropped for rank deficiency have no estimate and are
      # marked so, rather than appearing as a silent omission.
      term_status = case_when(
        is_reference ~ "reference",
        is.na(standard_error) ~ "not estimated",
        TRUE ~ "estimated"
      ),
      variable = ifelse(is.na(variable), term, variable),
      variable_label = ifelse(is.na(variable_label), term, variable_label),
      level = ifelse(is.na(level), "", level)
    ) %>%
    arrange(match(variable, c("(Intercept)", axes)), match(level, level)) %>%
    select(specification, specification_label, outcome, outcome_label,
           term, variable, variable_label, level, term_status, is_reference,
           estimate_log_odds, standard_error, z_value, p_value,
           odds_ratio, confidence_low, confidence_high)
}

# Convert a random-effect model into one row per stratum with additive-only and
# additive-plus-interaction predictions. Confidence intervals are approximate,
# matching the tutorial's warning that fixed and random-effect uncertainty are
# combined under a zero-covariance assumption.
make_stratum_predictions <- function(model_a, model_b, strata_data,
                                     axes, specification, specification_label,
                                     outcome, outcome_label) {
  # `formula(..., fixed.only = TRUE)` returns the fixed part directly and also
  # avoids relying on the former lme4::nobars location, which has moved to the
  # reformulas package in recent lme4 releases.
  fixed_formula_b <- formula(model_b, fixed.only = TRUE)
  fixed_formula_a <- formula(model_a, fixed.only = TRUE)

  matrix_b <- model.matrix(fixed_formula_b, data = strata_data)
  matrix_a <- model.matrix(fixed_formula_a, data = strata_data)

  beta_b <- fixef(model_b)
  beta_a <- fixef(model_a)

  # When lme4 drops rank-deficient columns, the coefficient vector is shorter
  # than the model matrix. Aligning explicitly prevents a dimension mismatch
  # that would otherwise surface as an opaque matrix multiplication error.
  if (!all(names(beta_b) %in% colnames(matrix_b))) {
    stop("Model B coefficients could not be aligned to the model matrix.")
  }
  matrix_b <- matrix_b[, names(beta_b), drop = FALSE]
  matrix_a <- matrix_a[, names(beta_a), drop = FALSE]

  fixed_eta_b <- as.numeric(matrix_b %*% beta_b)
  fixed_eta_a <- as.numeric(matrix_a %*% beta_a)

  covariance_b <- as.matrix(vcov(model_b))
  fixed_se_b <- sqrt(pmax(0, rowSums((matrix_b %*% covariance_b) * matrix_b)))

  # ranef() is computed once per model; the earlier version called it twice for
  # Model B, which doubled the cost on the largest specifications.
  ranef_b <- ranef(model_b, condVar = TRUE)$stratum
  ranef_a <- ranef(model_a, condVar = TRUE)$stratum
  posterior_variance_b <- attr(ranef_b, "postVar")

  random_b <- data.frame(
    level = rownames(ranef_b),
    random_effect_b = ranef_b[[1L]],
    random_se_b = sqrt(as.numeric(posterior_variance_b[1, 1, ]))
  )
  random_a <- data.frame(
    level = rownames(ranef_a),
    random_effect_a = ranef_a[[1L]]
  )

  aligned <- strata_data %>%
    mutate(level = as.character(stratum)) %>%
    left_join(random_b, by = "level") %>%
    left_join(random_a, by = "level")

  if (anyNA(aligned$random_effect_b) || anyNA(aligned$random_effect_a)) {
    stop("Random effects could not be aligned to every stratum.")
  }

  total_eta_b <- fixed_eta_b + aligned$random_effect_b
  total_se_b <- sqrt(fixed_se_b^2 + aligned$random_se_b^2)

  prediction <- aligned %>%
    mutate(
      specification = specification,
      specification_label = specification_label,
      outcome = outcome,
      outcome_label = outcome_label,
      observed_probability = events / n,
      null_fixed_probability = inv_logit(fixed_eta_a),
      null_total_probability = inv_logit(fixed_eta_a + random_effect_a),
      additive_probability = inv_logit(fixed_eta_b),
      total_probability = inv_logit(total_eta_b),
      total_probability_low = inv_logit(total_eta_b - 1.96 * total_se_b),
      total_probability_high = inv_logit(total_eta_b + 1.96 * total_se_b),
      interaction_log_odds = random_effect_b,
      interaction_log_odds_se = random_se_b,
      interaction_log_odds_low = random_effect_b - 1.96 * random_se_b,
      interaction_log_odds_high = random_effect_b + 1.96 * random_se_b,
      # This quantity is the absolute risk due to interaction (ARI) of Merlo,
      # Persmark and colleagues: the total predicted risk for the stratum minus
      # the risk predicted by the additive main effects alone. A positive value
      # means the stratum carries more risk than the simple addition of its
      # constituent positions implies. The literature name is carried in the
      # output alongside the descriptive one so the tables can be read directly
      # against published MAIHDA analyses.
      interaction_probability_difference =
        total_probability - additive_probability,
      interaction_probability_difference_low =
        inv_logit(fixed_eta_b + interaction_log_odds_low) - additive_probability,
      interaction_probability_difference_high =
        inv_logit(fixed_eta_b + interaction_log_odds_high) - additive_probability,
      absolute_risk = total_probability,
      absolute_risk_due_to_interaction = interaction_probability_difference,
      # The z statistic treats the conditional mode and its conditional
      # standard deviation as an approximate Wald quantity. This is the usual
      # applied convention and is what the interval on the figures shows.
      interaction_z = ifelse(random_se_b > 0, random_effect_b / random_se_b,
                             NA_real_),
      interaction_p_value = 2 * pnorm(-abs(interaction_z)),
      interaction_distinguishable =
        !is.na(interaction_log_odds_low) &
        (interaction_log_odds_low > 0 | interaction_log_odds_high < 0),
      prediction_rank = rank(total_probability, ties.method = "first"),
      interaction_rank = rank(interaction_probability_difference,
                              ties.method = "first")
    )

  # Benjamini-Hochberg control across the strata within this analysis. With
  # several hundred simultaneous comparisons the uncorrected flag alone will
  # identify strata by chance at a predictable rate, so both are carried
  # forward and the corrected one is used for the headline count.
  prediction$interaction_q_value <- if (all(is.na(prediction$interaction_p_value))) {
    NA_real_
  } else {
    p.adjust(prediction$interaction_p_value, method = "BH")
  }
  prediction$interaction_distinguishable_fdr <-
    !is.na(prediction$interaction_q_value) &
    prediction$interaction_q_value < FDR_LEVEL

  # --- Artefact diagnostics for the interaction residuals --------------------
  # An interaction that is "significant" can be an artefact of the probability
  # scale or of precision rather than a real departure from additivity. Three
  # checks travel with every stratum so a finding can be interrogated rather
  # than taken at face value.
  #
  # 1. Scale compression. Near a predicted probability of 0 or 1 the logistic
  #    curve is flat, so a given shift in log-odds maps to a small shift in
  #    probability, and conversely a stratum pinned near the ceiling cannot
  #    move upwards. A sizeable probability-scale ARI accompanied by a
  #    negligible log-odds residual is compression, not interaction.
  # 2. Precision-driven detection. Only a well-populated stratum has intervals
  #    tight enough to exclude zero, so significance tracks size. A stratum
  #    among the largest in the analysis warrants a check that the effect is
  #    substantively meaningful and not merely well measured.
  # 3. Marginal magnitude. An interval that only just excludes zero is far
  #    weaker evidence than one clear of it.
  prediction <- prediction %>%
    mutate(
      additive_probability_extreme =
        additive_probability > 0.80 | additive_probability < 0.05,
      # Ratio of the probability-scale effect to what the same log-odds shift
      # would produce at a predicted probability of 0.5, where the logistic
      # curve is steepest. Values well below 1 indicate compression.
      scale_compression_ratio = ifelse(
        abs(interaction_log_odds) > 1e-8,
        abs(interaction_probability_difference) /
          abs(inv_logit(interaction_log_odds) - 0.5),
        NA_real_
      ),
      possible_ceiling_artefact =
        interaction_distinguishable &
        additive_probability_extreme &
        !is.na(scale_compression_ratio) & scale_compression_ratio < 0.5,
      stratum_size_percentile = 100 * rank(n) / length(n),
      precision_driven = interaction_distinguishable &
        stratum_size_percentile >= 90,
      interval_excludes_zero_marginally = interaction_distinguishable &
        pmin(abs(interaction_log_odds_low), abs(interaction_log_odds_high)) <
          0.1 * abs(interaction_log_odds),
      artefact_flags = paste0(
        ifelse(possible_ceiling_artefact, "scale-compression;", ""),
        ifelse(precision_driven, "precision-driven;", ""),
        ifelse(interval_excludes_zero_marginally, "marginal-interval;", "")
      ),
      artefact_flags = ifelse(nzchar(artefact_flags),
                              sub(";$", "", artefact_flags), "")
    )

  prediction %>%
    select(
      specification, specification_label, outcome, outcome_label,
      all_of(axes), stratum, n, events, non_events,
      observed_probability, null_fixed_probability, null_total_probability,
      additive_probability, total_probability, total_probability_low,
      total_probability_high, interaction_log_odds, interaction_log_odds_se,
      interaction_log_odds_low, interaction_log_odds_high,
      interaction_probability_difference,
      interaction_probability_difference_low,
      interaction_probability_difference_high,
      absolute_risk, absolute_risk_due_to_interaction,
      additive_probability_extreme, scale_compression_ratio,
      possible_ceiling_artefact, stratum_size_percentile, precision_driven,
      interval_excludes_zero_marginally, artefact_flags,
      interaction_z, interaction_p_value, interaction_q_value,
      interaction_distinguishable, interaction_distinguishable_fdr,
      prediction_rank, interaction_rank
    )
}

# Stratum predictions from the posterior. Returns exactly the same columns as
# the maximum-likelihood version, so every downstream table, figure and
# manuscript sentence is indifferent to which engine produced them.
#
# The key difference is that the absolute risk and the absolute risk due to
# interaction are formed on the probability scale within each draw and only
# then summarised. This is what the tutorial's Bayesian code does, and it is
# why the Bayesian route needs no zero-covariance assumption: the joint
# uncertainty in the fixed and random parts is already carried by the draws.
make_stratum_predictions_bayesian <- function(model_a, model_b, strata_data,
                                              axes, specification,
                                              specification_label,
                                              outcome, outcome_label) {
  # Linear predictors on the log-odds scale, as draws x strata matrices.
  # posterior_linpred is used rather than posterior_epred because for a
  # binomial model with trials(n) the expectation is a count, not a
  # probability.
  linear_total_b <- brms::posterior_linpred(model_b, newdata = strata_data,
                                            re_formula = NULL)
  linear_fixed_b <- brms::posterior_linpred(model_b, newdata = strata_data,
                                            re_formula = NA)
  linear_total_a <- brms::posterior_linpred(model_a, newdata = strata_data,
                                            re_formula = NULL)
  linear_fixed_a <- brms::posterior_linpred(model_a, newdata = strata_data,
                                            re_formula = NA)

  probability_total <- inv_logit(linear_total_b)
  probability_fixed <- inv_logit(linear_fixed_b)
  # The absolute risk due to interaction, per draw.
  ari_draws <- probability_total - probability_fixed
  # The interaction residual on the log-odds scale, per draw.
  residual_draws <- linear_total_b - linear_fixed_b

  summarise_columns <- function(draws, point = TRUE) {
    centre <- if (point) {
      if (POSTERIOR_ROBUST) apply(draws, 2L, stats::median) else colMeans(draws)
    } else NULL
    lower <- apply(draws, 2L, stats::quantile, probs = 0.025)
    upper <- apply(draws, 2L, stats::quantile, probs = 0.975)
    list(centre = centre, lower = unname(lower), upper = unname(upper))
  }

  total_summary <- summarise_columns(probability_total)
  additive_summary <- summarise_columns(probability_fixed)
  ari_summary <- summarise_columns(ari_draws)
  residual_summary <- summarise_columns(residual_draws)
  null_total_summary <- summarise_columns(inv_logit(linear_total_a))
  null_fixed_summary <- summarise_columns(inv_logit(linear_fixed_a))

  # The Bayesian analogue of a two-sided significance test is the probability
  # of direction: the posterior probability that the effect has the sign of its
  # point estimate. It is reported directly, and is additionally mapped to a
  # two-sided pseudo p-value so that the same false discovery rate machinery
  # can be applied as under maximum likelihood, keeping the two engines
  # comparable.
  probability_of_direction <- apply(residual_draws, 2L, function(column) {
    max(mean(column > 0), mean(column < 0))
  })
  pseudo_p_value <- 2 * (1 - probability_of_direction)

  # Free the large matrices before assembling the output frame.
  rm(linear_total_b, linear_fixed_b, linear_total_a, linear_fixed_a,
     probability_total, probability_fixed, ari_draws, residual_draws)
  invisible(gc(verbose = FALSE))

  prediction <- strata_data %>%
    mutate(
      specification = specification,
      specification_label = specification_label,
      outcome = outcome,
      outcome_label = outcome_label,
      observed_probability = events / n,
      null_fixed_probability = null_fixed_summary$centre,
      null_total_probability = null_total_summary$centre,
      additive_probability = additive_summary$centre,
      total_probability = total_summary$centre,
      total_probability_low = total_summary$lower,
      total_probability_high = total_summary$upper,
      interaction_log_odds = residual_summary$centre,
      interaction_log_odds_se = (residual_summary$upper -
                                   residual_summary$lower) / (2 * 1.96),
      interaction_log_odds_low = residual_summary$lower,
      interaction_log_odds_high = residual_summary$upper,
      interaction_probability_difference = ari_summary$centre,
      interaction_probability_difference_low = ari_summary$lower,
      interaction_probability_difference_high = ari_summary$upper,
      absolute_risk = total_probability,
      absolute_risk_due_to_interaction = interaction_probability_difference,
      interaction_probability_of_direction = probability_of_direction,
      interaction_z = NA_real_,
      interaction_p_value = pseudo_p_value,
      # The canonical criterion in the tutorial: the 95% credible interval for
      # the stratum random effect excludes zero.
      interaction_distinguishable =
        interaction_log_odds_low > 0 | interaction_log_odds_high < 0,
      prediction_rank = rank(total_probability, ties.method = "first"),
      interaction_rank = rank(interaction_probability_difference,
                              ties.method = "first")
    )

  prediction$interaction_q_value <- p.adjust(prediction$interaction_p_value,
                                             method = "BH")
  prediction$interaction_distinguishable_fdr <-
    !is.na(prediction$interaction_q_value) &
    prediction$interaction_q_value < FDR_LEVEL

  # --- Artefact diagnostics for the interaction residuals --------------------
  # An interaction that is "significant" can be an artefact of the probability
  # scale or of precision rather than a real departure from additivity. Three
  # checks travel with every stratum so a finding can be interrogated rather
  # than taken at face value.
  #
  # 1. Scale compression. Near a predicted probability of 0 or 1 the logistic
  #    curve is flat, so a given shift in log-odds maps to a small shift in
  #    probability, and conversely a stratum pinned near the ceiling cannot
  #    move upwards. A sizeable probability-scale ARI accompanied by a
  #    negligible log-odds residual is compression, not interaction.
  # 2. Precision-driven detection. Only a well-populated stratum has intervals
  #    tight enough to exclude zero, so significance tracks size. A stratum
  #    among the largest in the analysis warrants a check that the effect is
  #    substantively meaningful and not merely well measured.
  # 3. Marginal magnitude. An interval that only just excludes zero is far
  #    weaker evidence than one clear of it.
  prediction <- prediction %>%
    mutate(
      additive_probability_extreme =
        additive_probability > 0.80 | additive_probability < 0.05,
      # Ratio of the probability-scale effect to what the same log-odds shift
      # would produce at a predicted probability of 0.5, where the logistic
      # curve is steepest. Values well below 1 indicate compression.
      scale_compression_ratio = ifelse(
        abs(interaction_log_odds) > 1e-8,
        abs(interaction_probability_difference) /
          abs(inv_logit(interaction_log_odds) - 0.5),
        NA_real_
      ),
      possible_ceiling_artefact =
        interaction_distinguishable &
        additive_probability_extreme &
        !is.na(scale_compression_ratio) & scale_compression_ratio < 0.5,
      stratum_size_percentile = 100 * rank(n) / length(n),
      precision_driven = interaction_distinguishable &
        stratum_size_percentile >= 90,
      interval_excludes_zero_marginally = interaction_distinguishable &
        pmin(abs(interaction_log_odds_low), abs(interaction_log_odds_high)) <
          0.1 * abs(interaction_log_odds),
      artefact_flags = paste0(
        ifelse(possible_ceiling_artefact, "scale-compression;", ""),
        ifelse(precision_driven, "precision-driven;", ""),
        ifelse(interval_excludes_zero_marginally, "marginal-interval;", "")
      ),
      artefact_flags = ifelse(nzchar(artefact_flags),
                              sub(";$", "", artefact_flags), "")
    )

  prediction %>%
    select(
      specification, specification_label, outcome, outcome_label,
      all_of(axes), stratum, n, events, non_events,
      observed_probability, null_fixed_probability, null_total_probability,
      additive_probability, total_probability, total_probability_low,
      total_probability_high, interaction_log_odds, interaction_log_odds_se,
      interaction_log_odds_low, interaction_log_odds_high,
      interaction_probability_difference,
      interaction_probability_difference_low,
      interaction_probability_difference_high,
      absolute_risk, absolute_risk_due_to_interaction,
      additive_probability_extreme, scale_compression_ratio,
      possible_ceiling_artefact, stratum_size_percentile, precision_driven,
      interval_excludes_zero_marginally, artefact_flags,
      interaction_probability_of_direction,
      interaction_z, interaction_p_value, interaction_q_value,
      interaction_distinguishable, interaction_distinguishable_fdr,
      prediction_rank, interaction_rank
    )
}

# Fixed effects from a brms fit, presented in the same shape as the maximum
# likelihood version. Credible intervals replace confidence intervals and the
# p-value column is left empty, since a posterior has no such quantity.
extract_fixed_effects_bayesian <- function(model, strata, axes, specification,
                                           specification_label, outcome,
                                           outcome_label) {
  summary_table <- as.data.frame(brms::fixef(model, robust = POSTERIOR_ROBUST,
                                             probs = c(0.025, 0.975)))
  summary_table$term <- rownames(summary_table)
  rownames(summary_table) <- NULL
  names(summary_table)[1:4] <- c("estimate_log_odds", "standard_error",
                                 "confidence_low_logit", "confidence_high_logit")

  # brms strips non-syntactic characters from term names; align them back to
  # the model matrix naming so the axis and level lookup matches.
  lookup <- build_term_lookup(strata, axes)
  lookup$term_clean <- make.names(lookup$term)
  summary_table$term_clean <- make.names(summary_table$term)

  estimated <- summary_table %>%
    left_join(lookup %>% select(-term), by = "term_clean") %>%
    mutate(
      odds_ratio = exp(estimate_log_odds),
      confidence_low = exp(confidence_low_logit),
      confidence_high = exp(confidence_high_logit),
      p_value = NA_real_,
      z_value = NA_real_,
      is_reference = ifelse(is.na(is_reference), FALSE, is_reference)
    )

  reference_rows <- lookup %>%
    filter(is_reference, !term_clean %in% estimated$term_clean) %>%
    mutate(
      estimate_log_odds = 0, standard_error = NA_real_, z_value = NA_real_,
      p_value = NA_real_, odds_ratio = 1, confidence_low = NA_real_,
      confidence_high = NA_real_
    ) %>%
    rename(term_original = term)

  bind_rows(
    estimated,
    reference_rows %>% mutate(term = term_original) %>% select(-term_original)
  ) %>%
    mutate(
      specification = specification,
      specification_label = specification_label,
      outcome = outcome,
      outcome_label = outcome_label,
      term = ifelse(is.na(term), term_clean, term),
      term_status = case_when(
        is_reference ~ "reference",
        is.na(standard_error) ~ "not estimated",
        TRUE ~ "estimated"
      ),
      variable = ifelse(is.na(variable), term, variable),
      variable_label = ifelse(is.na(variable_label), term, variable_label),
      level = ifelse(is.na(level), "", level)
    ) %>%
    arrange(match(variable, c("(Intercept)", axes)), match(level, level)) %>%
    select(specification, specification_label, outcome, outcome_label,
           term, variable, variable_label, level, term_status, is_reference,
           estimate_log_odds, standard_error, z_value, p_value,
           odds_ratio, confidence_low, confidence_high)
}

# =============================================================================
# 5. DATA PREPARATION AND THE MINIMUM CELL SIZE RULE
# =============================================================================

prepare_analysis_data <- function(data, axes, references, outcomes) {
  prepared <- data[, c(axes, outcomes), drop = FALSE]

  for (axis in axes) {
    # factor() also strips any ordered-factor class. This matters: an ordered
    # factor would otherwise be given polynomial contrasts by model.matrix,
    # producing .L and .Q terms in place of the interpretable level contrasts
    # the analysis is built around.
    values <- prepared[[axis]]
    if (is.ordered(values)) {
      values <- factor(as.character(values), levels = levels(values))
    }
    prepared[[axis]] <- droplevels(factor(values))

    reference <- unname(references[axis])
    if (!reference %in% levels(prepared[[axis]])) {
      stop("Reference level '", reference, "' is absent from ", axis,
           ". Observed levels: ",
           paste(levels(prepared[[axis]]), collapse = ", "))
    }
    prepared[[axis]] <- relevel(prepared[[axis]], ref = reference)
  }
  prepared
}

make_stratum_counts <- function(data, axes, outcome) {
  complete <- complete.cases(data[, c(axes, outcome), drop = FALSE])
  excluded <- sum(!complete)
  model_data <- droplevels(data[complete, c(axes, outcome), drop = FALSE])

  if (nrow(model_data) == 0L) {
    stop("No complete records remain for ", outcome, ".")
  }

  # interaction() creates a stable factor identifier while retaining each axis
  # separately for the additive main-effects model and output tables.
  model_data$stratum <- interaction(
    model_data[, axes, drop = FALSE], drop = TRUE, sep = " | "
  )

  strata <- model_data %>%
    group_by(across(all_of(c(axes, "stratum")))) %>%
    summarise(
      n = n(),
      events = sum(.data[[outcome]] == 1),
      non_events = sum(.data[[outcome]] == 0),
      .groups = "drop"
    ) %>%
    mutate(stratum = factor(stratum))

  # An internal consistency check on the collapse. If this ever fails the
  # binomial denominators are wrong and every downstream quantity is invalid.
  if (!isTRUE(all.equal(strata$n, strata$events + strata$non_events))) {
    stop("Stratum totals do not reconcile for ", outcome,
         "; events plus non-events does not equal n.")
    }

  list(
    data = strata,
    included_n = nrow(model_data),
    excluded_n = excluded,
    observed_strata = nrow(strata),
    possible_strata = prod(vapply(model_data[axes], nlevels, integer(1)))
  )
}

# --- The minimum cell size rule ---------------------------------------------
# Applied after the binomial collapse and before fitting. Everything removed is
# counted, so the analytic sample can be reconciled against the full cohort at
# any point.
apply_minimum_cell_rule <- function(strata, axes, references,
                                    min_n, min_events) {
  retained_flag <- strata$n >= min_n & strata$events >= min_events
  retained <- strata[retained_flag, , drop = FALSE]
  excluded <- strata[!retained_flag, , drop = FALSE]

  accounting <- data.frame(
    min_stratum_n = min_n,
    min_stratum_events = min_events,
    strata_before = nrow(strata),
    strata_retained = nrow(retained),
    strata_excluded = nrow(excluded),
    strata_retained_percent = if (nrow(strata) > 0) {
      100 * nrow(retained) / nrow(strata)
    } else NA_real_,
    individuals_before = sum(strata$n),
    individuals_retained = sum(retained$n),
    individuals_excluded = sum(excluded$n),
    individuals_retained_percent = if (sum(strata$n) > 0) {
      100 * sum(retained$n) / sum(strata$n)
    } else NA_real_,
    events_before = sum(strata$events),
    events_retained = sum(retained$events),
    events_excluded = sum(excluded$events),
    events_retained_percent = if (sum(strata$events) > 0) {
      100 * sum(retained$events) / sum(strata$events)
    } else NA_real_,
    excluded_median_n = if (nrow(excluded) > 0) median(excluded$n) else NA_real_,
    excluded_max_n = if (nrow(excluded) > 0) max(excluded$n) else NA_real_,
    retained_zero_event_strata = sum(retained$events == 0),
    retained_zero_non_event_strata = sum(retained$non_events == 0)
  )

  if (nrow(retained) == 0L) {
    return(list(data = retained, excluded = excluded, accounting = accounting,
                model_axes = character(0), reference_changes = data.frame(),
                dropped_axes = axes))
  }

  # Trimming can remove an entire category of an axis. Two consequences have to
  # be handled explicitly rather than left to surface as a modelling error.
  retained <- droplevels(retained)
  retained$stratum <- droplevels(factor(retained$stratum))

  reference_changes <- list()
  dropped_axes <- character(0)
  model_axes <- character(0)

  for (axis in axes) {
    levels_present <- levels(retained[[axis]])

    # An axis reduced to a single category carries no information once the
    # intercept is in the model. It is kept for labelling but removed from the
    # fixed-effects formula, where it would be aliased with the intercept.
    if (length(levels_present) < 2L) {
      dropped_axes <- c(dropped_axes, axis)
      next
    }
    model_axes <- c(model_axes, axis)

    intended_reference <- unname(references[axis])
    if (!intended_reference %in% levels_present) {
      # Falling back to the largest surviving category keeps the contrasts
      # well determined. This is recorded prominently because it changes the
      # interpretation of every odds ratio for that axis.
      level_sizes <- tapply(retained$n, retained[[axis]], sum)
      replacement <- names(which.max(level_sizes))
      reference_changes[[length(reference_changes) + 1L]] <- data.frame(
        variable = axis, intended_reference = intended_reference,
        used_reference = replacement
      )
      retained[[axis]] <- relevel(retained[[axis]], ref = replacement)
    } else {
      retained[[axis]] <- relevel(retained[[axis]], ref = intended_reference)
    }
  }

  list(
    data = retained,
    excluded = excluded,
    accounting = accounting,
    model_axes = model_axes,
    reference_changes = bind_rows(reference_changes),
    dropped_axes = dropped_axes
  )
}

stratum_size_summary <- function(strata, specification, specification_label,
                                 outcome, outcome_label,
                                 included_n, excluded_n, possible_strata,
                                 stage) {
  sizes <- strata$n
  if (length(sizes) == 0L) sizes <- NA_real_
  data.frame(
    specification = specification,
    specification_label = specification_label,
    outcome = outcome,
    outcome_label = outcome_label,
    stage = stage,
    included_n = included_n,
    excluded_n = excluded_n,
    possible_strata = possible_strata,
    observed_strata = nrow(strata),
    empty_strata = possible_strata - nrow(strata),
    minimum_n = min(sizes, na.rm = TRUE),
    q1_n = unname(quantile(sizes, 0.25, na.rm = TRUE)),
    median_n = median(sizes, na.rm = TRUE),
    mean_n = mean(sizes, na.rm = TRUE),
    q3_n = unname(quantile(sizes, 0.75, na.rm = TRUE)),
    maximum_n = max(sizes, na.rm = TRUE),
    strata_n_100_plus = sum(strata$n >= 100),
    strata_n_50_plus = sum(strata$n >= 50),
    strata_n_30_plus = sum(strata$n >= 30),
    strata_n_20_plus = sum(strata$n >= 20),
    strata_n_10_plus = sum(strata$n >= 10),
    strata_n_below_10 = sum(strata$n < 10),
    strata_no_events = sum(strata$events == 0),
    strata_no_non_events = sum(strata$non_events == 0)
  )
}

# =============================================================================
# 6. MODEL FITTING
# =============================================================================

build_model_formulae <- function(model_axes) {
  formula_a <- as.formula("cbind(events, non_events) ~ 1 + (1 | stratum)")
  formula_b <- if (length(model_axes) == 0L) {
    formula_a
  } else {
    as.formula(paste(
      "cbind(events, non_events) ~",
      paste(model_axes, collapse = " + "),
      "+ (1 | stratum)"
    ))
  }
  list(a = formula_a, b = formula_b)
}

# --- Partially adjusted single-axis models -----------------------------------
# One model per axis, each adding a single axis to the null model. The
# proportional change in variance from Model A gives that axis's individual
# contribution to the between-stratum variance. See the caveat in section 1:
# this is a secondary, descriptive output and invites single-axis reading if
# given more weight than the collective additive effect.
#
# The contributions do not sum to the Model B PCV, and are not meant to: the
# axes are correlated in the population, so their separate contributions
# overlap. Reporting them alongside the joint Model B figure is what makes the
# overlap visible, and the sum is deliberately not presented as a total.
fit_axis_decomposition <- function(strata, model_axes, variance_a,
                                   specification, specification_label,
                                   outcome, outcome_label) {
  if (!RUN_AXIS_DECOMPOSITION || length(model_axes) == 0L) return(NULL)

  bind_rows(lapply(model_axes, function(axis) {
    axis_formula <- as.formula(paste(
      "cbind(events, non_events) ~", axis, "+ (1 | stratum)"
    ))
    attempt <- tryCatch(
      fit_glmer_robust(axis_formula, strata, paste0("single-axis: ", axis)),
      error = function(condition) NULL
    )
    if (is.null(attempt)) {
      return(data.frame(
        specification = specification, specification_label = specification_label,
        outcome = outcome, outcome_label = outcome_label,
        axis = axis, axis_label = display_name_for(axis),
        variance_model_2 = NA_real_, vpc_model_2_percent = NA_real_,
        pcv_percent = NA_real_, singular = NA, status = "fit failed"
      ))
    }
    variance_2 <- random_intercept_variance(attempt$model)
    data.frame(
      specification = specification, specification_label = specification_label,
      outcome = outcome, outcome_label = outcome_label,
      axis = axis, axis_label = display_name_for(axis),
      variance_model_2 = variance_2,
      vpc_model_2_percent = 100 * logistic_vpc(variance_2),
      pcv_percent = if (is.finite(variance_a) && variance_a > 0) {
        100 * (variance_a - variance_2) / variance_a
      } else NA_real_,
      singular = attempt$singular,
      status = "fitted"
    )
  }))
}

fit_one_maihda <- function(strata, axes, model_axes, specification,
                           specification_label, outcome, outcome_label,
                           included_n, excluded_n, possible_strata,
                           accounting, eligibility_column = NA_character_,
                           observable_n = NA_integer_) {
  log_message("  Fitting ", specification, " / ", outcome,
              " (", nrow(strata), " strata, ", sum(strata$n), " individuals)")

  formulae <- build_model_formulae(model_axes)
  model_start_time <- Sys.time()

  if (ESTIMATION_ENGINE == "bayesian") {
    bayes_formula_a <- brms::brmsformula(
      events | trials(n) ~ 1 + (1 | stratum), family = binomial("logit")
    )
    bayes_formula_b <- if (length(model_axes) == 0L) {
      bayes_formula_a
    } else {
      brms::brmsformula(
        as.formula(paste("events | trials(n) ~",
                         paste(model_axes, collapse = " + "),
                         "+ (1 | stratum)")),
        family = binomial("logit")
      )
    }

    model_a <- fit_brms_model(bayes_formula_a, strata, "null_model",
                              paste0(outcome, " / Model A"))
    diagnostics_a <- brms_diagnostics(model_a, "A")
    sd_draws_a <- stratum_sd_draws(model_a)

    model_b <- fit_brms_model(bayes_formula_b, strata, "additive_model",
                              paste0(outcome, " / Model B"))
    diagnostics_b <- brms_diagnostics(model_b, "B")
    sd_draws_b <- stratum_sd_draws(model_b)

    summary_a <- bayesian_variance_summary(sd_draws_a)
    summary_b <- bayesian_variance_summary(sd_draws_b)
    pcv_summary <- bayesian_pcv_summary(sd_draws_a, sd_draws_b)

    variance_a <- summary_a$variance
    variance_b <- summary_b$variance
    pcv <- pcv_summary$pcv

    interval_a <- list(
      variance_low = summary_a$variance_low,
      variance_high = summary_a$variance_high,
      vpc_low = summary_a$vpc_low, vpc_high = summary_a$vpc_high,
      mor_low = summary_a$mor_low, mor_high = summary_a$mor_high,
      interval_method = "posterior credible interval"
    )
    interval_b <- list(
      variance_low = summary_b$variance_low,
      variance_high = summary_b$variance_high,
      vpc_low = summary_b$vpc_low, vpc_high = summary_b$vpc_high,
      mor_low = summary_b$mor_low, mor_high = summary_b$mor_high,
      interval_method = "posterior credible interval"
    )

    # Convergence is reported rather than assumed. A model that fails these
    # checks still produces output, but the failure is recorded in the metrics
    # and logged so it cannot pass unnoticed.
    for (diagnostics in list(diagnostics_a, diagnostics_b)) {
      if (!is.na(diagnostics$max_rhat) && diagnostics$max_rhat > MAX_RHAT) {
        log_message("    Model ", diagnostics$label, " R-hat is ",
                    formatC(diagnostics$max_rhat, format = "f", digits = 3),
                    ", above the ", MAX_RHAT, " threshold. Consider raising ",
                    "MAIHDA_MCMC_ITER.", level = "WARN")
      }
      if (!is.na(diagnostics$min_ess) &&
          diagnostics$min_ess < MIN_EFFECTIVE_SAMPLE_SIZE) {
        log_message("    Model ", diagnostics$label,
                    " minimum effective sample size is ",
                    round(diagnostics$min_ess), ", below the ",
                    MIN_EFFECTIVE_SAMPLE_SIZE, " threshold. Consider raising ",
                    "MAIHDA_MCMC_ITER.", level = "WARN")
      }
      if (!is.na(diagnostics$divergent_transitions) &&
          diagnostics$divergent_transitions > 0) {
        log_message("    Model ", diagnostics$label, " had ",
                    diagnostics$divergent_transitions,
                    " divergent transitions. Consider raising ",
                    "MAIHDA_MCMC_ADAPT_DELTA.", level = "WARN")
      }
    }
    rm(sd_draws_a, sd_draws_b)
  } else {
    fit_a <- fit_glmer_robust(formulae$a, strata, paste0(outcome, " / Model A"))
    fit_b <- fit_glmer_robust(formulae$b, strata, paste0(outcome, " / Model B"))
    model_a <- fit_a$model
    model_b <- fit_b$model

    variance_a <- random_intercept_variance(model_a)
    variance_b <- random_intercept_variance(model_b)
    pcv <- if (is.finite(variance_a) && variance_a > 0) {
      100 * (variance_a - variance_b) / variance_a
    } else {
      NA_real_
    }
    interval_a <- variance_interval(model_a)
    interval_b <- variance_interval(model_b)
    diagnostics_a <- list(max_rhat = NA_real_, min_ess = NA_real_,
                          divergent_transitions = NA_integer_)
    diagnostics_b <- diagnostics_a
    pcv_summary <- list(pcv = pcv, pcv_low = NA_real_, pcv_high = NA_real_)
    summary_a <- list(vpc = 100 * logistic_vpc(variance_a),
                      mor = median_odds_ratio(variance_a))
    summary_b <- list(vpc = 100 * logistic_vpc(variance_b),
                      mor = median_odds_ratio(variance_b))
  }

  # Under MCMC there is no boundary solution to report: the variance is a
  # posterior distribution bounded below by zero. The equivalent finding is a
  # credible interval whose lower limit sits essentially at zero, which is
  # reported here in the same interpretive terms.
  if (ESTIMATION_ENGINE == "bayesian") {
    if (!is.na(interval_b$vpc_high) && interval_b$vpc_high < 1) {
      log_message("    Model B between-stratum variance is close to zero ",
                  "(VPC 95% CrI up to ",
                  formatC(interval_b$vpc_high, format = "f", digits = 2),
                  "%). The additive main effects account for essentially all ",
                  "between-stratum variation, so little residual ",
                  "intersectional interaction is detectable.", level = "NOTE")
    }
    if (!is.na(interval_a$vpc_high) && interval_a$vpc_high < 1) {
      log_message("    Model A between-stratum variance is close to zero. ",
                  "This implies almost no detectable between-stratum ",
                  "variation and should be investigated before the results ",
                  "are used.", level = "WARN")
    }
  } else if (fit_a$singular) {
    log_message("    Model A variance is at the zero boundary. This implies no ",
                "detectable between-stratum variation and should be ",
                "investigated before the results are used.", level = "WARN")
  }
  if (ESTIMATION_ENGINE != "bayesian" && fit_b$singular) {
    log_message("    Model B variance is at the zero boundary",
                if (fit_b$boundary_confirmed) {
                  " (confirmed across all optimisers)"
                } else {
                  " (optimisers disagreed; see the fit attempts table)"
                },
                ". The additive main effects account for essentially all ",
                "between-stratum variation, so no residual intersectional ",
                "interaction is detected. PCV is at or near 100%.",
                level = "NOTE")
  }

  axis_decomposition <- fit_axis_decomposition(
    strata, model_axes, variance_a, specification, specification_label,
    outcome, outcome_label
  )
  if (!is.null(axis_decomposition)) {
    ranked_axes <- axis_decomposition[
      order(-axis_decomposition$pcv_percent), , drop = FALSE
    ]
    leading <- ranked_axes[1L, ]
    log_message("    Largest single-axis contribution: ", leading$axis_label,
                " (PCV ", formatC(leading$pcv_percent, format = "f", digits = 1),
                "%)")
  }

  predictions <- if (ESTIMATION_ENGINE == "bayesian") {
    make_stratum_predictions_bayesian(
      model_a, model_b, strata, axes, specification, specification_label,
      outcome, outcome_label
    )
  } else {
    make_stratum_predictions(
      model_a, model_b, strata, axes, specification, specification_label,
      outcome, outcome_label
    )
  }

  # AUC is computed from the prediction frame itself so that it cannot be
  # invalidated by a change in row ordering elsewhere.
  auc_for <- function(column) {
    weighted_auc(predictions$events, predictions$non_events,
                 predictions[[column]])
  }

  metrics <- data.frame(
    specification = specification,
    specification_label = specification_label,
    outcome = outcome,
    outcome_label = outcome_label,
    min_stratum_n = MIN_STRATUM_N,
    min_stratum_events = MIN_STRATUM_EVENTS,
    eligibility_column = eligibility_column,
    observable_n = observable_n,
    included_n = included_n,
    excluded_n = excluded_n,
    analysed_n = sum(strata$n),
    possible_strata = possible_strata,
    strata_before_rule = accounting$strata_before,
    observed_strata = nrow(strata),
    strata_excluded_by_rule = accounting$strata_excluded,
    individuals_excluded_by_rule = accounting$individuals_excluded,
    individuals_retained_percent = accounting$individuals_retained_percent,
    events_retained_percent = accounting$events_retained_percent,
    event_n = sum(strata$events),
    event_percent = 100 * sum(strata$events) / sum(strata$n),
    variance_model_a = variance_a,
    variance_model_a_low = interval_a$variance_low,
    variance_model_a_high = interval_a$variance_high,
    variance_model_b = variance_b,
    variance_model_b_low = interval_b$variance_low,
    variance_model_b_high = interval_b$variance_high,
    vpc_model_a_percent = 100 * logistic_vpc(variance_a),
    vpc_model_a_low = interval_a$vpc_low,
    vpc_model_a_high = interval_a$vpc_high,
    vpc_model_b_percent = 100 * logistic_vpc(variance_b),
    vpc_model_b_low = interval_b$vpc_low,
    vpc_model_b_high = interval_b$vpc_high,
    pcv_percent = pcv,
    pcv_low = pcv_summary$pcv_low,
    pcv_high = pcv_summary$pcv_high,
    mor_model_a = median_odds_ratio(variance_a),
    mor_model_a_low = interval_a$mor_low,
    mor_model_a_high = interval_a$mor_high,
    mor_model_b = median_odds_ratio(variance_b),
    interval_method = interval_a$interval_method,
    auc_model_a_total = auc_for("null_total_probability"),
    auc_model_a_fixed = auc_for("null_fixed_probability"),
    auc_model_b_total = auc_for("total_probability"),
    auc_model_b_fixed = auc_for("additive_probability"),
    log_likelihood_model_a = if (ESTIMATION_ENGINE == "bayesian") NA_real_ else as.numeric(logLik(model_a)),
    log_likelihood_model_b = if (ESTIMATION_ENGINE == "bayesian") NA_real_ else as.numeric(logLik(model_b)),
    aic_model_a = if (ESTIMATION_ENGINE == "bayesian") NA_real_ else AIC(model_a),
    aic_model_b = if (ESTIMATION_ENGINE == "bayesian") NA_real_ else AIC(model_b),
    n_distinguishable = sum(predictions$interaction_distinguishable),
    n_distinguishable_fdr = sum(predictions$interaction_distinguishable_fdr),
    # The ratio of the extremes is how the applied MAIHDA literature conveys
    # the width of the intersectional gradient ("the highest risk was ten times
    # the lowest"). The absolute difference alone understates a gradient in a
    # rare outcome and overstates one in a common outcome.
    lowest_absolute_risk = min(predictions$total_probability),
    highest_absolute_risk = max(predictions$total_probability),
    absolute_risk_difference = max(predictions$total_probability) -
      min(predictions$total_probability),
    absolute_risk_ratio = if (min(predictions$total_probability) > 0) {
      max(predictions$total_probability) / min(predictions$total_probability)
    } else NA_real_,
    estimation_engine = ESTIMATION_ENGINE,
    optimizer_model_a = if (ESTIMATION_ENGINE == "bayesian") "stan/hmc" else fit_a$optimizer,
    optimizer_model_b = if (ESTIMATION_ENGINE == "bayesian") "stan/hmc" else fit_b$optimizer,
    # Under MCMC there is no boundary solution to detect: a variance is a
    # posterior distribution bounded below by zero, not a point estimate that
    # can sit exactly on the bound. The singular flags are therefore FALSE and
    # the near-zero case is conveyed by the credible interval instead.
    singular_model_a = if (ESTIMATION_ENGINE == "bayesian") FALSE else fit_a$singular,
    singular_model_b = if (ESTIMATION_ENGINE == "bayesian") FALSE else fit_b$singular,
    boundary_confirmed_model_a = if (ESTIMATION_ENGINE == "bayesian") NA else fit_a$boundary_confirmed,
    boundary_confirmed_model_b = if (ESTIMATION_ENGINE == "bayesian") NA else fit_b$boundary_confirmed,
    convergence_model_a = if (ESTIMATION_ENGINE == "bayesian") "MCMC" else extract_convergence_message(model_a),
    convergence_model_b = if (ESTIMATION_ENGINE == "bayesian") "MCMC" else extract_convergence_message(model_b),
    max_rhat_model_a = diagnostics_a$max_rhat,
    max_rhat_model_b = diagnostics_b$max_rhat,
    min_ess_model_a = diagnostics_a$min_ess,
    min_ess_model_b = diagnostics_b$min_ess,
    divergent_model_a = diagnostics_a$divergent_transitions,
    divergent_model_b = diagnostics_b$divergent_transitions,
    fit_minutes = as.numeric(difftime(Sys.time(), model_start_time, units = "mins")),
    model_axes = paste(model_axes, collapse = " + "),
    axes_dropped = paste(setdiff(axes, model_axes), collapse = ", ")
  )

  fit_attempts <- if (ESTIMATION_ENGINE == "bayesian") {
    bind_rows(
      data.frame(model = "A", optimizer = "stan/hmc", status = "fitted",
                 log_likelihood = NA_real_, variance = variance_a,
                 singular = NA, detail = paste0(
                   "rhat ", formatC(diagnostics_a$max_rhat, format = "f", digits = 3),
                   "; min ESS ", round(diagnostics_a$min_ess),
                   "; divergences ", diagnostics_a$divergent_transitions)),
      data.frame(model = "B", optimizer = "stan/hmc", status = "fitted",
                 log_likelihood = NA_real_, variance = variance_b,
                 singular = NA, detail = paste0(
                   "rhat ", formatC(diagnostics_b$max_rhat, format = "f", digits = 3),
                   "; min ESS ", round(diagnostics_b$min_ess),
                   "; divergences ", diagnostics_b$divergent_transitions))
    )
  } else {
    bind_rows(
      fit_a$attempts %>% mutate(model = "A"),
      fit_b$attempts %>% mutate(model = "B")
    )
  } %>%
    mutate(specification = specification, outcome = outcome) %>%
    select(specification, outcome, model, everything())

  fixed_effects <- if (ESTIMATION_ENGINE == "bayesian") {
    extract_fixed_effects_bayesian(
      model_b, strata, model_axes, specification, specification_label,
      outcome, outcome_label
    )
  } else {
    extract_fixed_effects(
      model_b, strata, model_axes, specification, specification_label,
      outcome, outcome_label
    )
  }

  if (SAVE_MODELS) {
    model_file <- file.path(
      OUTPUT_ROOT, "models", paste0(specification, "__", outcome, ".rds")
    )
    saveRDS(
      list(
        model_a = model_a,
        model_b = model_b,
        axes = axes,
        model_axes = model_axes,
        specification = specification,
        outcome = outcome,
        min_stratum_n = MIN_STRATUM_N,
        min_stratum_events = MIN_STRATUM_EVENTS
      ),
      model_file
    )
    register_output(model_file)
  }

  # Everything needed downstream has now been extracted, so the fitted objects
  # are released here rather than at the end of the loop. Under MCMC these are
  # by far the largest objects in the session, and holding one a moment longer
  # than necessary is what turns a long run into a failed one.
  rm(model_a, model_b)
  invisible(gc(verbose = FALSE))

  list(
    metrics = metrics,
    fixed_effects = fixed_effects,
    predictions = predictions,
    fit_attempts = fit_attempts,
    axis_decomposition = axis_decomposition,
    strata = strata
  )
}

# A reduced fit used only by the threshold sweep. It returns the headline
# variance measures and nothing else, so the sweep stays cheap.
fit_metrics_only <- function(strata, axes, references, specification,
                             specification_label, outcome, outcome_label,
                             threshold) {
  ruled <- apply_minimum_cell_rule(
    strata, axes, references, threshold, MIN_STRATUM_EVENTS
  )
  base <- data.frame(
    specification = specification,
    specification_label = specification_label,
    outcome = outcome,
    outcome_label = outcome_label,
    min_stratum_n = threshold,
    strata_retained = ruled$accounting$strata_retained,
    individuals_retained = ruled$accounting$individuals_retained,
    individuals_retained_percent =
      ruled$accounting$individuals_retained_percent,
    events_retained_percent = ruled$accounting$events_retained_percent
  )

  if (nrow(ruled$data) < MIN_RETAINED_STRATA ||
      length(ruled$model_axes) == 0L) {
    return(bind_cols(base, data.frame(
      variance_model_a = NA_real_, variance_model_b = NA_real_,
      vpc_model_a_percent = NA_real_, vpc_model_b_percent = NA_real_,
      pcv_percent = NA_real_, mor_model_a = NA_real_, mor_model_b = NA_real_,
      singular_model_b = NA, status = "too few strata"
    )))
  }

  formulae <- build_model_formulae(ruled$model_axes)
  attempt <- tryCatch({
    fit_a <- fit_glmer_robust(formulae$a, ruled$data, "sensitivity A")
    fit_b <- fit_glmer_robust(formulae$b, ruled$data, "sensitivity B")
    variance_a <- random_intercept_variance(fit_a$model)
    variance_b <- random_intercept_variance(fit_b$model)
    data.frame(
      variance_model_a = variance_a,
      variance_model_b = variance_b,
      vpc_model_a_percent = 100 * logistic_vpc(variance_a),
      vpc_model_b_percent = 100 * logistic_vpc(variance_b),
      pcv_percent = if (is.finite(variance_a) && variance_a > 0) {
        100 * (variance_a - variance_b) / variance_a
      } else NA_real_,
      mor_model_a = median_odds_ratio(variance_a),
      mor_model_b = median_odds_ratio(variance_b),
      singular_model_b = fit_b$singular,
      status = "fitted"
    )
  }, error = function(condition) {
    data.frame(
      variance_model_a = NA_real_, variance_model_b = NA_real_,
      vpc_model_a_percent = NA_real_, vpc_model_b_percent = NA_real_,
      pcv_percent = NA_real_, mor_model_a = NA_real_, mor_model_b = NA_real_,
      singular_model_b = NA,
      status = paste("error:", conditionMessage(condition))
    )
  })

  bind_cols(base, attempt)
}

# =============================================================================
# 7. FIGURE FUNCTIONS
# =============================================================================

# The typeface is resolved once, against the fonts actually installed. Asking
# for Arial on a Linux cluster otherwise produces a warning for every panel and
# silently substitutes a different family anyway.
resolve_font_family <- function(preferred) {
  if (HAS_SYSTEMFONTS) {
    available <- tryCatch(unique(systemfonts::system_fonts()$family),
                          error = function(condition) character(0))
    found <- preferred[preferred %in% available]
    if (length(found) > 0L) return(found[1L])
  }
  ""
}
# ggplot2 renamed the scale transformation argument from `trans` to `transform`
# in version 3.5.0, and passing the wrong one is a hard error rather than a
# warning. Analysis machines routinely run whichever version their institution
# froze, so the argument is selected at run time rather than assumed.
GGPLOT2_USES_TRANSFORM <- utils::packageVersion("ggplot2") >= "3.5.0"
scale_size_sqrt_compatible <- function(range) {
  if (GGPLOT2_USES_TRANSFORM) {
    scale_size_continuous(range = range, transform = "sqrt")
  } else {
    scale_size_continuous(range = range, trans = "sqrt")
  }
}

PLOT_FONT <- resolve_font_family(PREFERRED_FONTS)
if (!nzchar(PLOT_FONT)) {
  log_message("None of the preferred typefaces are installed; using the ",
              "graphics device default.", level = "NOTE")
} else {
  log_message("Figure typeface: ", PLOT_FONT)
}

# ggprism supplies clean, publication-oriented typography. Its default axes are
# deliberately strong, so they are overridden with 0.3-point lines and ticks for
# the lighter journal finish. When ggprism is unavailable the same finish is
# built from theme_classic, so figures remain consistent either way rather than
# the run failing over a cosmetic dependency.
base_theme <- if (HAS_GGPRISM) {
  ggprism::theme_prism(
    base_size = 12, base_family = PLOT_FONT, base_fontface = "plain",
    border = FALSE
  )
} else {
  theme_classic(base_size = 12, base_family = PLOT_FONT) +
    theme(
      axis.text = element_text(colour = "black", face = "bold", size = 10),
      axis.title = element_text(colour = "black", face = "bold", size = 12),
      legend.text = element_text(size = 10),
      plot.margin = margin(8, 8, 8, 8)
    )
}

publication_theme <- base_theme +
  theme(
    axis.line = element_line(linewidth = 0.30, colour = "black"),
    axis.ticks = element_line(linewidth = 0.30, colour = "black"),
    axis.ticks.length = grid::unit(2.2, "pt"),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    plot.subtitle = element_text(size = 10, colour = "grey25"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.spacing = grid::unit(1, "lines")
  )

# Restrained Nature-style palette. These colours are distinct in greyscale and
# remain clear for common forms of colour-vision deficiency.
nature_colours <- c(
  vermillion = "#E64B35",
  blue = "#3C5488",
  teal = "#00A087",
  coral = "#F39B7F",
  lavender = "#8491B4",
  cyan = "#4DBBD5",
  charcoal = "#2F2F2F",
  light_grey = "#B8B8B8"
)
colour_observed <- unname(nature_colours["vermillion"])
colour_predicted <- unname(nature_colours["blue"])
colour_interaction <- unname(nature_colours["teal"])
colour_charcoal <- unname(nature_colours["charcoal"])

# A Bayesian interval is a credible interval and a maximum-likelihood one is a
# confidence interval. They are not the same object and should not be labelled
# as though they were, so the wording follows the engine throughout the
# figures, the Word tables and the generated manuscript.
INTERVAL_LABEL <- if (ESTIMATION_ENGINE == "bayesian") {
  "95% credible interval"
} else {
  "95% confidence interval"
}
INTERVAL_SHORT <- if (ESTIMATION_ENGINE == "bayesian") "95% CrI" else "95% CI"
ESTIMATE_LABEL <- if (ESTIMATION_ENGINE == "bayesian") {
  if (POSTERIOR_ROBUST) "posterior median" else "posterior mean"
} else {
  "maximum likelihood estimate"
}

add_condition_labels <- function(data) {
  specification <- unique(data$specification)
  if (length(specification) != 1L || !specification %in% names(SPECIFICATIONS)) {
    stop("A single known specification is required to label strata.")
  }
  axes <- SPECIFICATIONS[[specification]]$axes
  pieces <- lapply(axes, function(axis) {
    paste0(display_name_for(axis), " = ", as.character(data[[axis]]))
  })
  data$condition_label <- do.call(paste, c(pieces, sep = "; "))

  # A multi-line version is used inside figures. IMD parenthetical wording is
  # shortened because the quintile number already identifies the category. The
  # axes are distributed over three lines regardless of how many there are, so
  # this does not silently break if a specification changes shape.
  compact_values <- vapply(axes, function(axis) {
    value <- as.character(data[[axis]])
    if (axis %in% c("imd_quintile", "imd_3cat")) {
      value <- sub(" \\(.*\\)$", "", value)
    }
    paste0(display_name_for(axis), " ", value)
  }, character(nrow(data)))
  if (is.null(dim(compact_values))) {
    compact_values <- matrix(compact_values, nrow = nrow(data))
  }
  groups <- split(seq_along(axes), cut(seq_along(axes), breaks = 3,
                                       labels = FALSE))
  lines <- lapply(groups, function(indices) {
    apply(compact_values[, indices, drop = FALSE], 1L, paste, collapse = "; ")
  })
  data$condition_plot_label <- do.call(paste, c(lines, sep = "\n"))
  data
}

# Select equal numbers from the low and high ends while avoiding duplicated
# rows when very few strata are available.
select_plot_extremes <- function(data, ordering_variable,
                                 number_each_end = N_PLOT_LABEL) {
  if (nrow(data) == 0L) return(data)
  ordered <- data[order(data[[ordering_variable]]), , drop = FALSE]
  number_each_end <- min(number_each_end, ceiling(nrow(ordered) / 2))
  unique(bind_rows(head(ordered, number_each_end),
                   tail(ordered, number_each_end)))
}

# A shared placeholder keeps the multi-panel layouts intact when a panel has
# nothing to show, instead of failing or silently changing the grid.
empty_panel <- function(title, subtitle, message_text) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = message_text,
             family = PLOT_FONT, size = 4.2, colour = colour_charcoal,
             lineheight = 1.25) +
    xlim(-1, 1) + ylim(-1, 1) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    publication_theme +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank())
}

make_figures <- function(predictions, metrics_row) {
  predictions <- add_condition_labels(predictions)
  title_suffix <- paste(
    unique(predictions$specification_label),
    unique(predictions$outcome_label), sep = ": "
  )
  cell_rule_caption <- paste0(
    "Strata with n < ", metrics_row$min_stratum_n,
    if (metrics_row$min_stratum_events > 0) {
      paste0(" or fewer than ", metrics_row$min_stratum_events, " events")
    } else "",
    " excluded: ", metrics_row$observed_strata, " of ",
    metrics_row$strata_before_rule, " strata retained (",
    formatC(metrics_row$individuals_retained_percent, format = "f", digits = 1),
    "% of individuals)."
  )
  model_b_singular <- isTRUE(metrics_row$singular_model_b)

  # Panel A parallels the tutorial's observed and precision-weighted
  # distribution comparison. Strata are weighted equally in this plot because
  # the object of interest is the distribution across strata.
  distribution_data <- bind_rows(
    data.frame(
      probability = predictions$observed_probability,
      type = "Observed stratum proportion"
    ),
    data.frame(
      probability = predictions$null_total_probability,
      type = "Model A precision-weighted prediction"
    )
  )
  p_distribution <- ggplot(distribution_data,
                           aes(x = probability, fill = type)) +
    geom_histogram(
      aes(y = after_stat(count / sum(count))),
      bins = 35, position = "identity", alpha = 0.55, colour = "white"
    ) +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_fill_manual(values = c(
      "Observed stratum proportion" = colour_observed,
      "Model A precision-weighted prediction" = colour_predicted
    )) +
    labs(
      title = "A. Observed and precision-weighted stratum risks",
      subtitle = title_suffix,
      x = "Event probability", y = "Percentage of strata", fill = NULL,
      caption = cell_rule_caption
    ) + publication_theme +
    theme(plot.caption = element_text(size = 7.5, colour = "grey35",
                                      hjust = 0))

  # Panel B is the tutorial-style caterpillar plot of Model B predictions.
  ordered_prediction <- predictions %>% arrange(total_probability)
  prediction_labels <- select_plot_extremes(
    ordered_prediction, "total_probability"
  )
  p_caterpillar <- ggplot(
    ordered_prediction,
    aes(x = prediction_rank, y = total_probability)
  ) +
    geom_linerange(
      aes(ymin = total_probability_low, ymax = total_probability_high),
      colour = "grey65", linewidth = 0.25, alpha = 0.7
    ) +
    # Point size carries the stratum size, so the reader can see immediately
    # whether a position on the gradient rests on many people or few. Without
    # it every stratum looks equally well established.
    geom_point(aes(size = n), colour = colour_predicted, alpha = 0.75) +
    scale_size_sqrt_compatible(range = c(0.5, 3.6)) +
    ggrepel::geom_text_repel(
      data = prediction_labels,
      aes(label = condition_plot_label),
      size = LABEL_SIZE_CALLOUT, family = PLOT_FONT, colour = colour_charcoal,
      min.segment.length = 0, segment.size = 0.25,
      box.padding = 0.55, point.padding = 0.25,
      max.overlaps = Inf, seed = 20260814, show.legend = FALSE
    ) +
    scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
    labs(
      title = "B. Predicted risk by intersectional stratum",
      subtitle = paste0("Model B: additive main effects plus stratum random ",
                        "effect. Bars are ", INTERVAL_SHORT,
                        "; point size is stratum n."),
      x = "Stratum rank", y = "Predicted event probability", size = "Stratum n"
    ) + publication_theme

  # Panel C shows the stratum interaction residual on the probability scale.
  # Values above zero indicate higher risk than predicted by additive effects;
  # values below zero indicate lower risk than predicted additively.
  ordered_interaction <- predictions %>%
    arrange(interaction_probability_difference)

  if (model_b_singular) {
    # With the variance at the boundary every residual is zero. Drawing the
    # usual panel would show a flat line at zero with no explanation, which
    # reads as a plotting failure rather than as the result it actually is.
    p_interaction <- empty_panel(
      "C. Intersectional interaction residuals",
      title_suffix,
      paste0(
        "Model B between-stratum variance estimated at zero.\n",
        "Additive main effects account for all between-stratum\n",
        "variation, so no residual interaction is present.\n",
        "PCV = ", formatC(metrics_row$pcv_percent, format = "f", digits = 1),
        "%."
      )
    )
  } else {
    interaction_label_pool <- ordered_interaction %>%
      filter(interaction_distinguishable)
    if (nrow(interaction_label_pool) == 0L) {
      interaction_label_pool <- ordered_interaction
    }
    interaction_labels <- select_plot_extremes(
      interaction_label_pool, "interaction_probability_difference",
      number_each_end = 1L
    )
    p_interaction <- ggplot(
      ordered_interaction,
      aes(x = interaction_rank, y = interaction_probability_difference)
    ) +
      geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey25") +
      geom_linerange(
        aes(ymin = interaction_probability_difference_low,
            ymax = interaction_probability_difference_high,
            colour = interaction_distinguishable),
        linewidth = 0.25, alpha = 0.75
      ) +
      geom_point(aes(colour = interaction_distinguishable), size = 0.7) +
      ggrepel::geom_text_repel(
        data = interaction_labels,
        aes(label = condition_plot_label),
        size = LABEL_SIZE_CALLOUT, family = PLOT_FONT,
        colour = colour_charcoal,
        min.segment.length = 0, segment.size = 0.25,
        box.padding = 0.55, point.padding = 0.25,
        max.overlaps = Inf, seed = 20260814, show.legend = FALSE
      ) +
      scale_colour_manual(
        values = c(`FALSE` = "grey55", `TRUE` = colour_interaction),
        labels = c(`FALSE` = paste(INTERVAL_SHORT, "includes zero"),
                   `TRUE` = paste(INTERVAL_SHORT, "excludes zero")),
        drop = FALSE
      ) +
      scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
      labs(
        title = "C. Intersectional interaction residuals",
        subtitle = paste0("Absolute risk due to interaction, with ",
                          INTERVAL_SHORT),
        x = "Stratum rank", y = "Probability difference", colour = NULL
      ) + publication_theme
  }

  # Panel E presents the interactions that survive multiplicity control as a
  # named horizontal forest plot. Only strata significant after Benjamini-
  # Hochberg correction are shown: with several hundred strata compared
  # simultaneously, an uncorrected interval excluding zero is not evidence of
  # an interaction, and plotting those strata alongside the real ones invites
  # them to be read as findings. The uncorrected flag is retained in
  # all_stratum_predictions.csv for anyone who needs it.
  significant_interactions <- predictions %>%
    filter(interaction_distinguishable_fdr) %>%
    arrange(interaction_probability_difference) %>%
    mutate(
      flagged = nzchar(artefact_flags),
      condition_plot_label = factor(
        condition_plot_label, levels = unique(condition_plot_label)
      )
    )

  n_uncorrected <- sum(predictions$interaction_distinguishable)

  if (nrow(significant_interactions) > 0L) {
    p_significant <- ggplot(
      significant_interactions,
      aes(x = interaction_probability_difference, y = condition_plot_label,
          colour = interaction_probability_difference > 0)
    ) +
      geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey35") +
      geom_errorbar(
        aes(xmin = interaction_probability_difference_low,
            xmax = interaction_probability_difference_high),
        orientation = "y", width = 0, linewidth = 0.45
      ) +
      geom_point(aes(shape = flagged), size = 2.6) +
      scale_shape_manual(
        values = c(`FALSE` = 16, `TRUE` = 1),
        labels = c(`FALSE` = "No artefact flag",
                   `TRUE` = "Flagged: check before interpreting"),
        drop = FALSE
      ) +
      scale_colour_manual(
        values = c(`FALSE` = unname(nature_colours["blue"]),
                   `TRUE` = unname(nature_colours["vermillion"])),
        labels = c(`FALSE` = "Lower than additive expectation",
                   `TRUE` = "Higher than additive expectation"),
        drop = FALSE
      ) +
      scale_x_continuous(labels = percent_format(accuracy = 0.1),
                         expand = expansion(mult = 0.06)) +
      labs(
        title = "E. Intersectional interaction residuals surviving FDR correction",
        subtitle = paste0(
          nrow(significant_interactions), " of ", nrow(predictions),
          " strata significant at FDR < ", FDR_LEVEL,
          " (", n_uncorrected, " before correction)"
        ),
        x = "Difference from additive predicted probability", y = NULL,
        colour = NULL, shape = NULL,
        caption = if (any(significant_interactions$flagged)) {
          paste0("Open circles carry an artefact flag (scale compression, ",
                 "precision-driven, or a marginal interval).\nSee ",
                 "interaction_quality_checks.csv before interpreting them.")
        } else {
          "No stratum carries an artefact flag."
        }
      ) + publication_theme +
      theme(axis.text.y = element_text(size = LABEL_SIZE_AXIS,
                                       lineheight = 1.05),
            legend.position = "bottom", legend.box = "vertical",
            plot.caption = element_text(size = 8, colour = "grey35",
                                        hjust = 0))
  } else {
    p_significant <- empty_panel(
      "E. Intersectional interaction residuals surviving FDR correction",
      title_suffix,
      if (model_b_singular) {
        paste0("Model B between-stratum variance estimated at zero.\n",
               "No interaction residuals to display.")
      } else {
        paste0(
          "No stratum interaction is significant at FDR < ", FDR_LEVEL, ".\n",
          n_uncorrected, " of ", nrow(predictions),
          " strata had an uncorrected interval excluding zero,\n",
          "which is consistent with chance across this many comparisons."
        )
      }
    )
  }

  # Panel D illustrates partial pooling directly. Small strata often have
  # volatile observed proportions and are pulled more strongly towards the
  # model prediction.
  p_shrinkage <- ggplot(
    predictions,
    aes(x = observed_probability, y = total_probability, size = n)
  ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.45, colour = "grey40") +
    geom_point(alpha = 0.55, colour = colour_predicted) +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_size_sqrt_compatible(range = c(0.6, 4.2)) +
    labs(
      title = "D. Partial pooling of observed stratum risks",
      subtitle = paste0(
        "VPC A = ", round(metrics_row$vpc_model_a_percent, 2),
        "%; VPC B = ", round(metrics_row$vpc_model_b_percent, 2),
        "%; PCV = ", round(metrics_row$pcv_percent, 1), "%"
      ),
      x = "Observed event proportion", y = "Model B predicted probability",
      size = "Stratum n"
    ) + publication_theme

  combined_five <- (p_distribution | p_caterpillar) /
    (p_interaction | p_shrinkage) /
    p_significant +
    plot_layout(heights = c(
      1, 1, max(1.5, nrow(significant_interactions) / 5)
    )) +
    plot_annotation(title = title_suffix, tag_levels = NULL)

  list(
    distribution = p_distribution,
    caterpillar = p_caterpillar,
    interaction = p_interaction,
    shrinkage = p_shrinkage,
    significant = p_significant,
    combined_five = combined_five,
    n_significant = nrow(significant_interactions)
  )
}

# Figure writing is wrapped so that a single failed device cannot end a run
# that has already produced valid model results.
save_plot_pair <- function(plot, path_without_extension, width, height) {
  png_path <- paste0(path_without_extension, ".png")
  status <- tryCatch({
    ggsave(png_path, plot, width = width, height = height,
           dpi = FIGURE_DPI, bg = "white", limitsize = FALSE)
    register_output(png_path)
    TRUE
  }, error = function(condition) {
    log_message("Could not write ", png_path, ": ",
                conditionMessage(condition), level = "WARN")
    FALSE
  })

  if (WRITE_SVG) {
    svg_path <- paste0(path_without_extension, ".svg")
    tryCatch({
      ggsave(svg_path, plot, width = width, height = height, bg = "white",
             limitsize = FALSE)
      register_output(svg_path)
    }, error = function(condition) {
      log_message("Could not write ", svg_path, ": ",
                  conditionMessage(condition), level = "WARN")
    })
  }
  invisible(status)
}

# =============================================================================
# 8. RUN EVERY SPECIFICATION AND OUTCOME
# =============================================================================

all_metrics <- list()
all_fixed_effects <- list()
all_predictions <- list()
all_extremes <- list()
all_significant <- list()
all_size_summaries <- list()
all_retention <- list()
all_excluded_strata <- list()
all_reference_changes <- list()
all_fit_attempts <- list()
all_axis_decomposition <- list()
all_failures <- list()
all_figure_failures <- list()
sensitivity_results <- list()

total_analyses <- length(SPECIFICATIONS) * length(OUTCOME_LABELS)
run_index <- 0L

for (specification_name in names(SPECIFICATIONS)) {
  specification <- SPECIFICATIONS[[specification_name]]
  axes <- specification$axes

  log_message("Preparing specification: ", specification$label)
  prepared <- prepare_analysis_data(
    analysis_master, axes, specification$references, names(OUTCOME_LABELS)
  )

  for (outcome in names(OUTCOME_LABELS)) {
    run_index <- run_index + 1L
    outcome_label <- unname(OUTCOME_LABELS[outcome])
    result_key <- paste(specification_name, outcome, sep = "__")
    table_prefix <- file.path(
      OUTPUT_ROOT, "tables", "per_analysis",
      paste0(specification_name, "__", outcome)
    )

    log_message("[", run_index, "/", total_analyses, "] ",
                specification$label, " - ", outcome_label)

    # Each analysis is isolated. A failure on one outcome is recorded and the
    # run continues, rather than discarding every result computed so far.
    analysis_status <- tryCatch({
      resume_files <- paste0(table_prefix,
                             c("__metrics.csv", "__fixed_effects.csv",
                               "__stratum_predictions.csv",
                               "__stratum_retention.csv"))
      if (RESUME && all(file.exists(resume_files))) {
        # Every per-analysis output is restored, not just the metrics, so the
        # consolidated tables written at the end of a resumed run are complete
        # rather than containing only the outcomes fitted in this session.
        log_message("  Existing results found; restoring and skipping ",
                    "(resume mode).")
        all_metrics[[result_key]] <- read.csv(resume_files[1L])
        all_fixed_effects[[result_key]] <- read.csv(resume_files[2L])
        restored_predictions <- read.csv(resume_files[3L])
        all_predictions[[result_key]] <- restored_predictions
        all_retention[[result_key]] <- read.csv(resume_files[4L])

        restored_ranked <- restored_predictions %>% arrange(total_probability)
        all_extremes[[result_key]] <- bind_rows(
          head(restored_ranked, N_EXTREME) %>%
            mutate(extreme = "Lowest predicted risk"),
          tail(restored_ranked, N_EXTREME) %>%
            mutate(extreme = "Highest predicted risk")
        )
        all_significant[[result_key]] <- restored_predictions %>%
          filter(interaction_distinguishable_fdr)
        axis_file <- paste0(table_prefix, "__axis_decomposition.csv")
        if (file.exists(axis_file)) {
          all_axis_decomposition[[result_key]] <- read.csv(axis_file)
        }
        "skipped"
      } else {
        # The observability filter is applied per outcome, so each model uses
        # exactly the records for which its outcome could have been observed.
        # prepare_analysis_data preserves row order and row count, so the
        # eligibility index computed on analysis_master aligns directly.
        eligibility_column <- unname(OUTCOME_ELIGIBILITY[outcome])
        if (!is.na(eligibility_column)) {
          eligible_rows <- !is.na(analysis_master[[eligibility_column]]) &
            analysis_master[[eligibility_column]] == 1
          outcome_data <- prepared[eligible_rows, , drop = FALSE]
          log_message("  Observable records: ",
                      format(sum(eligible_rows), big.mark = ","), " of ",
                      format(nrow(prepared), big.mark = ","), " (",
                      formatC(100 * mean(eligible_rows), format = "f",
                              digits = 1), "%), using ", eligibility_column, ".")
        } else {
          outcome_data <- prepared
        }
        if (nrow(outcome_data) == 0L) {
          stop("No observable records remain for ", outcome, ".")
        }

        counted <- make_stratum_counts(outcome_data, axes, outcome)

        all_size_summaries[[paste0(result_key, "__before")]] <-
          stratum_size_summary(
            counted$data, specification_name, specification$label,
            outcome, outcome_label, counted$included_n, counted$excluded_n,
            counted$possible_strata, stage = "before minimum cell rule"
          )

        ruled <- apply_minimum_cell_rule(
          counted$data, axes, specification$references,
          MIN_STRATUM_N, MIN_STRATUM_EVENTS
        )

        retention_row <- bind_cols(
          data.frame(
            specification = specification_name,
            specification_label = specification$label,
            outcome = outcome,
            outcome_label = outcome_label
          ),
          ruled$accounting
        )
        all_retention[[result_key]] <- retention_row

        log_message("  Minimum cell rule: ", ruled$accounting$strata_retained,
                    " of ", ruled$accounting$strata_before, " strata retained (",
                    formatC(ruled$accounting$individuals_retained_percent,
                            format = "f", digits = 1),
                    "% of individuals, ",
                    formatC(ruled$accounting$events_retained_percent,
                            format = "f", digits = 1), "% of events).")

        if (nrow(ruled$excluded) > 0L) {
          # The label is captured before the pipe. Assigning a `specification`
          # column inside mutate() shadows the outer `specification` list, so
          # referring to `specification$label` in a later argument of the same
          # mutate would resolve to the new character column.
          current_specification_label <- specification$label
          all_excluded_strata[[result_key]] <- ruled$excluded %>%
            mutate(specification = specification_name,
                   specification_label = current_specification_label,
                   outcome = outcome, outcome_label = outcome_label,
                   exclusion_reason = case_when(
                     n < MIN_STRATUM_N & events < MIN_STRATUM_EVENTS ~
                       "below size and event thresholds",
                     n < MIN_STRATUM_N ~ "below size threshold",
                     TRUE ~ "below event threshold"
                   ))
        }
        if (nrow(ruled$reference_changes) > 0L) {
          all_reference_changes[[result_key]] <- ruled$reference_changes %>%
            mutate(specification = specification_name, outcome = outcome)
          for (change_index in seq_len(nrow(ruled$reference_changes))) {
            change <- ruled$reference_changes[change_index, ]
            log_message("  Reference level '", change$intended_reference,
                        "' for ", change$variable, " did not survive the ",
                        "minimum cell rule; using '", change$used_reference,
                        "' instead.", level = "WARN")
          }
        }
        if (length(ruled$dropped_axes) > 0L) {
          log_message("  Axes reduced to a single category and removed from ",
                      "the fixed effects: ",
                      paste(ruled$dropped_axes, collapse = ", "),
                      level = "WARN")
        }

        if (nrow(ruled$data) < MIN_RETAINED_STRATA) {
          stop("Only ", nrow(ruled$data), " strata survived the minimum cell ",
               "rule (minimum required: ", MIN_RETAINED_STRATA,
               "). Lower MAIHDA_MIN_STRATUM_N or use the collapsed ",
               "specification for this outcome.")
        }
        if (length(ruled$model_axes) == 0L) {
          stop("No axis retained more than one category after the minimum ",
               "cell rule; Model B is not identifiable.")
        }

        all_size_summaries[[paste0(result_key, "__after")]] <-
          stratum_size_summary(
            ruled$data, specification_name, specification$label,
            outcome, outcome_label, sum(ruled$data$n), counted$excluded_n,
            counted$possible_strata, stage = "after minimum cell rule"
          )

        fit <- fit_one_maihda(
          ruled$data, axes, ruled$model_axes, specification_name,
          specification$label, outcome, outcome_label,
          counted$included_n, counted$excluded_n, counted$possible_strata,
          ruled$accounting, eligibility_column, nrow(outcome_data)
        )

        all_metrics[[result_key]] <- fit$metrics
        all_fixed_effects[[result_key]] <- fit$fixed_effects
        all_fit_attempts[[result_key]] <- fit$fit_attempts
        all_axis_decomposition[[result_key]] <- fit$axis_decomposition
        labelled_predictions <- add_condition_labels(fit$predictions)
        all_predictions[[result_key]] <- labelled_predictions

        ranked <- labelled_predictions %>% arrange(total_probability)
        all_extremes[[result_key]] <- bind_rows(
          head(ranked, N_EXTREME) %>% mutate(extreme = "Lowest predicted risk"),
          tail(ranked, N_EXTREME) %>% mutate(extreme = "Highest predicted risk")
        )
        # Matches what panel E plots: only strata surviving FDR correction are
        # presented as findings. Every stratum, with both the corrected and
        # uncorrected flags, remains in all_stratum_predictions.csv.
        all_significant[[result_key]] <- labelled_predictions %>%
          filter(interaction_distinguishable_fdr)

        # Per-analysis tables are written alongside the consolidated tables so
        # that a single model can be inspected or shared on its own.
        write_csv_output(fit$metrics, paste0(table_prefix, "__metrics.csv"))
        write_csv_output(fit$fixed_effects,
                         paste0(table_prefix, "__fixed_effects.csv"))
        write_csv_output(labelled_predictions,
                         paste0(table_prefix, "__stratum_predictions.csv"))
        write_csv_output(retention_row,
                         paste0(table_prefix, "__stratum_retention.csv"))
        if (!is.null(fit$axis_decomposition)) {
          write_csv_output(fit$axis_decomposition,
                           paste0(table_prefix, "__axis_decomposition.csv"))
        }

        # Figures are isolated from the analysis result. The model output is
        # already written by this point, so a device or layout failure is
        # recorded as a figure problem and must not mark a completed analysis
        # as failed.
        figure_status <- tryCatch({
        figure_list <- make_figures(fit$predictions, fit$metrics)
        figure_prefix <- file.path(
          OUTPUT_ROOT, "figures", "individual",
          paste0(specification_name, "__", outcome)
        )
        save_plot_pair(figure_list$distribution,
                       paste0(figure_prefix, "__A_distribution"),
                       FIGURE_WIDTH, FIGURE_HEIGHT)
        save_plot_pair(figure_list$caterpillar,
                       paste0(figure_prefix, "__B_caterpillar"),
                       FIGURE_WIDTH, FIGURE_HEIGHT)
        save_plot_pair(figure_list$interaction,
                       paste0(figure_prefix, "__C_interactions"),
                       FIGURE_WIDTH, FIGURE_HEIGHT)
        save_plot_pair(figure_list$shrinkage,
                       paste0(figure_prefix, "__D_shrinkage"),
                       FIGURE_WIDTH, FIGURE_HEIGHT)
        # The height expands with the number of retained interactions so every
        # condition remains legible without reducing the text to an
        # impractical size. Each stratum label occupies three lines, so the
        # per-stratum allowance is generous; the ceiling is high enough that a
        # long list stays readable rather than being compressed to fit. The
        # floor is kept low because FDR correction routinely leaves one or two
        # strata, and a tall panel holding a single row is mostly whitespace.
        significant_height <- max(
          3.4, min(30, 2.8 + 0.5 * figure_list$n_significant)
        )
        save_plot_pair(figure_list$significant,
                       paste0(figure_prefix, "__E_significant_interactions"),
                       width = 10.5, height = significant_height)
        save_plot_pair(
          figure_list$combined_five,
          file.path(OUTPUT_ROOT, "figures", "multipanel",
                    paste0(specification_name, "__", outcome, "__five_panel")),
          width = 15, height = 12 + significant_height
        )
        rm(figure_list)
        "ok"
        }, error = function(condition) {
          log_message("  Figures could not be produced: ",
                      conditionMessage(condition), level = "WARN")
          all_figure_failures[[result_key]] <<- data.frame(
            specification = specification_name,
            outcome = outcome,
            message = conditionMessage(condition)
          )
          "failed"
        })

        rm(fit, labelled_predictions)
        "completed"
      }
    }, error = function(condition) {
      log_message("  FAILED: ", conditionMessage(condition), level = "ERROR")
      all_failures[[result_key]] <<- data.frame(
        specification = specification_name,
        specification_label = specification$label,
        outcome = outcome,
        outcome_label = outcome_label,
        message = conditionMessage(condition)
      )
      "failed"
    })

    # A projection after each analysis turns "is this going to finish tonight?"
    # into a number, which matters when a run is left unattended.
    elapsed_so_far <- as.numeric(difftime(Sys.time(), RUN_STARTED_AT,
                                          units = "mins"))
    if (run_index < total_analyses && run_index > 0L) {
      projected_total <- elapsed_so_far / run_index * total_analyses
      log_message("  Progress: ", run_index, "/", total_analyses, "; elapsed ",
                  formatC(elapsed_so_far, format = "f", digits = 1),
                  " min; projected total ",
                  formatC(projected_total, format = "f", digits = 0),
                  " min (finish about ",
                  format(RUN_STARTED_AT + projected_total * 60, "%H:%M"), ").")
    }

    memory_used <- gc(verbose = FALSE)
    if (REPORT_MEMORY) {
      # gc() reports megabytes in a column whose name varies with the R build,
      # so it is matched rather than assumed by position.
      megabyte_column <- grep("Mb", colnames(memory_used), value = TRUE)[1L]
      megabytes <- if (!is.na(megabyte_column)) {
        sum(memory_used[, megabyte_column], na.rm = TRUE)
      } else NA_real_
      log_message("  Memory in use: ",
                  formatC(megabytes, format = "f", digits = 0),
                  " MB of R objects.")
    }
  }
}

completed_analyses <- length(all_metrics)
if (completed_analyses == 0L) {
  # Writing the failure log before stopping means the reason is preserved even
  # though no results exist.
  if (length(all_failures) > 0L) {
    write_csv_output(bind_rows(all_failures),
                     file.path(OUTPUT_ROOT, "tables", "failed_analyses.csv"))
  }
  stop("Every analysis failed. See ", LOG_PATH, " and failed_analyses.csv.")
}
log_message("Model fitting complete: ", completed_analyses, " of ",
            total_analyses, " analyses produced results.")

# =============================================================================
# 9. MINIMUM CELL SIZE THRESHOLD SENSITIVITY
# =============================================================================
# The threshold is a judgement call, so it is worth showing what it costs and
# what it changes. This refits Models A and B for a single outcome across a
# range of thresholds and reports the retained sample alongside the headline
# variance measures.

if (RUN_THRESHOLD_SENSITIVITY && SENSITIVITY_OUTCOME %in% names(OUTCOME_LABELS)) {
  log_message("Running minimum cell size sensitivity sweep for ",
              SENSITIVITY_OUTCOME, " across thresholds: ",
              paste(SENSITIVITY_THRESHOLDS, collapse = ", "))
  for (specification_name in names(SPECIFICATIONS)) {
    specification <- SPECIFICATIONS[[specification_name]]
    prepared <- prepare_analysis_data(
      analysis_master, specification$axes, specification$references,
      names(OUTCOME_LABELS)
    )
    counted <- tryCatch(
      make_stratum_counts(prepared, specification$axes, SENSITIVITY_OUTCOME),
      error = function(condition) NULL
    )
    if (is.null(counted)) next

    for (threshold in SENSITIVITY_THRESHOLDS) {
      key <- paste(specification_name, threshold, sep = "__")
      sensitivity_results[[key]] <- tryCatch(
        fit_metrics_only(
          counted$data, specification$axes, specification$references,
          specification_name, specification$label, SENSITIVITY_OUTCOME,
          unname(OUTCOME_LABELS[SENSITIVITY_OUTCOME]), threshold
        ),
        error = function(condition) NULL
      )
    }
  }
  log_message("Sensitivity sweep complete.")
}

sensitivity_table <- if (length(sensitivity_results) > 0L) {
  bind_rows(sensitivity_results)
} else {
  data.frame()
}

# =============================================================================
# 10. CONSOLIDATED TABLES
# =============================================================================

metrics_table <- bind_rows(all_metrics)
fixed_effects_table <- bind_rows(all_fixed_effects)
predictions_table <- bind_rows(all_predictions)
extremes_table <- bind_rows(all_extremes)
significant_table <- bind_rows(all_significant)
size_summary_table <- bind_rows(all_size_summaries)
retention_table <- bind_rows(all_retention)
excluded_strata_table <- bind_rows(all_excluded_strata)
reference_changes_table <- bind_rows(all_reference_changes)
fit_attempts_table <- bind_rows(all_fit_attempts)
axis_decomposition_table <- bind_rows(all_axis_decomposition)
failures_table <- bind_rows(all_failures)
figure_failures_table <- bind_rows(all_figure_failures)

write_csv_output(metrics_table,
                 file.path(OUTPUT_ROOT, "tables", "all_model_metrics.csv"))
write_csv_output(fixed_effects_table,
                 file.path(OUTPUT_ROOT, "tables",
                           "all_fixed_effects_odds_ratios.csv"))
write_csv_output(predictions_table,
                 file.path(OUTPUT_ROOT, "tables",
                           "all_stratum_predictions.csv"))
write_csv_output(extremes_table,
                 file.path(OUTPUT_ROOT, "tables", "extreme_strata.csv"))
write_csv_output(significant_table,
                 file.path(OUTPUT_ROOT, "tables",
                           "distinguishable_interaction_strata.csv"))
write_csv_output(size_summary_table,
                 file.path(OUTPUT_ROOT, "tables", "stratum_size_summary.csv"))
write_csv_output(retention_table,
                 file.path(OUTPUT_ROOT, "tables", "stratum_retention.csv"))
write_csv_output(fit_attempts_table,
                 file.path(OUTPUT_ROOT, "tables", "model_fit_attempts.csv"))

# Every stratum whose interval excludes zero, corrected or not, with the
# information needed to judge whether the finding is real. Written even when
# empty, because "nothing to check" is itself a useful thing to be able to
# point at.
interaction_quality_checks <- predictions_table %>%
  filter(interaction_distinguishable | interaction_distinguishable_fdr) %>%
  transmute(
    specification_label, outcome_label, stratum = condition_label,
    n, events,
    additive_probability, total_probability,
    absolute_risk_due_to_interaction,
    interaction_log_odds, interaction_q_value,
    survives_fdr = interaction_distinguishable_fdr,
    additive_probability_extreme, scale_compression_ratio,
    possible_ceiling_artefact, stratum_size_percentile, precision_driven,
    interval_excludes_zero_marginally, artefact_flags,
    verdict = case_when(
      !interaction_distinguishable_fdr ~ "does not survive FDR correction",
      possible_ceiling_artefact ~
        "likely scale compression near the probability boundary",
      precision_driven & abs(absolute_risk_due_to_interaction) < 0.01 ~
        "detected by precision; effect is small in absolute terms",
      interval_excludes_zero_marginally ~ "interval only marginally excludes zero",
      TRUE ~ "no artefact flag raised"
    )
  ) %>%
  arrange(specification_label, outcome_label,
          desc(abs(absolute_risk_due_to_interaction)))

write_csv_output(interaction_quality_checks,
                 file.path(OUTPUT_ROOT, "tables",
                           "interaction_quality_checks.csv"))

flagged_total <- sum(nzchar(interaction_quality_checks$artefact_flags) &
                       interaction_quality_checks$survives_fdr)
if (flagged_total > 0L) {
  log_message(flagged_total, " FDR-surviving ",
              if (flagged_total == 1L) "interaction carries" else
                "interactions carry",
              " an artefact flag. See interaction_quality_checks.csv.",
              level = "NOTE")
}
if (nrow(axis_decomposition_table) > 0L) {
  write_csv_output(axis_decomposition_table,
                   file.path(OUTPUT_ROOT, "tables",
                             "axis_variance_decomposition.csv"))
}
if (nrow(excluded_strata_table) > 0L) {
  write_csv_output(excluded_strata_table,
                   file.path(OUTPUT_ROOT, "tables", "excluded_strata.csv"))
}
if (nrow(reference_changes_table) > 0L) {
  write_csv_output(reference_changes_table,
                   file.path(OUTPUT_ROOT, "tables",
                             "reference_level_changes.csv"))
}
if (nrow(failures_table) > 0L) {
  write_csv_output(failures_table,
                   file.path(OUTPUT_ROOT, "tables", "failed_analyses.csv"))
}
if (nrow(figure_failures_table) > 0L) {
  write_csv_output(figure_failures_table,
                   file.path(OUTPUT_ROOT, "tables", "failed_figures.csv"))
}
if (nrow(sensitivity_table) > 0L) {
  write_csv_output(sensitivity_table,
                   file.path(OUTPUT_ROOT, "tables",
                             "minimum_cell_sensitivity.csv"))
}

# =============================================================================
# 11. APA-FORMATTED WORD TABLES
# =============================================================================

# The Word exports use officer for document structure and flextable for stable
# table layout. The visual treatment follows APA conventions: table number in
# bold, title in italics, no vertical rules, three restrained horizontal rules,
# repeated headers and a compact note beneath each table. Letter landscape is
# used because the model summaries are genuinely wide tabular material.

format_number <- function(x, digits = 2) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits,
                               big.mark = ","))
}

format_integer <- function(x) {
  ifelse(is.na(x), "", formatC(round(as.numeric(x)), format = "d",
                               big.mark = ","))
}

format_p_value <- function(x) {
  ifelse(is.na(x), "",
         ifelse(x < 0.001, "< .001",
                sub("^0", "", formatC(x, format = "f", digits = 3))))
}

# flextable validates opts_word strictly and the accepted keys have changed
# across releases, so an option that is merely cosmetic can otherwise abort the
# whole Word export. The richest supported set is used and the script falls
# back rather than failing. Header rows repeat across pages by default in Word
# output, so nothing important is lost by the fallback.
apply_word_table_properties <- function(table) {
  option_sets <- list(
    list(split = FALSE, keep_with_next = TRUE),
    list(split = FALSE),
    list()
  )
  for (options in option_sets) {
    attempt <- try(
      flextable::set_table_properties(
        table, layout = "autofit", width = 1, opts_word = options
      ),
      silent = TRUE
    )
    if (!inherits(attempt, "try-error")) return(attempt)
  }
  flextable::set_table_properties(table, layout = "autofit", width = 1)
}

make_apa_flextable <- function(data, font_size = 9) {
  table <- flextable::flextable(data)
  table <- flextable::theme_booktabs(table, bold_header = TRUE)
  table <- flextable::font(table, fontname = "Times New Roman", part = "all")
  table <- flextable::fontsize(table, size = font_size, part = "all")
  table <- flextable::align(table, align = "left", part = "all")
  descriptive_columns <- intersect(
    names(data),
    c("Outcome", "Specification", "Term", "Variable", "Level", "Stratum",
      "Extreme", "Status")
  )
  value_columns <- setdiff(names(data), descriptive_columns)
  if (length(value_columns) > 0L) {
    table <- flextable::align(
      table, j = value_columns, align = "right", part = "body"
    )
  }
  table <- flextable::valign(table, valign = "center", part = "all")
  table <- flextable::padding(
    table, padding.top = 3, padding.bottom = 3,
    padding.left = 4, padding.right = 4, part = "all"
  )
  table <- apply_word_table_properties(table)
  table <- flextable::autofit(table)
  flextable::fit_to_width(table, max_width = 9)
}

new_apa_table_document <- function() {
  document <- officer::read_docx()
  landscape_section <- officer::prop_section(
    page_size = officer::page_size(
      orient = "landscape", width = 11, height = 8.5
    ),
    page_margins = officer::page_mar(
      top = 1, bottom = 1, left = 1, right = 1,
      header = 0.492, footer = 0.492
    )
  )
  officer::body_set_default_section(document, landscape_section)
}

add_apa_table <- function(document, table_number, table_title, table_data,
                          table_note, add_page_break = FALSE,
                          font_size = 9) {
  if (add_page_break) {
    document <- officer::body_add_break(document)
    # A short spacer prevents LibreOffice and Word from pinning the following
    # table number against the printable page boundary after an explicit break.
    document <- officer::body_add_par(
      document, intToUtf8(160L), style = "Normal"
    )
  }

  number_style <- officer::fp_text(
    font.family = "Times New Roman", font.size = 12, bold = TRUE
  )
  title_style <- officer::fp_text(
    font.family = "Times New Roman", font.size = 12, italic = TRUE
  )
  note_label_style <- officer::fp_text(
    font.family = "Times New Roman", font.size = 9, italic = TRUE
  )
  note_style <- officer::fp_text(
    font.family = "Times New Roman", font.size = 9
  )

  document <- officer::body_add_fpar(
    document,
    officer::fpar(officer::ftext(paste("Table", table_number), number_style),
                  fp_p = officer::fp_par(
                    padding.top = 4, padding.bottom = 2,
                    keep_with_next = TRUE
                  ))
  )
  document <- officer::body_add_fpar(
    document,
    officer::fpar(officer::ftext(table_title, title_style),
                  fp_p = officer::fp_par(
                    padding.bottom = 6, keep_with_next = TRUE
                  ))
  )
  document <- flextable::body_add_flextable(
    document, make_apa_flextable(table_data, font_size = font_size)
  )
  document <- officer::body_add_fpar(
    document,
    officer::fpar(
      officer::ftext("Note. ", note_label_style),
      officer::ftext(table_note, note_style),
      fp_p = officer::fp_par(padding.top = 5, padding.bottom = 6)
    )
  )
  document
}

write_apa_table_document <- function(path, table_number, table_title,
                                     table_data, table_note,
                                     font_size = 9) {
  document <- new_apa_table_document()
  document <- add_apa_table(
    document, table_number, table_title, table_data, table_note,
    font_size = font_size
  )
  print(document, target = path)
  register_output(path)
  invisible(path)
}

# Prepare presentation versions separately from the analysis tables. This
# preserves full-precision CSV outputs while giving Word readers sensible
# rounding, descriptive headings and conventional p-value formatting.
word_cohort_characteristics <- cohort_characteristics %>%
  transmute(
    Variable = variable_label,
    Level = ifelse(is_reference, paste0(level, " (ref)"), level),
    `n` = format_integer(n),
    `%` = format_number(percent, 1),
    `Missing n` = format_integer(missing_n),
    `Missing %` = format_number(missing_percent, 2)
  )

word_outcome_counts <- outcome_counts %>%
  transmute(
    Outcome = outcome_label,
    `Complete N` = format_integer(complete_n),
    `Missing N` = format_integer(missing_n),
    Events = format_integer(event_n),
    `Non-events` = format_integer(non_event_n),
    `Event %` = format_number(event_percent, 2)
  )

word_event_rates <- descriptive_event_rates %>%
  transmute(
    Outcome = outcome_label,
    Variable = variable_label,
    Level = ifelse(is_reference, paste0(level, " (ref)"), level),
    `n` = format_integer(n),
    Events = format_integer(events),
    `Event %` = format_number(event_percent, 2),
    `95% CI` = paste0(format_number(ci_low_percent, 2), " to ",
                      format_number(ci_high_percent, 2))
  )

word_stratum_distribution <- descriptive_stratum_distribution %>%
  transmute(
    Specification = specification_label,
    `Possible strata` = format_integer(possible_strata),
    `Observed strata` = format_integer(observed_strata),
    `Empty strata` = format_integer(empty_strata),
    `Median n` = format_number(median_n, 1),
    `IQR` = paste0(format_number(q1_n, 0), " to ", format_number(q3_n, 0)),
    `Minimum n` = format_integer(minimum_n),
    `Maximum n` = format_integer(maximum_n),
    `n < 10` = format_integer(strata_below_10),
    `10-29` = format_integer(strata_10_to_29),
    `30-99` = format_integer(strata_30_to_99),
    `100+` = format_integer(strata_100_plus)
  )

word_retention <- retention_table %>%
  transmute(
    Specification = specification_label,
    Outcome = outcome_label,
    `Strata before` = format_integer(strata_before),
    `Strata retained` = format_integer(strata_retained),
    `Strata excluded` = format_integer(strata_excluded),
    `Individuals retained` = format_integer(individuals_retained),
    `Individuals retained %` = format_number(individuals_retained_percent, 1),
    `Events retained %` = format_number(events_retained_percent, 1),
    `Largest excluded n` = format_integer(excluded_max_n)
  )

word_axis_decomposition <- if (nrow(axis_decomposition_table) > 0L) {
  axis_decomposition_table %>%
    transmute(
      Outcome = outcome_label,
      Axis = axis_label,
      `VPC (%)` = format_number(vpc_model_2_percent, 2),
      `PCV (%)` = format_number(pcv_percent, 2),
      Status = ifelse(status == "fitted", "", status)
    )
} else {
  NULL
}

word_model_metrics <- metrics_table %>%
  transmute(
    Specification = specification_label,
    Outcome = outcome_label,
    `Strata` = format_integer(observed_strata),
    `Analysed N` = format_integer(analysed_n),
    `Event %` = format_number(event_percent, 2),
    `VPC A (%)` = format_number(vpc_model_a_percent, 2),
    `VPC A interval` = ifelse(
      is.na(vpc_model_a_low), "",
      paste0(format_number(vpc_model_a_low, 2), " to ",
             format_number(vpc_model_a_high, 2))
    ),
    `VPC B (%)` = format_number(vpc_model_b_percent, 2),
    `PCV (%)` = format_number(pcv_percent, 2),
    `PCV interval` = ifelse(
      is.na(pcv_low), "",
      paste0(format_number(pcv_low, 1), " to ", format_number(pcv_high, 1))
    ),
    `MOR A` = format_number(mor_model_a, 2),
    `MOR B` = format_number(mor_model_b, 2),
    `AUC additive` = format_number(auc_model_b_fixed, 3),
    `AUC total` = format_number(auc_model_b_total, 3),
    `Interactions (FDR)` = format_integer(n_distinguishable_fdr),
    Status = ifelse(singular_model_b, "Model B variance at zero", "Estimated")
  )

word_fixed_effects <- fixed_effects_table %>%
  transmute(
    Specification = specification_label,
    Outcome = outcome_label,
    Variable = variable_label,
    Level = level,
    OR = ifelse(is_reference, "1.00 (ref)", format_number(odds_ratio, 2)),
    `Interval` = ifelse(
      is_reference | is.na(confidence_low), "",
      paste0(format_number(confidence_low, 2), " to ",
             format_number(confidence_high, 2))
    ),
    `p` = format_p_value(p_value)
  ) %>%
  # A posterior has no p-value, so the column is dropped rather than left as a
  # row of blanks implying a test was performed.
  { if (ESTIMATION_ENGINE == "bayesian") select(., -p) else . }

word_extremes <- extremes_table %>%
  transmute(
    Specification = specification_label,
    Outcome = outcome_label,
    Extreme = extreme,
    Stratum = condition_label,
    `n` = format_integer(n),
    Events = format_integer(events),
    `Observed %` = format_number(100 * observed_probability, 2),
    `Predicted %` = format_number(100 * total_probability, 2),
    `95% CI` = paste0(format_number(100 * total_probability_low, 2), " to ",
                      format_number(100 * total_probability_high, 2))
  )

word_significant <- if (nrow(significant_table) > 0L) {
  significant_table %>%
    transmute(
      Specification = specification_label,
      Outcome = outcome_label,
      Stratum = condition_label,
      `n` = format_integer(n),
      `Additive %` = format_number(100 * additive_probability, 2),
      `Total %` = format_number(100 * total_probability, 2),
      `Difference (pp)` = format_number(
        100 * interaction_probability_difference, 2
      ),
      `95% CI (pp)` = paste0(
        format_number(100 * interaction_probability_difference_low, 2),
        " to ",
        format_number(100 * interaction_probability_difference_high, 2)
      ),
      `q` = format_p_value(interaction_q_value)
    ) %>%
    { if (ESTIMATION_ENGINE == "bayesian" &&
          "interaction_probability_of_direction" %in% names(significant_table)) {
        mutate(., `Prob. of direction` = format_number(
          100 * significant_table$interaction_probability_of_direction, 1))
      } else . }
} else {
  data.frame(
    Result = paste0(
      "No stratum interaction was significant at FDR < ", FDR_LEVEL, "."
    )
  )
}

word_sensitivity <- if (nrow(sensitivity_table) > 0L) {
  sensitivity_table %>%
    transmute(
      Specification = specification_label,
      Outcome = outcome_label,
      `Minimum n` = format_integer(min_stratum_n),
      `Strata retained` = format_integer(strata_retained),
      `Individuals retained %` = format_number(individuals_retained_percent, 1),
      `Events retained %` = format_number(events_retained_percent, 1),
      `VPC A (%)` = format_number(vpc_model_a_percent, 2),
      `VPC B (%)` = format_number(vpc_model_b_percent, 2),
      `PCV (%)` = format_number(pcv_percent, 2),
      `MOR A` = format_number(mor_model_a, 2)
    )
} else {
  NULL
}

cell_rule_note <- paste0(
  "Strata with fewer than ", MIN_STRATUM_N, " individuals",
  if (MIN_STRATUM_EVENTS > 0) {
    paste0(", or fewer than ", MIN_STRATUM_EVENTS, " events,")
  } else "",
  " were excluded before model fitting."
)

axis_list_text <- paste(
  display_name_for(SPECIFICATIONS[[1L]]$axes), collapse = ", "
)

word_tables <- list(
  list(number = 1, title = "Cohort Characteristics",
       data = word_cohort_characteristics,
       note = paste0(
         "Percentages are of the records with that variable observed. ",
         "Intersectional strata are defined by ", axis_list_text,
         ". Categories marked (ref) are the reference categories in Model B."
       ),
       filename = "Table_1_cohort_characteristics.docx", font_size = 9),
  list(number = 2, title = "Outcome Frequencies",
       data = word_outcome_counts,
       note = "Percentages use all complete observations for each outcome, before the minimum cell size rule was applied.",
       filename = "Table_2_outcome_frequencies.docx", font_size = 10),
  list(number = 3, title = "Crude Event Rates by Category",
       data = word_event_rates,
       note = paste(
         "Unadjusted event rates with Wilson score 95% confidence intervals.",
         "The Wilson interval is used because several categories are small or",
         "have low event counts, where the Wald interval performs poorly.",
         "These are marginal rates and take no account of the other axes."
       ),
       filename = "Table_3_event_rates_by_category.docx", font_size = 7.5),
  list(number = 4, title = "Intersectional Stratum Size Distribution",
       data = word_stratum_distribution,
       note = paste(
         "Distribution of stratum sizes before the minimum cell size rule.",
         "Empty strata are combinations that are possible in principle but",
         "were not observed in the cohort."
       ),
       filename = "Table_4_stratum_size_distribution.docx", font_size = 9),
  list(number = 5, title = "Stratum Retention Under the Minimum Cell Size Rule",
       data = word_retention,
       note = paste(
         cell_rule_note,
         "Zero-event strata were retained, because a stratum with no events is",
         "an observation of low risk rather than missing information."
       ),
       filename = "Table_5_stratum_retention.docx", font_size = 8),
  list(number = 6, title = "MAIHDA Model Summary",
       data = word_model_metrics,
       note = paste(
         "VPC = variance partition coefficient; PCV = proportional change in",
         "variance; MOR = median odds ratio; AUC = area under the receiver",
         "operating characteristic curve; FDR = false discovery rate. Model A",
         "is the null model and Model B contains additive main effects plus",
         "the stratum random intercept. A Model B variance estimated at zero",
         "indicates that additive main effects account for all between-stratum",
         "variation, so no residual intersectional interaction was detected."
       ),
       filename = "Table_6_model_summary.docx", font_size = 8),
  list(number = 7, title = "Additive Fixed Effects from Model B",
       data = word_fixed_effects,
       note = "OR = odds ratio; CI = confidence interval. Reference categories are shown explicitly with an odds ratio of 1.00.",
       filename = "Table_7_fixed_effects.docx", font_size = 8),
  list(number = 8, title = "Lowest- and Highest-Risk Intersectional Strata",
       data = word_extremes,
       note = paste0("The ", N_EXTREME,
                     " lowest and highest Model B predicted-risk strata are shown for each analysis. ",
                     cell_rule_note),
       filename = "Table_8_extreme_strata.docx", font_size = 7.5),
  list(number = 9,
       title = "Intersectional Interaction Residuals Surviving FDR Correction",
       data = word_significant,
       note = paste0(
         "Only strata significant after Benjamini-Hochberg correction across ",
         "the strata within each analysis are shown, because several hundred ",
         "strata are compared simultaneously and an uncorrected interval ",
         "excluding zero is not on its own evidence of an interaction. ",
         "Difference is the Model B total predicted probability minus its ",
         "additive-only prediction, expressed in percentage points. q is the ",
         "adjusted p value; all strata are reported with both the corrected ",
         "and uncorrected flags in all_stratum_predictions.csv."
       ),
       filename = "Table_9_significant_interactions.docx", font_size = 7.5)
)

word_diagnostics <- if (ESTIMATION_ENGINE == "bayesian" &&
                        any(!is.na(metrics_table$max_rhat_model_a))) {
  metrics_table %>%
    transmute(
      Outcome = outcome_label,
      `R-hat A` = format_number(max_rhat_model_a, 3),
      `R-hat B` = format_number(max_rhat_model_b, 3),
      `Min ESS A` = format_integer(min_ess_model_a),
      `Min ESS B` = format_integer(min_ess_model_b),
      `Divergences A` = format_integer(divergent_model_a),
      `Divergences B` = format_integer(divergent_model_b),
      `Fit minutes` = format_number(fit_minutes, 1),
      Status = ifelse(
        pmax(max_rhat_model_a, max_rhat_model_b, na.rm = TRUE) > MAX_RHAT |
          pmin(min_ess_model_a, min_ess_model_b, na.rm = TRUE) <
            MIN_EFFECTIVE_SAMPLE_SIZE,
        "Check", "OK"
      )
    )
} else {
  NULL
}

if (!is.null(word_diagnostics)) {
  word_tables[[length(word_tables) + 1L]] <- list(
    number = length(word_tables) + 1L,
    title = "Markov Chain Monte Carlo Convergence Diagnostics",
    data = word_diagnostics,
    note = paste0(
      "R-hat compares within- and between-chain variance and should be below ",
      MAX_RHAT, "; effective sample size should exceed ",
      MIN_EFFECTIVE_SAMPLE_SIZE,
      "; divergent transitions should be zero. Models were fitted with ",
      MCMC_CHAINS, " chains of ", MCMC_ITERATIONS, " iterations (",
      MCMC_WARMUP, " warmup), adapt_delta ", MCMC_ADAPT_DELTA,
      ". Rows marked Check did not meet a threshold and should be refitted ",
      "with more iterations before the results are relied upon."
    ),
    filename = "Table_12_mcmc_diagnostics.docx", font_size = 8
  )
}

if (!is.null(word_axis_decomposition)) {
  word_tables[[length(word_tables) + 1L]] <- list(
    number = length(word_tables) + 1L,
    title = "Contribution of Each Axis to Between-Stratum Variance",
    data = word_axis_decomposition,
    note = paste(
      "Each row reports a partially adjusted model containing the null model",
      "plus a single axis. PCV is the proportional reduction in between-stratum",
      "variance attributable to that axis alone. Contributions overlap and do",
      "not sum to the fully adjusted PCV, because the axes are correlated in",
      "the population. This table is secondary and descriptive: attending to",
      "individual axis contributions invites a reversion to single-axis",
      "thinking, which is counter to the purpose of an intersectional",
      "analysis. It should not displace the collective additive effect."
    ),
    filename = "Table_10_axis_variance_decomposition.docx", font_size = 8
  )
}

if (!is.null(word_sensitivity)) {
  word_tables[[length(word_tables) + 1L]] <- list(
    number = length(word_tables) + 1L,
    title = "Sensitivity to the Minimum Cell Size Threshold",
    data = word_sensitivity,
    note = paste(
      "Models A and B refitted at each candidate minimum stratum size for a",
      "single outcome. A threshold that materially changes VPC or PCV while",
      "retaining a similar share of the cohort indicates that the result is",
      "sensitive to the trimming rule and should be reported as such."
    ),
    filename = "Table_11_minimum_cell_sensitivity.docx", font_size = 8
  )
}

word_directory <- file.path(OUTPUT_ROOT, "tables", "word")
for (table_definition in word_tables) {
  status <- tryCatch({
    write_apa_table_document(
      path = file.path(word_directory, table_definition$filename),
      table_number = table_definition$number,
      table_title = table_definition$title,
      table_data = table_definition$data,
      table_note = table_definition$note,
      font_size = table_definition$font_size
    )
    TRUE
  }, error = function(condition) {
    log_message("Could not write ", table_definition$filename, ": ",
                conditionMessage(condition), level = "WARN")
    FALSE
  })
}

# A combined table book is convenient for review and submission preparation.
tryCatch({
  combined_word_document <- new_apa_table_document()
  for (table_index in seq_along(word_tables)) {
    table_definition <- word_tables[[table_index]]
    combined_word_document <- add_apa_table(
      combined_word_document,
      table_definition$number, table_definition$title,
      table_definition$data, table_definition$note,
      add_page_break = table_index > 1L,
      font_size = table_definition$font_size
    )
  }
  combined_path <- file.path(word_directory, "MAIHDA_APA_tables.docx")
  print(combined_word_document, target = combined_path)
  register_output(combined_path)
}, error = function(condition) {
  log_message("Could not write the combined table book: ",
              conditionMessage(condition), level = "WARN")
})

# =============================================================================
# 12. CROSS-OUTCOME AND CROSS-SPECIFICATION FIGURES
# =============================================================================

metrics_long <- metrics_table %>%
  select(specification_label, outcome_label,
         vpc_model_a_percent, vpc_model_b_percent) %>%
  pivot_longer(
    cols = c(vpc_model_a_percent, vpc_model_b_percent),
    names_to = "model", values_to = "vpc_percent"
  ) %>%
  mutate(
    model = recode(
      model,
      vpc_model_a_percent = "Model A: null",
      vpc_model_b_percent = "Model B: additive"
    )
  )

p_vpc <- ggplot(
  metrics_long,
  aes(x = vpc_percent, y = reorder(outcome_label, vpc_percent),
      colour = model, shape = specification_label)
) +
  geom_point(size = 2.5, position = position_dodge(width = 0.55)) +
  facet_wrap(~ specification_label, scales = "free_y") +
  scale_colour_manual(values = c(
    "Model A: null" = colour_observed,
    "Model B: additive" = colour_predicted
  )) +
  labs(
    title = "A. Variance partition coefficients across outcomes",
    x = "Latent-scale VPC (%)", y = NULL, colour = NULL, shape = NULL
  ) + publication_theme + theme(legend.position = "bottom")

auc_long <- metrics_table %>%
  select(specification_label, outcome_label,
         auc_model_a_total, auc_model_b_total, auc_model_b_fixed) %>%
  pivot_longer(
    cols = starts_with("auc_"), names_to = "prediction", values_to = "auc"
  ) %>%
  mutate(prediction = recode(
    prediction,
    auc_model_a_total = "Model A: stratum",
    auc_model_b_total = "Model B: additive + interaction",
    auc_model_b_fixed = "Model B: additive only"
  ))

p_auc <- ggplot(
  auc_long,
  aes(x = auc, y = reorder(outcome_label, auc), colour = prediction)
) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey55") +
  geom_point(size = 2.3, position = position_dodge(width = 0.55)) +
  facet_wrap(~ specification_label, scales = "free_y") +
  scale_x_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, 0.1)) +
  scale_colour_manual(values = c(
    "Model A: stratum" = colour_observed,
    "Model B: additive + interaction" = colour_predicted,
    "Model B: additive only" = colour_interaction
  )) +
  labs(
    title = "B. Discriminatory accuracy across outcomes",
    x = "Area under the ROC curve", y = NULL, colour = NULL
  ) + publication_theme + theme(legend.position = "bottom")

p_pcv <- ggplot(
  metrics_table,
  aes(x = pcv_percent, y = reorder(outcome_label, pcv_percent),
      colour = specification_label)
) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey55") +
  geom_point(size = 2.5) +
  labs(
    title = "C. Proportional change in between-stratum variance",
    x = "PCV from Model A to Model B (%)", y = NULL, colour = NULL
  ) + publication_theme

mor_long <- metrics_table %>%
  select(specification_label, outcome_label, mor_model_a, mor_model_b) %>%
  pivot_longer(starts_with("mor_"), names_to = "model", values_to = "mor") %>%
  mutate(model = recode(model,
                        mor_model_a = "Model A: null",
                        mor_model_b = "Model B: additive"))

p_mor <- ggplot(
  mor_long,
  aes(x = mor, y = reorder(outcome_label, mor), colour = model,
      shape = specification_label)
) +
  geom_vline(xintercept = 1, linewidth = 0.35, colour = "grey55") +
  geom_point(size = 2.5, position = position_dodge(width = 0.55)) +
  labs(
    title = "D. Median odds ratios across outcomes",
    x = "Median odds ratio", y = NULL, colour = NULL, shape = NULL
  ) + publication_theme

# The four panel-specific headings already identify the content clearly. An
# additional overall heading competes for space once all twelve outcome labels
# are present, so the combined version deliberately uses the panel headings.
summary_multipanel <- (p_vpc | p_auc) / (p_pcv | p_mor)

save_plot_pair(p_vpc,
               file.path(OUTPUT_ROOT, "figures", "individual", "summary__vpc"),
               width = 12, height = 8)
save_plot_pair(p_auc,
               file.path(OUTPUT_ROOT, "figures", "individual", "summary__auc"),
               width = 12, height = 8)
save_plot_pair(p_pcv,
               file.path(OUTPUT_ROOT, "figures", "individual", "summary__pcv"),
               width = 9, height = 7)
save_plot_pair(p_mor,
               file.path(OUTPUT_ROOT, "figures", "individual", "summary__mor"),
               width = 10, height = 7)
save_plot_pair(summary_multipanel,
               file.path(OUTPUT_ROOT, "figures", "multipanel",
                         "all_outcomes__summary"),
               width = 18, height = 14)

# --- Minimum cell size figures ----------------------------------------------
# The trimming rule deserves its own visual account, because a reader's first
# question about any trimmed intersectional analysis is what was removed.

retention_long <- retention_table %>%
  select(specification_label, outcome_label,
         `Strata retained` = strata_retained_percent,
         `Individuals retained` = individuals_retained_percent,
         `Events retained` = events_retained_percent) %>%
  pivot_longer(cols = c(`Strata retained`, `Individuals retained`,
                        `Events retained`),
               names_to = "quantity", values_to = "percent")

p_retention <- ggplot(
  retention_long,
  aes(x = percent, y = reorder(outcome_label, percent), colour = quantity)
) +
  geom_point(size = 2.4, position = position_dodge(width = 0.55)) +
  facet_wrap(~ specification_label, scales = "free_y") +
  scale_x_continuous(limits = c(0, 100)) +
  scale_colour_manual(values = c(
    "Strata retained" = colour_interaction,
    "Individuals retained" = colour_predicted,
    "Events retained" = colour_observed
  )) +
  labs(
    title = "Retention under the minimum cell size rule",
    subtitle = cell_rule_note,
    x = "Percentage retained", y = NULL, colour = NULL
  ) + publication_theme

save_plot_pair(p_retention,
               file.path(OUTPUT_ROOT, "figures", "individual",
                         "summary__stratum_retention"),
               width = 12, height = 8)

# --- Which axis drives the between-stratum variance -------------------------
# The single most useful decomposition in the whole analysis: Model B tells you
# that between-stratum variation is largely additive, and this tells you what
# it is additive in.
if (nrow(axis_decomposition_table) > 0L) {
  p_axis_decomposition <- ggplot(
    axis_decomposition_table %>% filter(status == "fitted"),
    aes(x = pcv_percent, y = reorder(axis_label, pcv_percent),
        colour = outcome_label)
  ) +
    geom_point(size = 2.4, alpha = 0.85) +
    scale_x_continuous(limits = c(0, NA)) +
    labs(
      title = "Contribution of each axis to between-stratum variance",
      subtitle = paste(
        "Proportional change in variance when each axis alone is added to the",
        "null model.\nContributions overlap and are not additive, because the",
        "axes are correlated in the population."
      ),
      x = "PCV from the null model (%)", y = NULL, colour = NULL
    ) + publication_theme +
    theme(legend.position = "right",
          legend.text = element_text(size = 8))

  save_plot_pair(p_axis_decomposition,
                 file.path(OUTPUT_ROOT, "figures", "individual",
                           "summary__axis_variance_decomposition"),
                 width = 12, height = 7)
}

# --- VPC in the context of the published literature -------------------------
# A VPC is difficult to judge in isolation. Plotting the observed values
# against the range reported across published MAIHDA health analyses gives the
# reader an immediate sense of whether this cohort is unusual.
p_vpc_benchmark <- ggplot(
  metrics_table,
  aes(x = vpc_model_a_percent, y = reorder(outcome_label, vpc_model_a_percent))
) +
  annotate("rect", xmin = VPC_BENCHMARK_LOW, xmax = VPC_BENCHMARK_HIGH,
           ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.5) +
  geom_vline(xintercept = VPC_BENCHMARK_MEDIAN, linetype = "dashed",
             linewidth = 0.4, colour = "grey35") +
  {
    if (all(is.na(metrics_table$vpc_model_a_low))) {
      NULL
    } else {
      geom_errorbar(aes(xmin = vpc_model_a_low, xmax = vpc_model_a_high),
                    orientation = "y", width = 0, linewidth = 0.4,
                    colour = colour_predicted)
    }
  } +
  geom_point(size = 2.6, colour = colour_predicted) +
  labs(
    title = "Between-stratum variation against published MAIHDA analyses",
    subtitle = paste0(
      "Shaded band: range of VPCs reported across published MAIHDA health ",
      "analyses (", VPC_BENCHMARK_LOW, "% to ", VPC_BENCHMARK_HIGH,
      "%).\nDashed line: their median (", VPC_BENCHMARK_MEDIAN,
      "%). Bars are profile-likelihood 95% intervals where available."
    ),
    x = "Model A VPC (%)", y = NULL
  ) + publication_theme

save_plot_pair(p_vpc_benchmark,
               file.path(OUTPUT_ROOT, "figures", "individual",
                         "summary__vpc_benchmark"),
               width = 11, height = 7)

# --- Variance components with their intervals -------------------------------
# The single clearest advantage of the Bayesian fit is that the VPC and PCV
# arrive with intervals, so they can be plotted as estimates with uncertainty
# rather than as bare points. Drawn whenever intervals exist, which under
# maximum likelihood means whenever profiling succeeded.
if (any(!is.na(metrics_table$vpc_model_a_low))) {
  variance_components <- bind_rows(
    metrics_table %>%
      transmute(outcome_label,
                quantity = "VPC, Model A (null)",
                estimate = vpc_model_a_percent,
                low = vpc_model_a_low, high = vpc_model_a_high),
    metrics_table %>%
      transmute(outcome_label,
                quantity = "VPC, Model B (additive)",
                estimate = vpc_model_b_percent,
                low = vpc_model_b_low, high = vpc_model_b_high),
    metrics_table %>%
      transmute(outcome_label,
                quantity = "PCV, A to B",
                estimate = pcv_percent,
                low = pcv_low, high = pcv_high)
  )

  p_variance_components <- ggplot(
    variance_components,
    aes(x = estimate, y = reorder(outcome_label, estimate))
  ) +
    geom_errorbar(aes(xmin = low, xmax = high), orientation = "y",
                  width = 0, linewidth = 0.45, colour = colour_predicted) +
    geom_point(size = 2.4, colour = colour_predicted) +
    facet_wrap(~ quantity, scales = "free_x") +
    labs(
      title = "Variance components with uncertainty",
      subtitle = paste0(ESTIMATE_LABEL, " and ", INTERVAL_LABEL,
                        ". Note the differing horizontal scales."),
      x = "Percent", y = NULL
    ) + publication_theme

  save_plot_pair(p_variance_components,
                 file.path(OUTPUT_ROOT, "figures", "individual",
                           "summary__variance_components"),
                 width = 13, height = 7)
}

# --- MCMC diagnostics -------------------------------------------------------
# Convergence is a property of the run, not an assumption, so it is plotted
# rather than left buried in a column. Anything the wrong side of the
# thresholds is coloured so it cannot be missed at a glance.
if (ESTIMATION_ENGINE == "bayesian" &&
    any(!is.na(metrics_table$max_rhat_model_a))) {
  diagnostics_long <- bind_rows(
    metrics_table %>%
      transmute(outcome_label, model = "Model A",
                rhat = max_rhat_model_a, ess = min_ess_model_a,
                divergences = divergent_model_a),
    metrics_table %>%
      transmute(outcome_label, model = "Model B",
                rhat = max_rhat_model_b, ess = min_ess_model_b,
                divergences = divergent_model_b)
  )

  p_rhat <- ggplot(diagnostics_long,
                   aes(x = rhat, y = reorder(outcome_label, rhat),
                       colour = rhat > MAX_RHAT, shape = model)) +
    geom_vline(xintercept = MAX_RHAT, linetype = "dashed", linewidth = 0.4,
               colour = "grey45") +
    geom_point(size = 2.4, position = position_dodge(width = 0.5)) +
    scale_colour_manual(values = c(`FALSE` = colour_interaction,
                                   `TRUE` = colour_observed),
                        labels = c(`FALSE` = "Within threshold",
                                   `TRUE` = "Above threshold"),
                        drop = FALSE) +
    labs(title = "A. Convergence (R-hat)",
         subtitle = paste0("Dashed line at ", MAX_RHAT),
         x = "Maximum R-hat", y = NULL, colour = NULL, shape = NULL) +
    publication_theme

  p_ess <- ggplot(diagnostics_long,
                  aes(x = ess, y = reorder(outcome_label, ess),
                      colour = ess < MIN_EFFECTIVE_SAMPLE_SIZE, shape = model)) +
    geom_vline(xintercept = MIN_EFFECTIVE_SAMPLE_SIZE, linetype = "dashed",
               linewidth = 0.4, colour = "grey45") +
    geom_point(size = 2.4, position = position_dodge(width = 0.5)) +
    scale_colour_manual(values = c(`FALSE` = colour_interaction,
                                   `TRUE` = colour_observed),
                        labels = c(`FALSE` = "Adequate",
                                   `TRUE` = "Below threshold"),
                        drop = FALSE) +
    labs(title = "B. Effective sample size",
         subtitle = paste0("Dashed line at ", MIN_EFFECTIVE_SAMPLE_SIZE),
         x = "Minimum effective sample size", y = NULL, colour = NULL,
         shape = NULL) +
    publication_theme

  p_diagnostics <- p_rhat | p_ess
  save_plot_pair(p_diagnostics,
                 file.path(OUTPUT_ROOT, "figures", "individual",
                           "summary__mcmc_diagnostics"),
                 width = 14, height = 7)
}

# --- Descriptive figures ----------------------------------------------------
# Drawn here rather than in section 3B because the publication theme is defined
# with the other figure code.

# Levels are keyed by variable before being made a factor. Two axes can share a
# category name -- "Other" appears in more than one classification -- and a
# plain level factor would silently merge them across facets.
LEVEL_KEY_SEPARATOR <- ""
key_levels <- function(data) {
  data %>%
    mutate(
      level_key = factor(
        paste(variable, level, sep = LEVEL_KEY_SEPARATOR),
        levels = rev(unique(paste(variable, level, sep = LEVEL_KEY_SEPARATOR)))
      ),
      # Facets follow the order the axes are declared in the specification,
      # not alphabetical order, so every descriptive figure reads in the same
      # sequence as the tables and the model output.
      variable_label = factor(variable_label,
                              levels = display_name_for(descriptive_axes))
    )
}
strip_level_key <- function(x) {
  sub(paste0(".*", LEVEL_KEY_SEPARATOR), "", x)
}

p_cohort <- ggplot(
  key_levels(cohort_characteristics),
  aes(x = percent, y = level_key, fill = is_reference)
) +
  geom_col(width = 0.7) +
  geom_text(aes(label = format_integer(n)), hjust = -0.15, size = 3,
            family = PLOT_FONT, colour = colour_charcoal) +
  facet_wrap(~ variable_label, scales = "free_y", ncol = 2) +
  scale_y_discrete(labels = strip_level_key) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  scale_fill_manual(values = c(`FALSE` = colour_predicted,
                               `TRUE` = colour_interaction),
                    labels = c(`FALSE` = "Category",
                               `TRUE` = "Reference category"),
                    drop = FALSE) +
  labs(
    title = "Cohort composition by stratum axis",
    subtitle = paste0("N = ", format(nrow(analysis_master), big.mark = ","),
                      "; bars show the percentage within each variable"),
    x = "Percentage of cohort", y = NULL, fill = NULL
  ) + publication_theme

save_plot_pair(p_cohort,
               file.path(OUTPUT_ROOT, "figures", "individual",
                         "descriptive__cohort_composition"),
               width = 11, height = 9)

# The stratum size distribution is shown on a log scale because the counts span
# several orders of magnitude and a linear axis would compress everything
# except the largest strata into a single bar.
p_stratum_sizes <- ggplot(descriptive_strata, aes(x = n)) +
  geom_histogram(bins = 40, fill = colour_predicted, colour = "white",
                 alpha = 0.85) +
  geom_vline(xintercept = MIN_STRATUM_N, linetype = "dashed",
             linewidth = 0.5, colour = unname(nature_colours["vermillion"])) +
  scale_x_log10(labels = label_number(accuracy = 1)) +
  facet_wrap(~ specification_label, scales = "free_y") +
  labs(
    title = "Distribution of intersectional stratum sizes",
    subtitle = paste0(
      "Dashed line marks the minimum cell size threshold (n = ",
      MIN_STRATUM_N, "). Horizontal axis is on a log scale."
    ),
    x = "Stratum size (n)", y = "Number of strata"
  ) + publication_theme

save_plot_pair(p_stratum_sizes,
               file.path(OUTPUT_ROOT, "figures", "individual",
                         "descriptive__stratum_size_distribution"),
               width = 10, height = 7)

# Crude event rates give the reader the marginal picture the MAIHDA model is
# then decomposing. Restricted to the composite outcomes so the panel stays
# readable; the full set is in descriptive_event_rates_by_category.csv.
headline_outcomes <- intersect(
  c("event_365d_any", "event_90d_any"), names(OUTCOME_LABELS)
)
if (length(headline_outcomes) == 0L) {
  headline_outcomes <- names(OUTCOME_LABELS)[1L]
}

p_event_rates <- ggplot(
  key_levels(descriptive_event_rates %>%
               filter(outcome %in% headline_outcomes)),
  aes(x = event_percent, y = level_key, colour = outcome_label)
) +
  geom_errorbar(aes(xmin = ci_low_percent, xmax = ci_high_percent),
                orientation = "y", width = 0, linewidth = 0.45,
                position = position_dodge(width = 0.5)) +
  geom_point(size = 2.2, position = position_dodge(width = 0.5)) +
  facet_wrap(~ variable_label, scales = "free_y", ncol = 2) +
  scale_y_discrete(labels = strip_level_key) +
  scale_colour_manual(values = unname(
    nature_colours[c("vermillion", "blue", "teal", "lavender")]
  )[seq_along(headline_outcomes)]) +
  labs(
    title = "Crude event rates by category",
    subtitle = "Unadjusted marginal rates with Wilson score 95% intervals",
    x = "Event rate (%)", y = NULL, colour = NULL
  ) + publication_theme

save_plot_pair(p_event_rates,
               file.path(OUTPUT_ROOT, "figures", "individual",
                         "descriptive__event_rates_by_category"),
               width = 11, height = 9)

descriptive_multipanel <- (p_cohort | p_event_rates) / (p_stratum_sizes |
                                                          p_retention)
save_plot_pair(descriptive_multipanel,
               file.path(OUTPUT_ROOT, "figures", "multipanel",
                         "descriptives__summary"),
               width = 20, height = 17)

if (nrow(sensitivity_table) > 0L) {
  sensitivity_long <- sensitivity_table %>%
    filter(status == "fitted") %>%
    select(specification_label, min_stratum_n,
           `VPC A (%)` = vpc_model_a_percent,
           `VPC B (%)` = vpc_model_b_percent,
           `PCV (%)` = pcv_percent,
           `Individuals retained (%)` = individuals_retained_percent) %>%
    pivot_longer(cols = -c(specification_label, min_stratum_n),
                 names_to = "quantity", values_to = "value")

  if (nrow(sensitivity_long) > 0L) {
    p_sensitivity <- ggplot(
      sensitivity_long,
      aes(x = min_stratum_n, y = value, colour = specification_label)
    ) +
      geom_line(linewidth = 0.5) +
      geom_point(size = 1.8) +
      geom_vline(xintercept = MIN_STRATUM_N, linetype = "dashed",
                 linewidth = 0.4, colour = "grey45") +
      facet_wrap(~ quantity, scales = "free_y") +
      labs(
        title = "Sensitivity to the minimum cell size threshold",
        subtitle = paste0(
          "Outcome: ", unname(OUTCOME_LABELS[SENSITIVITY_OUTCOME]),
          ". The dashed line marks the threshold used in the main analysis (n >= ",
          MIN_STRATUM_N, ")."
        ),
        x = "Minimum stratum size", y = NULL, colour = NULL
      ) + publication_theme

    save_plot_pair(p_sensitivity,
                   file.path(OUTPUT_ROOT, "figures", "individual",
                             "summary__minimum_cell_sensitivity"),
                   width = 12, height = 8)
  }
}

# =============================================================================
# 13. RUN MANIFEST AND SESSION INFORMATION
# =============================================================================

elapsed_minutes <- as.numeric(
  difftime(Sys.time(), RUN_STARTED_AT, units = "mins")
)

manifest <- data.frame(
  field = c(
    "input_file", "output_root", "test_mode", "n_specifications",
    "n_outcomes", "n_analyses_attempted", "n_analyses_completed",
    "n_analyses_failed", "min_stratum_n", "min_stratum_events",
    "estimation_engine", "mcmc_chains", "mcmc_iterations", "mcmc_warmup",
    "mcmc_adapt_delta", "prior_fixed_effects", "max_rhat", "min_ess",
    "fdr_level", "nAGQ", "seed", "ggprism_available", "svg_written",
    "plot_font", "elapsed_minutes", "r_version", "lme4_version",
    "completed_at_utc"
  ),
  value = c(
    DATA_PATH, OUTPUT_ROOT, TEST_MODE, length(SPECIFICATIONS),
    length(OUTCOME_LABELS), total_analyses, completed_analyses,
    nrow(failures_table), MIN_STRATUM_N, MIN_STRATUM_EVENTS,
    ESTIMATION_ENGINE, MCMC_CHAINS, MCMC_ITERATIONS, MCMC_WARMUP,
    MCMC_ADAPT_DELTA, PRIOR_FIXED_EFFECTS, MAX_RHAT, MIN_EFFECTIVE_SAMPLE_SIZE,
    FDR_LEVEL, N_AGQ, 20260814, HAS_GGPRISM, WRITE_SVG,
    if (nzchar(PLOT_FONT)) PLOT_FONT else "device default",
    formatC(elapsed_minutes, format = "f", digits = 1),
    R.version.string, as.character(utils::packageVersion("lme4")),
    format(Sys.time(), tz = "UTC")
  )
)
write_csv_output(manifest, file.path(OUTPUT_ROOT, "tables",
                                     "run_manifest.csv"))

# The output manifest makes an incomplete run immediately obvious and gives
# anyone reviewing the results a complete index of what was produced.
existing_outputs <- unique(output_registry$paths)
output_manifest <- data.frame(
  file = sub(paste0("^", OUTPUT_ROOT, .Platform$file.sep), "", existing_outputs),
  exists = file.exists(existing_outputs),
  size_kb = round(file.size(existing_outputs) / 1024, 1)
) %>%
  arrange(file)
write.csv(output_manifest,
          file.path(OUTPUT_ROOT, "OUTPUT_MANIFEST.csv"), row.names = FALSE)

capture.output(sessionInfo(),
               file = file.path(OUTPUT_ROOT, "logs", "sessionInfo.txt"))

singular_b <- metrics_table$outcome_label[metrics_table$singular_model_b]
completion_lines <- c(
  "MAIHDA analysis completed.",
  paste("Input:", DATA_PATH),
  paste("Output:", OUTPUT_ROOT),
  paste("Specifications:", paste(names(SPECIFICATIONS), collapse = ", ")),
  paste("Outcomes:", paste(names(OUTCOME_LABELS), collapse = ", ")),
  paste("Analyses completed:", completed_analyses, "of", total_analyses),
  paste("Analyses failed:", nrow(failures_table)),
  paste0("Minimum cell rule: n >= ", MIN_STRATUM_N,
         ", events >= ", MIN_STRATUM_EVENTS),
  paste("Files written:", nrow(output_manifest)),
  paste0("Elapsed: ", formatC(elapsed_minutes, format = "f", digits = 1),
         " minutes"),
  paste("Completed UTC:", format(Sys.time(), tz = "UTC"))
)
if (length(singular_b) > 0L) {
  completion_lines <- c(
    completion_lines, "",
    "Model B variance estimated at zero for the following analyses.",
    "This is a substantive result, not an error: additive main effects",
    "account for all between-stratum variation and no residual",
    "intersectional interaction was detected.",
    paste0("  - ", singular_b)
  )
}
writeLines(completion_lines,
           con = file.path(OUTPUT_ROOT, "logs", "completion.txt"))

for (line in completion_lines) log_message(line)
if (nrow(failures_table) > 0L) {
  log_message("Some analyses failed. See failed_analyses.csv.", level = "WARN")
}
log_message("Results written to: ", OUTPUT_ROOT)
