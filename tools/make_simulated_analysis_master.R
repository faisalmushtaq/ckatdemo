# =============================================================================
# Simulated analysis_master for MAIHDA development and testing
# =============================================================================
#
# PURPOSE
# -------
# Generates a `simulated_analysis_master.rds` with the same schema that
# run_maihda_analysis.R expects, so the analysis pipeline can be exercised end
# to end without access to patient data.
#
# The simulation deliberately contains BOTH of the situations a MAIHDA run has
# to cope with, because each exercises a different code path:
#
#   * Outcomes with genuine intersectional interaction. A small number of
#     stratum combinations receive an extra log-odds shift on top of the
#     additive linear predictor. Model B should retain non-zero between-stratum
#     variance and produce distinguishable interaction residuals.
#
#   * Outcomes generated from a purely additive linear predictor. Model B has
#     nothing left to explain, so its random-intercept variance is estimated at
#     the boundary and lme4 reports a singular fit. This is the correct answer,
#     not a failure, and the pipeline must report it as such.
#
# The cohort is also made deliberately sparse in the tail: rare ethnicity and
# frailty combinations are given small marginal probabilities so that the
# six-way stratification produces many very small cells. This is what the
# minimum-cell-size rule in the analysis script has to handle.
#
# USAGE
# -----
#   Rscript tools/make_simulated_analysis_master.R [output_path] [n_rows]
#
# Defaults to simulated_analysis_master.rds in the working directory with
# 60,000 rows, which is large enough to be realistic and small enough that a
# full 12-outcome run finishes in a few minutes.
# =============================================================================

suppressPackageStartupMessages(library(dplyr))

arguments <- commandArgs(trailingOnly = TRUE)
output_path <- if (length(arguments) >= 1L) arguments[1L] else
  "simulated_analysis_master.rds"
n_rows <- if (length(arguments) >= 2L) as.integer(arguments[2L]) else 60000L

set.seed(20260814)

# -----------------------------------------------------------------------------
# 1. Marginal distributions
# -----------------------------------------------------------------------------
# Probabilities are chosen to mimic an older UK primary-care opioid cohort:
# skewed towards younger-old age bands, majority White ethnicity, and a frailty
# distribution weighted towards fit and mildly frail.

draw <- function(values, probabilities) {
  factor(sample(values, n_rows, replace = TRUE, prob = probabilities),
         levels = values)
}

age_band <- draw(
  c("65-69", "70-74", "75-79", "80-84", "85+"),
  c(0.26, 0.24, 0.20, 0.16, 0.14)
)
sex <- draw(c("Male", "Female"), c(0.44, 0.56))
ethnicity <- draw(
  c("White", "Asian", "Black", "Mixed", "Other"),
  c(0.855, 0.070, 0.035, 0.018, 0.022)
)
imd_quintile <- draw(
  c("1 (Least deprived)", "2", "3", "4", "5 (Most deprived)"),
  c(0.18, 0.19, 0.20, 0.21, 0.22)
)
opioid_strength <- draw(c("Weak", "Moderate", "Strong"), c(0.52, 0.30, 0.18))
efi_category <- draw(
  c("Fit", "Mild", "Moderate", "Severe"),
  c(0.30, 0.36, 0.24, 0.10)
)

analysis_master <- tibble(
  patient_id = sprintf("SIM%07d", seq_len(n_rows)),
  age_band = age_band,
  sex = sex,
  ethnicity = ethnicity,
  imd_quintile = imd_quintile,
  opioid_strength = opioid_strength,
  efi_category = efi_category
)

# -----------------------------------------------------------------------------
# 2. Collapsed axes
# -----------------------------------------------------------------------------
# These mirror the collapsed variables the analysis script expects, including
# the exact reference-level wording.

analysis_master <- analysis_master %>%
  mutate(
    age_3cat = factor(
      case_when(
        age_band %in% c("65-69", "70-74") ~ "65-74",
        age_band %in% c("75-79", "80-84") ~ "75-84",
        TRUE ~ "85+"
      ),
      levels = c("65-74", "75-84", "85+")
    ),
    ethnicity_4cat = factor(
      case_when(
        ethnicity == "White" ~ "White",
        ethnicity == "Asian" ~ "Asian",
        ethnicity == "Black" ~ "Black",
        TRUE ~ "Other"
      ),
      levels = c("White", "Asian", "Black", "Other")
    ),
    imd_3cat = factor(
      case_when(
        imd_quintile %in% c("1 (Least deprived)", "2") ~ "Q1-Q2",
        imd_quintile == "3" ~ "Q3",
        TRUE ~ "Q4-Q5"
      ),
      levels = c("Q1-Q2", "Q3", "Q4-Q5")
    ),
    efi_category_3 = factor(
      case_when(
        efi_category == "Fit" ~ "Fit",
        efi_category == "Mild" ~ "Mild",
        TRUE ~ "Moderate-Severe"
      ),
      levels = c("Fit", "Mild", "Moderate-Severe")
    )
  )

# -----------------------------------------------------------------------------
# 3. Additive linear predictor
# -----------------------------------------------------------------------------
# Effect sizes are on the log-odds scale and are broadly consistent with the
# direction and magnitude reported in opioid safety literature: risk rises with
# age and frailty, is modestly higher in women and in more deprived quintiles,
# and rises with opioid strength.

effect_of <- function(variable, effects) unname(effects[as.character(variable)])

additive_linear_predictor <- with(analysis_master,
  effect_of(age_band, c(`65-69` = 0.00, `70-74` = 0.18, `75-79` = 0.42,
                        `80-84` = 0.70, `85+` = 1.05)) +
  effect_of(sex, c(Male = 0.00, Female = 0.16)) +
  effect_of(ethnicity, c(White = 0.00, Asian = -0.14, Black = -0.08,
                         Mixed = 0.05, Other = -0.05)) +
  effect_of(imd_quintile, c(`1 (Least deprived)` = 0.00, `2` = 0.09,
                            `3` = 0.17, `4` = 0.27, `5 (Most deprived)` = 0.38)) +
  effect_of(opioid_strength, c(Weak = 0.00, Moderate = 0.30, Strong = 0.62)) +
  effect_of(efi_category, c(Fit = 0.00, Mild = 0.34, Moderate = 0.72,
                            Severe = 1.15))
)

# -----------------------------------------------------------------------------
# 4. True intersectional interaction
# -----------------------------------------------------------------------------
# Genuine MAIHDA interactions are sparse: most strata sit on the additive
# surface and a handful depart from it. Three departures are injected here,
# each chosen to be substantively plausible and each large enough to be
# detectable at this sample size.

interaction_linear_predictor <- with(analysis_master,
  ifelse(age_band == "85+" & efi_category == "Severe" &
           opioid_strength == "Strong", 0.85, 0) +
  ifelse(sex == "Female" & imd_quintile == "5 (Most deprived)" &
           efi_category %in% c("Moderate", "Severe"), 0.45, 0) +
  ifelse(ethnicity == "Asian" & age_band %in% c("80-84", "85+"), -0.40, 0)
)

# A small amount of unstructured stratum-level noise keeps the between-stratum
# variance from being exactly the sum of the three injected terms, which is
# more realistic and gives the caterpillar plots a natural spread.
stratum_key <- interaction(
  analysis_master$age_band, analysis_master$sex, analysis_master$ethnicity,
  analysis_master$imd_quintile, analysis_master$opioid_strength,
  analysis_master$efi_category, drop = TRUE
)
stratum_noise <- rnorm(nlevels(stratum_key), mean = 0, sd = 0.16)
interaction_linear_predictor <- interaction_linear_predictor +
  stratum_noise[as.integer(stratum_key)]

# -----------------------------------------------------------------------------
# 5. Outcomes
# -----------------------------------------------------------------------------
# Intercepts set the marginal event rate for each outcome. Composite outcomes
# are built from their components so that `any` is genuinely the union of the
# component events, exactly as in the real analysis_master.

simulate_binary <- function(intercept, include_interaction) {
  eta <- intercept + additive_linear_predictor +
    if (include_interaction) interaction_linear_predictor else 0
  rbinom(n_rows, size = 1L, prob = plogis(eta))
}

# Falls and delirium carry the injected interactions. Fractures and death are
# generated additively so that the pipeline is also tested on outcomes whose
# Model B is legitimately singular.
analysis_master <- analysis_master %>%
  mutate(
    event_90d_fall      = simulate_binary(-3.10, TRUE),
    event_365d_fall     = simulate_binary(-1.95, TRUE),
    event_90d_delirium  = simulate_binary(-4.05, TRUE),
    event_365d_delirium = simulate_binary(-2.90, TRUE),
    event_90d_fracture  = simulate_binary(-4.30, FALSE),
    event_365d_fracture = simulate_binary(-3.25, FALSE),
    event_90d_death     = simulate_binary(-4.60, FALSE),
    event_365d_death    = simulate_binary(-3.40, FALSE)
  ) %>%
  mutate(
    event_90d_nonfatal = pmax(event_90d_fall, event_90d_fracture,
                              event_90d_delirium),
    event_365d_nonfatal = pmax(event_365d_fall, event_365d_fracture,
                               event_365d_delirium),
    event_90d_any = pmax(event_90d_nonfatal, event_90d_death),
    event_365d_any = pmax(event_365d_nonfatal, event_365d_death)
  )

# -----------------------------------------------------------------------------
# 6. Realistic missingness
# -----------------------------------------------------------------------------
# Ethnicity and deprivation are the variables most often incomplete in primary
# care extracts. A small amount of missingness exercises the complete-case
# handling and the exclusion accounting in the analysis script.

introduce_missing <- function(x, proportion) {
  x[sample.int(n_rows, floor(proportion * n_rows))] <- NA
  x
}
analysis_master$ethnicity <- introduce_missing(analysis_master$ethnicity, 0.018)
analysis_master$ethnicity_4cat <- ifelse(
  is.na(analysis_master$ethnicity), NA,
  as.character(analysis_master$ethnicity_4cat)
)
analysis_master$ethnicity_4cat <- factor(
  analysis_master$ethnicity_4cat, levels = c("White", "Asian", "Black", "Other")
)
analysis_master$imd_quintile <- introduce_missing(
  analysis_master$imd_quintile, 0.009
)
analysis_master$imd_3cat <- ifelse(
  is.na(analysis_master$imd_quintile), NA,
  as.character(analysis_master$imd_3cat)
)
analysis_master$imd_3cat <- factor(
  analysis_master$imd_3cat, levels = c("Q1-Q2", "Q3", "Q4-Q5")
)

saveRDS(analysis_master, output_path)

outcome_names <- grep("^event_", names(analysis_master), value = TRUE)
message("Wrote ", output_path)
message("Rows: ", nrow(analysis_master), "; columns: ", ncol(analysis_master))
message("Observed six-way strata: ", nlevels(stratum_key))
for (outcome in outcome_names) {
  message(sprintf("  %-22s events = %6d (%.2f%%)", outcome,
                  sum(analysis_master[[outcome]]),
                  100 * mean(analysis_master[[outcome]])))
}
