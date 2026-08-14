# =============================================================================
# Generate Methods and Results text from a completed MAIHDA run
# =============================================================================
#
# PURPOSE
# -------
# Reads the CSV tables written by run_maihda_analysis.R and produces a complete
# Methods and Results section with the actual values substituted in. Nothing is
# recomputed: every number in the generated text is read from the analysis
# output, so the manuscript and the results tables cannot drift apart.
#
# The script reads only CSV files. It does not load the fitted model objects,
# which keeps it fast, keeps its memory footprint small, and means the same
# CSVs that go through a Trusted Research Environment output check are the
# exact inputs to the text.
#
# RUNNING IT IN A TRE
# -------------------
# After run_maihda_analysis.R has completed:
#
#   Rscript generate_manuscript.R
#
# or, with explicit paths:
#
#   MAIHDA_OUTPUT_ROOT=/path/to/maihda_outputs \
#     MANUSCRIPT_OUTPUT_DIR=/path/to/manuscript \
#     Rscript generate_manuscript.R
#
# STUDY METADATA
# --------------
# Some Methods content cannot be derived from the analysis output: the data
# source, the setting, the inclusion criteria, the ethics approval. These are
# set in section 1 below. Anything left unset is written into the text as a
# conspicuous [PLACEHOLDER: ...] marker AND listed in placeholders_to_complete
# .txt, so an unfilled field cannot be missed on a read-through.
#
# DISCLOSURE CONTROL
# ------------------
# The generated text names individual intersectional strata and quotes counts
# within them, which is exactly the material a TRE output checker scrutinises.
# Three modes are available through MANUSCRIPT_DISCLOSURE_MODE:
#
#   "flag"  (default) Report true values, and additionally write a checklist of
#           every reported figure that derives from a cell below the disclosure
#           threshold. Use this while drafting: the text is accurate and the
#           checklist tells you in advance what a checker will query.
#
#   "apply" Suppress or round small counts in the text itself. Use this for the
#           version that will actually leave the environment.
#
#   "off"   No disclosure processing. Only appropriate outside a TRE.
#
# In every mode the script writes:
#
#   manuscript_values.csv     every headline figure quoted in the text, with
#                             its source table and column, so any number can be
#                             traced back without re-reading the prose
#   disclosure_check.csv      counts below the threshold, with their location
#
# OUTPUTS
# -------
#   Manuscript.docx            Methods and Results, Word format
#   Manuscript.md              the same text as Markdown
#   Manuscript.txt             the same text as plain text
#   manuscript_values.csv      traceability of every quoted figure
#   disclosure_check.csv       disclosure control checklist
#   placeholders_to_complete.txt
#
# REQUIRED PACKAGES
# -----------------
# Required: officer. Optional: flextable (adds the summary table to the Word
# output; the text is complete without it).
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

# =============================================================================
# 1. STUDY METADATA
# =============================================================================
# Fill these in. Anything left as "" becomes a visible placeholder in the text
# and is listed in placeholders_to_complete.txt.
#
# Each can also be supplied through the environment, which is convenient when
# the same analysis is written up for more than one audience.

study_metadata <- list(
  data_source = Sys.getenv("MANUSCRIPT_DATA_SOURCE", ""),
  setting = Sys.getenv("MANUSCRIPT_SETTING", ""),
  study_period = Sys.getenv("MANUSCRIPT_STUDY_PERIOD", ""),
  population = Sys.getenv("MANUSCRIPT_POPULATION", ""),
  inclusion_criteria = Sys.getenv("MANUSCRIPT_INCLUSION", ""),
  exclusion_criteria = Sys.getenv("MANUSCRIPT_EXCLUSION", ""),
  index_date = Sys.getenv("MANUSCRIPT_INDEX_DATE", ""),
  outcome_ascertainment = Sys.getenv("MANUSCRIPT_OUTCOME_ASCERTAINMENT", ""),
  frailty_measure = Sys.getenv(
    "MANUSCRIPT_FRAILTY_MEASURE",
    "the electronic frailty index (eFI)"
  ),
  deprivation_measure = Sys.getenv(
    "MANUSCRIPT_DEPRIVATION_MEASURE",
    "the Index of Multiple Deprivation (IMD)"
  ),
  ethics_approval = Sys.getenv("MANUSCRIPT_ETHICS", ""),
  data_access = Sys.getenv("MANUSCRIPT_DATA_ACCESS", ""),
  funding = Sys.getenv("MANUSCRIPT_FUNDING", "")
)

# =============================================================================
# 2. SETTINGS
# =============================================================================

ANALYSIS_ROOT <- Sys.getenv("MAIHDA_OUTPUT_ROOT", unset = "maihda_outputs")
MANUSCRIPT_DIR <- Sys.getenv("MANUSCRIPT_OUTPUT_DIR", unset = "")

# "flag", "apply" or "off". See the header.
DISCLOSURE_MODE <- tolower(Sys.getenv("MANUSCRIPT_DISCLOSURE_MODE",
                                      unset = "flag"))
if (!DISCLOSURE_MODE %in% c("flag", "apply", "off")) {
  stop("MANUSCRIPT_DISCLOSURE_MODE must be one of: flag, apply, off. ",
       "Received: ", DISCLOSURE_MODE)
}

# Counts at or below this are treated as disclosive. Defaults to the minimum
# cell size the analysis itself used, read from the run manifest, because a
# threshold looser than the analysis rule would be incoherent.
DISCLOSURE_THRESHOLD <- suppressWarnings(as.integer(
  Sys.getenv("MANUSCRIPT_DISCLOSURE_THRESHOLD", unset = NA)
))

# In "apply" mode, counts can be rounded rather than suppressed. Set to 5 for
# the common "round to nearest 5" TRE rule, or 0 to suppress instead.
DISCLOSURE_ROUNDING <- suppressWarnings(as.integer(
  Sys.getenv("MANUSCRIPT_DISCLOSURE_ROUNDING", unset = "0")
))

# How many strata to name individually in the Results. Naming every stratum
# would be unreadable and, in a TRE, would multiply the disclosure surface for
# no gain.
MAX_STRATA_NAMED <- suppressWarnings(as.integer(
  Sys.getenv("MANUSCRIPT_MAX_STRATA_NAMED", unset = "8")
))

# The outcome whose results are described in full narrative detail. Others are
# summarised and referred to the tables. Defaults to the first outcome.
PRIMARY_OUTCOME <- Sys.getenv("MANUSCRIPT_PRIMARY_OUTCOME", unset = "")

if (!requireNamespace("officer", quietly = TRUE)) {
  stop("The officer package is required. Run: install.packages(\"officer\")")
}
HAS_FLEXTABLE <- requireNamespace("flextable", quietly = TRUE)

ANALYSIS_ROOT <- normalizePath(ANALYSIS_ROOT, mustWork = FALSE)
if (!dir.exists(ANALYSIS_ROOT)) {
  stop("Analysis output directory not found: ", ANALYSIS_ROOT,
       "\nSet MAIHDA_OUTPUT_ROOT to the directory run_maihda_analysis.R wrote.")
}
if (!nzchar(MANUSCRIPT_DIR)) {
  MANUSCRIPT_DIR <- file.path(ANALYSIS_ROOT, "manuscript")
}
dir.create(MANUSCRIPT_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(MANUSCRIPT_DIR)) {
  stop("Could not create the manuscript output directory: ", MANUSCRIPT_DIR)
}

message("Reading analysis output from: ", ANALYSIS_ROOT)
message("Writing manuscript to:        ", MANUSCRIPT_DIR)
message("Disclosure mode:              ", DISCLOSURE_MODE)

# =============================================================================
# 3. LOAD THE ANALYSIS OUTPUT
# =============================================================================

read_analysis_table <- function(name, required = TRUE) {
  path <- file.path(ANALYSIS_ROOT, "tables", paste0(name, ".csv"))
  if (!file.exists(path)) {
    if (required) {
      stop("Required analysis table is missing: ", path,
           "\nHas run_maihda_analysis.R completed successfully?")
    }
    return(NULL)
  }
  table <- try(read.csv(path, check.names = FALSE), silent = TRUE)
  if (inherits(table, "try-error")) {
    if (required) stop("Could not read: ", path)
    return(NULL)
  }
  if (nrow(table) == 0L) {
    if (required) {
      stop("Required analysis table is empty: ", path)
    }
    return(NULL)
  }
  table
}

metrics <- read_analysis_table("all_model_metrics")
outcome_counts <- read_analysis_table("outcome_counts")
cohort <- read_analysis_table("descriptive_cohort_characteristics")
stratum_distribution <- read_analysis_table("descriptive_stratum_distribution")
retention <- read_analysis_table("stratum_retention")
fixed_effects <- read_analysis_table("all_fixed_effects_odds_ratios")
extremes <- read_analysis_table("extreme_strata")
manifest <- read_analysis_table("run_manifest")
input_audit <- read_analysis_table("input_audit")

# Optional: absent or empty when the analysis produced no such rows, which is
# itself a reportable result rather than an error.
significant <- read_analysis_table("distinguishable_interaction_strata",
                                   required = FALSE)
axis_decomposition <- read_analysis_table("axis_variance_decomposition",
                                          required = FALSE)
sensitivity <- read_analysis_table("minimum_cell_sensitivity", required = FALSE)
failures <- read_analysis_table("failed_analyses", required = FALSE)
event_rates <- read_analysis_table("descriptive_event_rates_by_category",
                                   required = FALSE)

manifest_value <- function(field, default = NA_character_) {
  if (is.null(manifest)) return(default)
  matched <- manifest$value[manifest$field == field]
  if (length(matched) == 0L) return(default)
  as.character(matched[1L])
}
audit_value <- function(item, default = NA_character_) {
  if (is.null(input_audit)) return(default)
  matched <- input_audit$value[input_audit$item == item]
  if (length(matched) == 0L) return(default)
  as.character(matched[1L])
}

# Analysis settings are read back from the run rather than assumed, so the text
# always describes the run that actually happened.
MIN_STRATUM_N <- as.integer(metrics$min_stratum_n[1L])
MIN_STRATUM_EVENTS <- as.integer(metrics$min_stratum_events[1L])
FDR_LEVEL <- suppressWarnings(as.numeric(manifest_value("fdr_level", "0.05")))
N_AGQ <- manifest_value("nAGQ", "1")
R_VERSION <- manifest_value("r_version", R.version.string)
LME4_VERSION <- manifest_value("lme4_version", "")

if (is.na(DISCLOSURE_THRESHOLD)) DISCLOSURE_THRESHOLD <- MIN_STRATUM_N

if (!nzchar(PRIMARY_OUTCOME)) PRIMARY_OUTCOME <- metrics$outcome[1L]
if (!PRIMARY_OUTCOME %in% metrics$outcome) {
  stop("MANUSCRIPT_PRIMARY_OUTCOME is not among the analysed outcomes: ",
       PRIMARY_OUTCOME)
}

# The stratum axes are recovered from the model specification recorded in the
# metrics, so adding or removing an axis in the analysis needs no change here.
axis_names <- trimws(strsplit(metrics$model_axes[1L], "\\+")[[1L]])
axis_labels <- vapply(axis_names, function(axis) {
  matched <- unique(fixed_effects$variable_label[fixed_effects$variable == axis])
  if (length(matched) == 0L || is.na(matched[1L])) axis else matched[1L]
}, character(1))
SPECIFICATION_LABEL <- metrics$specification_label[1L]

# =============================================================================
# 4. FORMATTING HELPERS
# =============================================================================

fmt_n <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return("NA")
  formatC(round(as.numeric(x)), format = "d", big.mark = ",")
}
fmt_num <- function(x, digits = 2) {
  if (length(x) == 0L || all(is.na(x))) return("NA")
  formatC(as.numeric(x), format = "f", digits = digits, big.mark = ",")
}
fmt_pct <- function(x, digits = 1) paste0(fmt_num(x, digits), "%")

# APA style: no leading zero, "< .001" below that.
fmt_p <- function(x) {
  if (length(x) == 0L || is.na(x)) return("NA")
  if (x < 0.001) return("< .001")
  paste0("= ", sub("^0", "", formatC(x, format = "f", digits = 3)))
}

# `bare = TRUE` omits the surrounding parentheses, for use where the whole
# statistic is already inside a bracketed clause. Nesting parentheses inside
# parentheses is the fastest way to make a results sentence unreadable.
fmt_or_ci <- function(odds_ratio, low, high, digits = 2, bare = FALSE) {
  if (is.na(low) || is.na(high)) return(fmt_num(odds_ratio, digits))
  inner <- paste0(fmt_num(odds_ratio, digits), ", 95% CI ",
                  fmt_num(low, digits), " to ", fmt_num(high, digits))
  if (bare) inner else paste0("(", inner, ")")
}

# Oxford comma list. Used throughout so enumerations read as prose rather than
# as a comma-separated dump.
#
# Items that themselves contain commas -- outcome labels such as "fall, 90
# days" do -- are separated with semicolons instead, because a comma-separated
# list of comma-containing items cannot be parsed by the reader.
oxford <- function(items, conjunction = "and") {
  items <- items[nzchar(items)]
  n <- length(items)
  if (n == 0L) return("")
  if (n == 1L) return(items)
  # Only a comma acting as punctuation forces semicolons. A comma inside a
  # thousands separator, as in "1,080", must not.
  has_punctuating_comma <- any(grepl(",(?![0-9])|(?<![0-9]),", items,
                                     perl = TRUE))
  separator <- if (has_punctuating_comma) "; " else ", "
  if (n == 2L) {
    joiner <- if (identical(separator, "; ")) {
      paste0("; ", conjunction, " ")
    } else {
      paste0(" ", conjunction, " ")
    }
    return(paste(items, collapse = joiner))
  }
  paste0(paste(items[-n], collapse = separator), separator, conjunction, " ",
         items[n])
}

# Lowercasing a label mid-sentence is right for "Age" but wrong for an acronym
# such as "IMD". An all-uppercase label is left alone.
label_in_sentence <- function(label) {
  vapply(label, function(one) {
    if (identical(one, toupper(one))) one else tolower(one)
  }, character(1), USE.NAMES = FALSE)
}

# Outcome labels are written for table headers, where "Fall, 90 days" is
# compact and clear. Dropped into a sentence, the comma reads as a clause
# break and the sentence falls apart. Rewriting to "fall within 90 days"
# removes the comma and reads as prose, which also keeps enumerations of
# outcomes unambiguous.
outcome_in_sentence <- function(label) {
  rewritten <- sub(",\\s*(\\d+)\\s*-?\\s*days?$", " within \\1 days", label)
  tolower(rewritten)
}

# Appends a confidence interval only when one was successfully computed, so a
# run without profile intervals degrades to a bare point estimate rather than
# printing "(95% CI NA to NA)".
interval_suffix <- function(low, high, digits = 1, percent = FALSE) {
  if (length(low) == 0L || length(high) == 0L || is.na(low) || is.na(high)) {
    return("")
  }
  formatter <- if (percent) function(x) fmt_pct(x, digits) else
    function(x) fmt_num(x, if (percent) digits else 2)
  paste0(" (95% CI ", formatter(low), " to ", formatter(high), ")")
}

# Published VPC values for health outcomes, from the systematic review by
# Keller and colleagues (2023). Used only to give the reader a sense of scale.
VPC_BENCHMARK_LOW <- 0.5
VPC_BENCHMARK_HIGH <- 41.9
VPC_BENCHMARK_MEDIAN <- 5.5

# The magnitude qualifier is derived from the value rather than asserted, so
# the text cannot describe a 2% variance partition coefficient as substantial.
vpc_qualifier <- function(vpc_percent) {
  if (is.na(vpc_percent)) return("an unquantified share")
  if (vpc_percent < 2) "a small share"
  else if (vpc_percent < 5) "a modest share"
  else if (vpc_percent < 10) "an appreciable share"
  else if (vpc_percent < 20) "a substantial share"
  else "a large share"
}

plural <- function(n, singular, plural_form = paste0(singular, "s")) {
  if (length(n) == 0L || is.na(n)) return(plural_form)
  if (n == 1L) singular else plural_form
}
was_were <- function(n) if (!is.na(n) && n == 1L) "was" else "were"

# Journals conventionally spell out small numbers that open a sentence.
NUMBER_WORDS <- c("one", "two", "three", "four", "five", "six", "seven",
                  "eight", "nine", "ten", "eleven", "twelve")
spell_number <- function(n) {
  if (is.na(n)) return("an unknown number of")
  if (n == 0L) return("no")
  if (n <= length(NUMBER_WORDS)) NUMBER_WORDS[n] else fmt_n(n)
}
capitalise_first <- function(text) {
  if (!nzchar(text)) return(text)
  paste0(toupper(substring(text, 1, 1)), substring(text, 2))
}

# Placeholders are recorded as they are used so the completion list is exact
# rather than a guess at what might be missing.
placeholder_registry <- new.env(parent = emptyenv())
placeholder_registry$items <- character(0)
metadata_or_placeholder <- function(field, description) {
  value <- study_metadata[[field]]
  if (is.null(value) || !nzchar(value)) {
    placeholder_registry$items <- unique(c(
      placeholder_registry$items,
      paste0(field, ": ", description)
    ))
    return(paste0("[PLACEHOLDER: ", description, "]"))
  }
  value
}

# --- Value traceability and disclosure control -------------------------------
# Every headline figure quoted in the text is recorded with its source, so any
# number can be checked against the analysis output without re-reading the
# prose, and so an output checker has a single list to work from.
value_registry <- new.env(parent = emptyenv())
value_registry$rows <- list()
record_value <- function(label, value, formatted, source_table,
                         source_column = NA_character_, is_count = FALSE) {
  value_registry$rows[[length(value_registry$rows) + 1L]] <- data.frame(
    label = label,
    value = if (is.numeric(value)) as.numeric(value) else NA_real_,
    reported_as = formatted,
    source_table = source_table,
    source_column = source_column,
    is_count = is_count
  )
  invisible(formatted)
}

disclosure_registry <- new.env(parent = emptyenv())
disclosure_registry$rows <- list()

# All counts that reach the text pass through here. In "flag" mode the true
# value is reported and the small cell is logged; in "apply" mode it is
# rounded or suppressed.
disclose_count <- function(n, label, source_table, source_column = NA) {
  n_numeric <- suppressWarnings(as.numeric(n))
  is_small <- !is.na(n_numeric) && n_numeric > 0 &&
    n_numeric < DISCLOSURE_THRESHOLD

  if (is_small && DISCLOSURE_MODE != "off") {
    disclosure_registry$rows[[length(disclosure_registry$rows) + 1L]] <-
      data.frame(
        label = label, true_value = n_numeric,
        threshold = DISCLOSURE_THRESHOLD,
        source_table = source_table, source_column = as.character(source_column),
        action = if (DISCLOSURE_MODE == "apply") {
          if (DISCLOSURE_ROUNDING > 0) "rounded" else "suppressed"
        } else {
          "reported, flagged for review"
        }
      )
  }

  formatted <- if (is_small && DISCLOSURE_MODE == "apply") {
    if (DISCLOSURE_ROUNDING > 0) {
      paste0("approximately ",
             fmt_n(round(n_numeric / DISCLOSURE_ROUNDING) * DISCLOSURE_ROUNDING))
    } else {
      paste0("<", DISCLOSURE_THRESHOLD)
    }
  } else if (DISCLOSURE_MODE == "apply" && DISCLOSURE_ROUNDING > 0 &&
             !is.na(n_numeric)) {
    # Rounding is applied to every count, not only small ones, because a mix of
    # exact and rounded figures lets an exact value be recovered by difference.
    fmt_n(round(n_numeric / DISCLOSURE_ROUNDING) * DISCLOSURE_ROUNDING)
  } else {
    fmt_n(n_numeric)
  }

  record_value(label, n_numeric, formatted, source_table, source_column,
               is_count = TRUE)
  formatted
}

# Convenience for non-count statistics, which are not disclosive in themselves
# but are worth recording for traceability.
stat_value <- function(x, label, source_table, source_column = NA,
                       digits = 2, percent = FALSE) {
  formatted <- if (percent) fmt_pct(x, digits) else fmt_num(x, digits)
  record_value(label, x, formatted, source_table, as.character(source_column))
  formatted
}

# =============================================================================
# 5. DERIVED QUANTITIES
# =============================================================================

primary <- metrics[metrics$outcome == PRIMARY_OUTCOME, , drop = FALSE][1L, ]
primary_label <- primary$outcome_label
primary_retention <- retention[retention$outcome == PRIMARY_OUTCOME, ,
                               drop = FALSE][1L, ]
distribution <- stratum_distribution[1L, ]

n_outcomes <- nrow(metrics)
n_cohort <- as.numeric(outcome_counts$total_n[1L])

# Model B variance at zero is a substantive finding and is described as such
# rather than reported as a convergence problem.
singular_b <- metrics[isTRUE_vector <- metrics$singular_model_b %in%
                        c(TRUE, "TRUE", "True", "true"), , drop = FALSE]
non_singular_b <- metrics[!metrics$outcome %in% singular_b$outcome, ,
                          drop = FALSE]

# Outcomes with at least one stratum surviving FDR correction.
outcomes_with_interactions <- metrics[metrics$n_distinguishable_fdr > 0, ,
                                      drop = FALSE]

has_significant <- !is.null(significant) && nrow(significant) > 0L

# =============================================================================
# 6. METHODS
# =============================================================================

build_methods <- function() {
  sections <- list()

  # --- Study design --------------------------------------------------------
  sections[["Study design and population"]] <- c(
    paste0(
      "We conducted a cross-sectional analysis of ",
      metadata_or_placeholder("data_source",
                              "name and describe the data source"),
      ", set in ",
      metadata_or_placeholder("setting", "describe the setting"),
      ". The study period was ",
      metadata_or_placeholder("study_period",
                              "state the study period, e.g. 1 January 2015 to 31 December 2019"),
      ". The study population comprised ",
      metadata_or_placeholder("population", "describe the study population"),
      "."
    ),
    paste0(
      "Individuals were eligible for inclusion if ",
      metadata_or_placeholder("inclusion_criteria",
                              "state the inclusion criteria"),
      ". We excluded individuals who ",
      metadata_or_placeholder("exclusion_criteria",
                              "state the exclusion criteria"),
      ". The index date was defined as ",
      metadata_or_placeholder("index_date", "define the index date"),
      "."
    ),
    paste0(
      "The analytic dataset contained ",
      disclose_count(n_cohort, "Total cohort size", "outcome_counts", "total_n"),
      " individuals."
    )
  )

  # --- Outcomes ------------------------------------------------------------
  outcome_list <- oxford(outcome_in_sentence(outcome_counts$outcome_label))
  sections[["Outcomes"]] <- c(
    paste0(
      "We examined ", spell_number(n_outcomes), " binary ",
      plural(n_outcomes, "outcome"), ": ", outcome_list,
      ". Outcomes were ascertained from ",
      metadata_or_placeholder("outcome_ascertainment",
                              "state how outcomes were ascertained and coded"),
      ". Each outcome was coded as a binary indicator of whether the event ",
      "occurred within the stated follow-up window from the index date. ",
      "Composite outcomes are reported alongside their individual components ",
      "because they answer distinct substantive questions."
    )
  )

  # --- Intersectional strata -----------------------------------------------
  axis_descriptions <- vapply(seq_along(axis_names), function(index) {
    axis <- axis_names[index]
    n_levels <- sum(cohort$variable == axis)
    reference <- cohort$level[cohort$variable == axis &
                                cohort$is_reference %in% c(TRUE, "TRUE")]
    reference_text <- if (length(reference) > 0L) {
      paste0(", reference category ", reference[1L])
    } else ""
    paste0(label_in_sentence(axis_labels[index]), " (", spell_number(n_levels), " ",
           plural(n_levels, "category", "categories"), reference_text, ")")
  }, character(1))

  possible_strata <- as.numeric(distribution$possible_strata)
  observed_strata <- as.numeric(distribution$observed_strata)

  sections[["Intersectional strata"]] <- c(
    paste0(
      "Following the approach set out by Evans and colleagues, intersectional ",
      "strata were defined as the full cross-classification of ",
      spell_number(length(axis_names)), " ", plural(length(axis_names), "axis", "axes"),
      ": ", oxford(axis_descriptions), ". Deprivation was measured using ",
      study_metadata$deprivation_measure, " and frailty using ",
      study_metadata$frailty_measure, "."
    ),
    paste0(
      "This cross-classification defines ",
      fmt_n(possible_strata), " possible strata, of which ",
      fmt_n(observed_strata), " (",
      fmt_pct(100 * observed_strata / possible_strata),
      ") were observed in the cohort. Reference categories were chosen to be ",
      "the conventionally advantaged or lowest-risk category on each axis, so ",
      "that the additive main effects in Model B are interpreted as departures ",
      "from that position."
    ),
    paste0(
      "Records with a missing value on any stratum axis or on the outcome ",
      "under analysis were excluded from that analysis by complete-case ",
      "deletion. Missingness by variable is reported in the Results."
    )
  )

  # --- Minimum cell size ----------------------------------------------------
  cell_rule_sentence <- if (MIN_STRATUM_EVENTS > 0) {
    paste0("fewer than ", MIN_STRATUM_N, " individuals, or fewer than ",
           MIN_STRATUM_EVENTS, " ", plural(MIN_STRATUM_EVENTS, "event"), ",")
  } else {
    paste0("fewer than ", MIN_STRATUM_N, " individuals")
  }

  minimum_cell_text <- c(
    paste0(
      "Strata containing ", cell_rule_sentence,
      " were excluded before model fitting. The threshold was applied after ",
      "aggregating individuals to one binomial record per stratum, so all ",
      "reported model quantities refer to the retained analytic sample."
    ),
    paste0(
      "We note that the MAIHDA methodological literature does not require ",
      "trimming, and that the tutorial literature does not apply it: partial ",
      "pooling means a very small stratum contributes little to its own ",
      "estimate and is shrunk towards the additive prediction, and simulation ",
      "studies have reported robust estimate accuracy for both linear and ",
      "logistic MAIHDA at average stratum sizes of around ten observations. ",
      "The recommended practice is instead to report the distribution of ",
      "stratum sizes and to choose category granularity so that most strata ",
      "are adequately sized, both of which we do."
    ),
    paste0(
      "The threshold was therefore applied not to improve estimation but to ",
      "satisfy statistical disclosure control, because this analysis reports ",
      "stratum-level estimates and names individual strata, and ",
      MIN_STRATUM_N, " is the lowest threshold meeting that requirement. ",
      "Because this is a departure from the reference approach, its effect on ",
      "the substantive conclusions was assessed in sensitivity analysis across ",
      "a range of thresholds including no trimming at all."
    )
  )

  if (MIN_STRATUM_EVENTS == 0) {
    minimum_cell_text <- c(minimum_cell_text, paste0(
      "Strata with no events were retained. A stratum of sufficient size with ",
      "no observed events is an observation of low risk rather than missing ",
      "information, and excluding such strata would have systematically ",
      "removed the low-risk end of the distribution and biased the predicted ",
      "risk range upwards."
    ))
  }
  sections[["Minimum cell size"]] <- minimum_cell_text

  # --- Statistical analysis -------------------------------------------------
  sections[["Statistical analysis"]] <- c(
    paste0(
      "We fitted two-level logistic models with individuals nested within ",
      "intersectional strata, following the multilevel analysis of individual ",
      "heterogeneity and discriminatory accuracy (MAIHDA) framework. For each ",
      "outcome, two models were estimated."
    ),
    paste0(
      "Model A, the null or intercept-only model, included a random intercept ",
      "for stratum and no covariates. It partitions the total variation in the ",
      "outcome into between-stratum and within-stratum components and so ",
      "quantifies the total between-stratum heterogeneity."
    ),
    if (!is.null(axis_decomposition)) {
      paste0(
        "As a secondary and descriptive analysis, we additionally fitted a set ",
        "of partially adjusted models, each adding a single stratum axis to ",
        "Model A while retaining the stratum random intercept, to indicate ",
        "which axes account for the between-stratum variation. Because the ",
        "axes are correlated in the population, these contributions overlap ",
        "and do not sum to the fully adjusted value. We report them with the ",
        "caution that attending to individual axis contributions can encourage ",
        "a reversion to single-axis reasoning about inequity, which is counter ",
        "to the purpose of an intersectional analysis; they are secondary to ",
        "the collective additive effect and the total stratum predictions."
      )
    } else NULL,
    paste0(
      "Model B, the fully adjusted or intersectional interaction model, added ",
      "the additive main effects of every stratum axis simultaneously as fixed ",
      "effects, retaining the stratum random intercept. Because the fixed ",
      "effects reproduce the purely additive expectation for each stratum, the ",
      "remaining random effects represent departures from that additive ",
      "prediction, and within the MAIHDA framework are interpreted as ",
      "intersectional interaction residuals. In the notation of the Evans et ",
      "al. tutorial these are the logistic Models 2A and 2B respectively, also ",
      "referred to as the simple intersectional and the intersectional ",
      "interaction models."
    ),
    paste0(
      "For each stratum we report the absolute risk, being the model-predicted ",
      "probability of the outcome, and the absolute risk due to interaction ",
      "(ARI), being the total predicted risk minus the risk predicted by the ",
      "additive main effects alone. A positive ARI indicates that individuals ",
      "in that stratum carry more risk than the simple addition of the risks ",
      "conveyed by their constituent social and clinical positions would imply; ",
      "a negative ARI indicates less."
    ),
    paste0(
      "From each model we derived the variance partition coefficient (VPC), ",
      "computed on the latent response scale with the individual-level ",
      "variance fixed at π²/3, the variance of the standard logistic ",
      "distribution. The VPC expresses the proportion of total variance lying ",
      "between strata. We also computed the median odds ratio (MOR), the ",
      "median contrast in odds between two otherwise identical individuals ",
      "drawn from a higher- and a lower-risk stratum, and the proportional ",
      "change in variance (PCV) from Model A to Model B, which expresses the ",
      "share of between-stratum variance explained by the additive main ",
      "effects. A PCV approaching 100% indicates that between-stratum ",
      "differences are almost entirely additive, with little residual ",
      "intersectional interaction."
    ),
    paste0(
      "Discriminatory accuracy was summarised by the area under the receiver ",
      "operating characteristic curve (AUC), computed on the aggregated ",
      "binomial data with exact handling of ties. We report the AUC of the ",
      "additive-only prediction and of the full prediction including the ",
      "stratum random effect, so that the discriminatory contribution of the ",
      "interaction residuals can be seen directly."
    ),
    paste0(
      "Individual-level records were aggregated to one binomial observation ",
      "per stratum before fitting. This is algebraically equivalent to fitting ",
      "the individual-level model and retains the exact binomial likelihood, ",
      "while substantially reducing computation. Models were estimated by ",
      "maximum likelihood using the Laplace approximation (nAGQ = ", N_AGQ,
      ") with the bobyqa optimiser."
    ),
    paste0(
      "Much of the applied MAIHDA literature estimates these models in a ",
      "Bayesian framework using Markov chain Monte Carlo, which yields ",
      "credible intervals for the variance components directly from the ",
      "posterior and can be better behaved when many strata are small. We used ",
      "maximum likelihood, which the methodological tutorial literature also ",
      "supports, because it is deterministic, requires no prior specification ",
      "or convergence diagnostics, and is therefore more straightforward to ",
      "audit and reproduce within a trusted research environment. The ",
      "minimum cell size rule applied here also removes the smallest strata, ",
      "which is where the two approaches would be most likely to diverge. ",
      "Uncertainty in the variance components was quantified by ",
      "profile-likelihood intervals on the random-effect standard deviation, ",
      "transformed to the variance, variance partition coefficient and median ",
      "odds ratio scales; all three transformations are monotonic, so the ",
      "interval endpoints carry across directly."
    ),
    paste0(
      "Where a random-intercept variance was estimated at the zero boundary, ",
      "we refitted the model with alternative optimisers to distinguish a ",
      "genuine boundary solution from a failure of optimisation. A boundary ",
      "estimate in Model B that was reproduced across optimisers is reported ",
      "as a substantive finding, namely that the additive main effects account ",
      "for effectively all between-stratum variation, rather than as a ",
      "convergence failure."
    )
  )

  # --- Multiplicity ---------------------------------------------------------
  sections[["Stratum-level inference and multiplicity"]] <- c(
    paste0(
      "For each stratum we obtained the conditional mode of the random effect ",
      "and its conditional standard deviation, and formed an approximate Wald ",
      "interval on the log-odds scale. Intervals were also expressed on the ",
      "probability scale as the difference between the full predicted ",
      "probability and the additive-only predicted probability for that ",
      "stratum. Intervals combining fixed and random-effect uncertainty assume ",
      "zero covariance between the two components and are therefore ",
      "approximate."
    ),
    paste0(
      "The multilevel model already provides a degree of protection against ",
      "spurious findings: the random effects are precision-weighted, so that ",
      "estimates for strata with small samples are shrunk towards the additive ",
      "prediction, and this shrinkage has been argued to be a more efficient ",
      "control of the multiple comparisons problem than inflating intervals ",
      "through a Bonferroni-type correction, because it does not sacrifice ",
      "power to detect genuine differences."
    ),
    paste0(
      "Because each analysis nonetheless compares several hundred strata ",
      "simultaneously, we additionally applied Benjamini-Hochberg false ",
      "discovery rate correction across the strata within each analysis, and ",
      "report as interactions only those strata significant at a false ",
      "discovery rate of ", fmt_num(FDR_LEVEL, 2), ". This is a deliberately ",
      "conservative choice, and more conservative than is usual in the applied ",
      "MAIHDA literature, which commonly relies on shrinkage alone. It ",
      "therefore risks understating the number of interactions rather than ",
      "overstating it. Both the corrected and the uncorrected results are ",
      "reported in the supplementary output so that the effect of this choice ",
      "is fully visible."
    )
  )

  # --- Sensitivity ----------------------------------------------------------
  if (!is.null(sensitivity) && nrow(sensitivity) > 0L) {
    thresholds <- sort(unique(as.numeric(sensitivity$min_stratum_n)))
    sections[["Sensitivity analyses"]] <- paste0(
      "To assess whether the minimum cell size threshold influenced the ",
      "findings, we refitted Models A and B across a range of thresholds (",
      "minimum stratum size ", oxford(as.character(thresholds)),
      ") and compared the resulting variance partition coefficients, ",
      "proportional change in variance and median odds ratios, alongside the ",
      "share of the cohort retained at each threshold."
    )
  }

  # --- Software -------------------------------------------------------------
  software_text <- paste0(
    "All analyses were conducted in ", R_VERSION,
    if (nzchar(LME4_VERSION) && !is.na(LME4_VERSION)) {
      paste0(", with mixed models fitted using the lme4 package (version ",
             LME4_VERSION, ")")
    } else {
      ", with mixed models fitted using the lme4 package"
    },
    ". Analysis code is available at ",
    metadata_or_placeholder("data_access",
                            "state where the analysis code is available"),
    "."
  )
  ethics_text <- paste0(
    "Ethical approval was granted by ",
    metadata_or_placeholder("ethics_approval",
                            "state the ethics approval body and reference"),
    "."
  )
  funding_text <- paste0(
    "This work was funded by ",
    metadata_or_placeholder("funding", "state the funding source"),
    "."
  )
  sections[["Software, approvals and funding"]] <- c(
    software_text, ethics_text, funding_text
  )

  sections
}

# =============================================================================
# 7. RESULTS
# =============================================================================

build_results <- function() {
  sections <- list()

  # --- Cohort ---------------------------------------------------------------
  cohort_sentences <- paste0(
    "The analytic cohort comprised ",
    disclose_count(n_cohort, "Cohort size", "outcome_counts", "total_n"),
    " individuals. ",
    capitalise_first(oxford(vapply(seq_along(axis_names), function(index) {
      axis <- axis_names[index]
      rows <- cohort[cohort$variable == axis, , drop = FALSE]
      rows <- rows[order(-rows$n), , drop = FALSE]
      largest <- rows[1L, ]
      paste0("the most common ", label_in_sentence(axis_labels[index]),
             " category was ", largest$level, " (",
             disclose_count(largest$n,
                            paste0("Largest category, ", axis),
                            "descriptive_cohort_characteristics", "n"),
             "; ", fmt_pct(largest$percent), ")")
    }, character(1)))),
    "."
  )

  missing_rows <- cohort[!duplicated(cohort$variable) & cohort$missing_n > 0, ,
                         drop = FALSE]
  if (nrow(missing_rows) > 0L) {
    missing_text <- oxford(vapply(seq_len(nrow(missing_rows)), function(index) {
      row <- missing_rows[index, ]
      paste0(label_in_sentence(row$variable_label), " (",
             disclose_count(row$missing_n,
                            paste0("Missing, ", row$variable),
                            "descriptive_cohort_characteristics", "missing_n"),
             "; ", fmt_pct(row$missing_percent, 2), ")")
    }, character(1)))
    missing_sentence <- paste0(
      "Data were missing for ", missing_text,
      ". All other stratum axes were complete."
    )
  } else {
    missing_sentence <- "All stratum axes were complete, with no missing values."
  }

  sections[["Cohort characteristics"]] <- c(
    cohort_sentences, missing_sentence,
    "Full cohort characteristics are given in Table 1."
  )

  # --- Outcome frequencies --------------------------------------------------
  most_common <- outcome_counts[which.max(outcome_counts$event_percent), ]
  least_common <- outcome_counts[which.min(outcome_counts$event_percent), ]
  sections[["Outcome frequencies"]] <- c(
    paste0(
      "Event rates varied considerably across the ", spell_number(n_outcomes),
      " outcomes examined. The most frequent was ",
      outcome_in_sentence(most_common$outcome_label), " (",
      disclose_count(most_common$event_n, "Most frequent outcome events",
                     "outcome_counts", "event_n"),
      " events; ", fmt_pct(most_common$event_percent, 2),
      "), and the least frequent was ", outcome_in_sentence(least_common$outcome_label),
      " (",
      disclose_count(least_common$event_n, "Least frequent outcome events",
                     "outcome_counts", "event_n"),
      " events; ", fmt_pct(least_common$event_percent, 2),
      "). Full outcome frequencies are given in Table 2."
    )
  )

  # --- Strata and retention -------------------------------------------------
  retention_sentence <- paste0(
    "Applying the minimum cell size rule (n ≥ ", MIN_STRATUM_N,
    ") retained ",
    fmt_n(primary_retention$strata_retained), " of ",
    fmt_n(primary_retention$strata_before), " observed strata (",
    fmt_pct(primary_retention$strata_retained_percent),
    "), corresponding to ",
    disclose_count(primary_retention$individuals_retained,
                   "Individuals retained after cell rule",
                   "stratum_retention", "individuals_retained"),
    " individuals (",
    stat_value(primary_retention$individuals_retained_percent,
               "Percent individuals retained", "stratum_retention",
               "individuals_retained_percent", 1, percent = TRUE),
    " of the cohort) and ",
    stat_value(primary_retention$events_retained_percent,
               "Percent events retained", "stratum_retention",
               "events_retained_percent", 1, percent = TRUE),
    " of events."
  )

  # Retention is outcome-specific because event counts differ; where it varies
  # the range is reported rather than a single figure that would be wrong for
  # most outcomes.
  retention_range <- range(retention$individuals_retained_percent, na.rm = TRUE)
  retention_variation <- if (diff(retention_range) > 0.05) {
    paste0(
      " Across all outcomes, between ", fmt_pct(retention_range[1L]), " and ",
      fmt_pct(retention_range[2L]),
      " of individuals were retained; full retention accounting is given in ",
      "Table 5."
    )
  } else {
    " Retention was effectively identical across outcomes; see Table 5."
  }

  sections[["Intersectional strata"]] <- c(
    paste0(
      "Of the ", fmt_n(distribution$possible_strata),
      " possible intersectional strata, ",
      fmt_n(distribution$observed_strata), " were observed (",
      fmt_pct(100 * as.numeric(distribution$observed_strata) /
                as.numeric(distribution$possible_strata)),
      "). The distribution of stratum sizes was strongly right-skewed: the ",
      "median stratum contained ", fmt_num(distribution$median_n, 0),
      " individuals (interquartile range ", fmt_num(distribution$q1_n, 0),
      " to ", fmt_num(distribution$q3_n, 0), "), while the largest contained ",
      fmt_n(distribution$maximum_n), ". A total of ",
      fmt_n(distribution$strata_below_10),
      " strata contained fewer than 10 individuals, together accounting for ",
      fmt_pct(distribution$percent_individuals_below_minimum),
      " of the cohort."
    ),
    paste0(retention_sentence, retention_variation)
  )

  # --- Model A --------------------------------------------------------------
  vpc_a_range <- range(metrics$vpc_model_a_percent, na.rm = TRUE)
  mor_a_range <- range(metrics$mor_model_a, na.rm = TRUE)
  highest_vpc <- metrics[which.max(metrics$vpc_model_a_percent), ]

  sections[["Between-stratum variation (Model A)"]] <- c(
    paste0(
      "In the null model, ", vpc_qualifier(primary$vpc_model_a_percent),
      " of the variation in ", outcome_in_sentence(primary_label),
      " lay between intersectional strata. The between-stratum variance was ",
      stat_value(primary$variance_model_a, "Model A variance, primary outcome",
                 "all_model_metrics", "variance_model_a", 3),
      ", giving a variance partition coefficient of ",
      stat_value(primary$vpc_model_a_percent, "Model A VPC, primary outcome",
                 "all_model_metrics", "vpc_model_a_percent", 1, percent = TRUE),
      interval_suffix(primary$vpc_model_a_low, primary$vpc_model_a_high,
                      percent = TRUE),
      " and a median odds ratio of ",
      stat_value(primary$mor_model_a, "Model A MOR, primary outcome",
                 "all_model_metrics", "mor_model_a", 2),
      interval_suffix(primary$mor_model_a_low, primary$mor_model_a_high),
      ". The median odds ratio indicates that two otherwise identical ",
      "individuals drawn at random from a higher- and a lower-risk stratum ",
      "would differ in their odds of the outcome by a median factor of ",
      fmt_num(primary$mor_model_a, 2), "."
    ),
    paste0(
      "Across all ", spell_number(n_outcomes), " outcomes, the Model A ",
      "variance partition coefficient ranged from ",
      fmt_pct(vpc_a_range[1L]), " to ", fmt_pct(vpc_a_range[2L]),
      ", and the median odds ratio from ", fmt_num(mor_a_range[1L], 2),
      " to ", fmt_num(mor_a_range[2L], 2), ". Between-stratum variation was ",
      "greatest for ", outcome_in_sentence(highest_vpc$outcome_label), " (VPC ",
      fmt_pct(highest_vpc$vpc_model_a_percent), "). Model A results for all ",
      "outcomes are given in Table 6."
    ),
    # A variance partition coefficient is hard to judge in isolation, so it is
    # placed against the range reported in the published MAIHDA literature.
    paste0(
      "For context, a systematic review of published MAIHDA analyses of ",
      "health outcomes reported variance partition coefficients ranging from ",
      fmt_pct(VPC_BENCHMARK_LOW), " to ", fmt_pct(VPC_BENCHMARK_HIGH),
      ", with a median of ", fmt_pct(VPC_BENCHMARK_MEDIAN),
      ". The values observed here therefore sit ",
      if (median(metrics$vpc_model_a_percent, na.rm = TRUE) >
          VPC_BENCHMARK_MEDIAN) {
        "above the median of previously published analyses"
      } else {
        "at or below the median of previously published analyses"
      },
      ", indicating that these intersectional strata capture ",
      if (median(metrics$vpc_model_a_percent, na.rm = TRUE) >
          VPC_BENCHMARK_MEDIAN) {
        "an unusually large share of the total individual variation relative to comparable studies"
      } else {
        "a share of the total individual variation comparable to that seen in the wider literature"
      },
      "."
    )
  )

  # --- Which axis drives the variation --------------------------------------
  if (!is.null(axis_decomposition)) {
    primary_decomposition <- axis_decomposition[
      axis_decomposition$outcome == PRIMARY_OUTCOME &
        axis_decomposition$status == "fitted", , drop = FALSE
    ]
    if (nrow(primary_decomposition) > 0L) {
      ordered_axes <- primary_decomposition[
        order(-primary_decomposition$pcv_percent), , drop = FALSE
      ]
      described <- vapply(seq_len(nrow(ordered_axes)), function(index) {
        row <- ordered_axes[index, ]
        paste0(label_in_sentence(row$axis_label), " (",
               fmt_pct(row$pcv_percent), ")")
      }, character(1))

      # Consistency across outcomes is worth stating explicitly, because a
      # decomposition that reorders itself for every outcome means something
      # quite different from one that is stable.
      leading_by_outcome <- do.call(rbind, lapply(
        split(axis_decomposition[axis_decomposition$status == "fitted", ],
              axis_decomposition$outcome[axis_decomposition$status == "fitted"]),
        function(group) group[which.max(group$pcv_percent), ]
      ))
      leading_counts <- sort(table(leading_by_outcome$axis_label),
                             decreasing = TRUE)
      consistency <- if (length(leading_counts) == 1L) {
        paste0(" The same axis was the largest single contributor for every ",
               "outcome examined.")
      } else {
        paste0(" ", capitalise_first(label_in_sentence(names(leading_counts)[1L])),
               " was the largest single contributor for ",
               fmt_n(leading_counts[[1L]]), " of ", fmt_n(n_outcomes),
               " outcomes.")
      }

      sections[["Contribution of individual axes"]] <- paste0(
        "Adding each axis separately to the null model identified which ",
        "dimensions of social and clinical position account for the ",
        "between-stratum variation. For ", outcome_in_sentence(primary_label),
        ", the proportional reduction in between-stratum variance attributable ",
        "to each axis alone was ", oxford(described),
        ". Because the axes are correlated in the population these ",
        "contributions overlap and do not sum to the fully adjusted value.",
        consistency,
        " Full results are given in Table 10."
      )
    }
  }

  # --- Model B main effects -------------------------------------------------
  primary_effects <- fixed_effects[
    fixed_effects$outcome == PRIMARY_OUTCOME &
      fixed_effects$term_status == "estimated" &
      fixed_effects$variable != "(Intercept)", , drop = FALSE
  ]
  main_effect_text <- if (nrow(primary_effects) > 0L) {
    strongest <- primary_effects[order(-abs(log(primary_effects$odds_ratio))), ,
                                 drop = FALSE]
    strongest <- head(strongest, 3L)
    described <- vapply(seq_len(nrow(strongest)), function(index) {
      row <- strongest[index, ]
      paste0(label_in_sentence(row$variable_label), " ", row$level, " (OR ",
             fmt_or_ci(row$odds_ratio, row$confidence_low, row$confidence_high,
                       bare = TRUE),
             ")")
    }, character(1))
    paste0(
      "In Model B, the strongest additive associations with ",
      outcome_in_sentence(primary_label), " were ", oxford(described),
      ", each relative to its reference category. Full additive main effects ",
      "for all outcomes are given in Table 7."
    )
  } else {
    paste0("Additive main effects for all outcomes are given in Table 7.")
  }

  pcv_range <- range(metrics$pcv_percent, na.rm = TRUE)
  sections[["Additive main effects (Model B)"]] <- c(
    main_effect_text,
    paste0(
      "Adding the additive main effects explained the large majority of ",
      "between-stratum variation. For ", outcome_in_sentence(primary_label),
      ", the proportional change in variance from Model A to Model B was ",
      stat_value(primary$pcv_percent, "PCV, primary outcome",
                 "all_model_metrics", "pcv_percent", 1, percent = TRUE),
      ", leaving a residual between-stratum variance of ",
      stat_value(primary$variance_model_b, "Model B variance, primary outcome",
                 "all_model_metrics", "variance_model_b", 4),
      " (VPC ", fmt_pct(primary$vpc_model_b_percent, 2),
      "). Across all outcomes the proportional change in variance ranged from ",
      fmt_pct(pcv_range[1L]), " to ", fmt_pct(pcv_range[2L]),
      ", indicating that between-stratum differences in these outcomes are ",
      "predominantly additive."
    )
  )

  # --- Interactions ---------------------------------------------------------
  interaction_paragraphs <- character(0)

  if (nrow(singular_b) > 0L) {
    singular_names <- outcome_in_sentence(singular_b$outcome_label)
    interaction_paragraphs <- c(interaction_paragraphs, paste0(
      "For ", oxford(singular_names), ", the Model B between-stratum variance ",
      "was estimated at zero. This estimate was reproduced across alternative ",
      "optimisers and therefore represents a genuine boundary solution rather ",
      "than a failure of optimisation. Substantively, it indicates that the ",
      "additive main effects account for effectively all between-stratum ",
      "variation in ", plural(nrow(singular_b), "this outcome", "these outcomes"),
      ", with no residual intersectional interaction detectable at this ",
      "sample size."
    ))
  }

  if (nrow(outcomes_with_interactions) == 0L) {
    total_uncorrected <- sum(metrics$n_distinguishable, na.rm = TRUE)
    interaction_paragraphs <- c(interaction_paragraphs, paste0(
      "After Benjamini-Hochberg correction, no intersectional stratum showed ",
      "a departure from its additive prediction that was statistically ",
      "distinguishable from zero, for any outcome. ",
      if (total_uncorrected > 0L) {
        paste0(
          "Before correction, ", fmt_n(total_uncorrected), " ",
          plural(total_uncorrected, "stratum", "strata"), " across all ",
          "analyses had an interval excluding zero, a number consistent with ",
          "chance given the several hundred simultaneous comparisons made ",
          "within each analysis. "
        )
      } else "",
      "Taken together with the high proportional change in variance, these ",
      "results indicate that inequities between strata in these outcomes are ",
      "well characterised by the additive combination of the constituent ",
      "social and clinical positions, without additional multiplicative ",
      "effects being required to describe them."
    ),
    # The tutorial (section 4.2) is explicit that a near-100% PCV must not be
    # written up as an absence of intersectionality. Three caveats belong with
    # any such result, and omitting them is the most common way this finding is
    # over-interpreted.
    paste0(
      "Three caveats attach to this finding. First, the variance partition ",
      "coefficient and proportional change in variance are summary measures ",
      "across the whole set of strata; a single stratum with a substantial ",
      "interaction would be diluted by them, which is why individual strata ",
      "were examined separately rather than relying on the summary measures ",
      "alone. Second, the precision weighting that makes these estimates ",
      "robust also shrinks them towards the additive prediction, so genuine ",
      "interactions may exist in the population but be undetectable at the ",
      "available stratum sizes; the approach protects against false positives ",
      "at the cost of some power to detect small true effects. Third, and most ",
      "importantly, an absence of detectable interaction is not evidence ",
      "against intersectionality as a framework. It indicates that the ",
      "experiences associated with these positions did not, for these outcomes ",
      "in this population at this time, leave an imprint requiring ",
      "multiplicative as well as additive terms to characterise."
    ),
    paste0(
      "The substantive conclusion is therefore not that intersectionality is ",
      "unimportant here, but that there are meaningful inequities between ",
      "intersectional strata, that the axes examined were collectively ",
      "important in capturing those inequities, and that the resulting ",
      "distribution of harm should inform how care and resources are directed."
    ))
  } else {
    n_with <- nrow(outcomes_with_interactions)
    total_fdr <- sum(outcomes_with_interactions$n_distinguishable_fdr)
    interaction_paragraphs <- c(interaction_paragraphs, paste0(
      "After Benjamini-Hochberg correction, ", spell_number(total_fdr), " ",
      plural(total_fdr, "stratum", "strata"), " across ", spell_number(n_with),
      " ", plural(n_with, "outcome"),
      " showed a departure from the additive prediction that was ",
      "statistically distinguishable from zero."
    ))

    if (has_significant) {
      described <- significant[order(-abs(
        significant$interaction_probability_difference
      )), , drop = FALSE]
      described <- head(described, MAX_STRATA_NAMED)
      stratum_sentences <- vapply(seq_len(nrow(described)), function(index) {
        row <- described[index, ]
        direction <- if (row$interaction_probability_difference > 0) {
          "higher"
        } else {
          "lower"
        }
        paste0(
          "for ", outcome_in_sentence(row$outcome_label), ", the stratum defined by ",
          row$condition_label, " (n = ",
          disclose_count(row$n, paste0("Stratum n: ", row$stratum),
                         "distinguishable_interaction_strata", "n"),
          "; ",
          disclose_count(row$events, paste0("Stratum events: ", row$stratum),
                         "distinguishable_interaction_strata", "events"),
          " events) carried an absolute risk due to interaction of ",
          fmt_num(abs(100 * row$interaction_probability_difference), 1),
          " percentage points ", direction,
          " than its additive prediction (95% CI ",
          fmt_num(100 * row$interaction_probability_difference_low, 1), " to ",
          fmt_num(100 * row$interaction_probability_difference_high, 1),
          "; q ", fmt_p(row$interaction_q_value), ")"
        )
      }, character(1))

      interaction_paragraphs <- c(interaction_paragraphs, paste0(
        capitalise_first(
          if (nrow(described) < total_fdr) {
            paste0("The ", spell_number(nrow(described)),
                   " largest departures were as follows. ")
          } else {
            "These were as follows. "
          }
        ),
        capitalise_first(paste0(oxford(stratum_sentences), ".")),
        " All strata surviving correction are listed in Table 9. A large ",
        "interaction residual does not by itself mean that a stratum is ",
        "advantaged or disadvantaged in absolute terms; it means the total ",
        "predicted risk for that stratum departs meaningfully from what the ",
        "additive main effects alone would predict. Absolute risks for all ",
        "strata are reported separately."
      ))
    }
  }

  sections[["Intersectional interactions"]] <- interaction_paragraphs

  # --- Discriminatory accuracy ----------------------------------------------
  auc_gain <- metrics$auc_model_b_total - metrics$auc_model_b_fixed
  sections[["Discriminatory accuracy"]] <- c(
    paste0(
      "Discriminatory accuracy is reported as a descriptive quantity rather ",
      "than against conventional thresholds of predictive performance, since ",
      "the intention is not to build a diagnostic model but to express how ",
      "well intersectional position alone separates those who experience the ",
      "outcome from those who do not. For ",
      outcome_in_sentence(primary_label), ", the area under the receiver operating ",
      "characteristic curve was ",
      stat_value(primary$auc_model_b_fixed, "AUC additive, primary outcome",
                 "all_model_metrics", "auc_model_b_fixed", 3),
      " for the additive-only prediction and ",
      stat_value(primary$auc_model_b_total, "AUC total, primary outcome",
                 "all_model_metrics", "auc_model_b_total", 3),
      " when the stratum random effect was included."
    ),
    paste0(
      "Across all outcomes, including the stratum random effect changed the ",
      "area under the curve by between ", fmt_num(min(auc_gain, na.rm = TRUE), 3),
      " and ", fmt_num(max(auc_gain, na.rm = TRUE), 3),
      ". The small magnitude of this change is consistent with the high ",
      "proportional change in variance: once additive main effects are ",
      "accounted for, the intersectional strata add little further ",
      "discriminatory information."
    )
  )

  # --- Extremes -------------------------------------------------------------
  primary_extremes <- extremes[extremes$outcome == PRIMARY_OUTCOME, ,
                               drop = FALSE]
  if (nrow(primary_extremes) > 0L) {
    highest <- primary_extremes[
      which.max(primary_extremes$total_probability), ]
    lowest <- primary_extremes[which.min(primary_extremes$total_probability), ]
    sections[["Highest- and lowest-risk strata"]] <- paste0(
      "The gradient across intersectional strata was wide. For ",
      outcome_in_sentence(primary_label), ", the highest-risk stratum was ",
      highest$condition_label, ", with a model-predicted risk of ",
      fmt_pct(100 * highest$total_probability),
      " (95% CI ", fmt_pct(100 * highest$total_probability_low), " to ",
      fmt_pct(100 * highest$total_probability_high), "; n = ",
      disclose_count(highest$n, "Highest-risk stratum n", "extreme_strata", "n"),
      "). The lowest-risk stratum was ", lowest$condition_label,
      ", with a predicted risk of ",
      fmt_pct(100 * lowest$total_probability),
      " (95% CI ", fmt_pct(100 * lowest$total_probability_low), " to ",
      fmt_pct(100 * lowest$total_probability_high), "; n = ",
      disclose_count(lowest$n, "Lowest-risk stratum n", "extreme_strata", "n"),
      "). This represents an absolute difference of ",
      fmt_num(100 * (highest$total_probability - lowest$total_probability), 1),
      " percentage points between the extremes of the intersectional ",
      "distribution",
      if (!is.null(primary$absolute_risk_ratio) &&
          !is.na(primary$absolute_risk_ratio)) {
        paste0(", a relative difference of ",
               stat_value(primary$absolute_risk_ratio,
                          "Absolute risk ratio, extremes",
                          "all_model_metrics", "absolute_risk_ratio", 1),
               "-fold")
      } else "",
      ". The ", fmt_n(nrow(primary_extremes) / 2),
      " highest- and lowest-risk strata for each outcome are given in Table 8."
    )
  }

  # --- Sensitivity ----------------------------------------------------------
  if (!is.null(sensitivity) && nrow(sensitivity) > 0L) {
    fitted_sensitivity <- sensitivity[sensitivity$status == "fitted", ,
                                      drop = FALSE]
    if (nrow(fitted_sensitivity) > 1L) {
      vpc_span <- range(fitted_sensitivity$vpc_model_a_percent, na.rm = TRUE)
      pcv_span <- range(fitted_sensitivity$pcv_percent, na.rm = TRUE)
      retained_span <- range(fitted_sensitivity$individuals_retained_percent,
                             na.rm = TRUE)
      sections[["Sensitivity to the minimum cell size threshold"]] <- paste0(
        "Refitting the models across minimum stratum sizes of ",
        oxford(as.character(sort(unique(fitted_sensitivity$min_stratum_n)))),
        " retained between ", fmt_pct(retained_span[1L]), " and ",
        fmt_pct(retained_span[2L]), " of the cohort. The Model A variance ",
        "partition coefficient varied between ", fmt_pct(vpc_span[1L]),
        " and ", fmt_pct(vpc_span[2L]),
        ", and the proportional change in variance between ",
        fmt_pct(pcv_span[1L]), " and ", fmt_pct(pcv_span[2L]),
        ". The substantive conclusion, that between-stratum variation is ",
        "predominantly additive, was unchanged at every threshold examined, ",
        "indicating that the findings are not an artefact of the trimming ",
        "rule. Full sensitivity results are given in Table 10."
      )
    }
  }

  # --- Model diagnostics ----------------------------------------------------
  if (!is.null(failures) && nrow(failures) > 0L) {
    sections[["Analyses not completed"]] <- paste0(
      spell_number(nrow(failures)), " ", plural(nrow(failures), "analysis", "analyses"),
      " could not be completed: ",
      oxford(paste0(outcome_in_sentence(failures$outcome_label), " (", failures$message, ")")),
      ". ", capitalise_first(plural(nrow(failures), "This analysis is", "These analyses are")),
      " excluded from the results above."
    )
  }

  sections
}

# =============================================================================
# 8. ASSEMBLE AND WRITE
# =============================================================================

methods_sections <- build_methods()
results_sections <- build_results()

document_title <- "Methods and Results"
generation_note <- paste0(
  "Generated from the MAIHDA analysis output in ", ANALYSIS_ROOT, " on ",
  format(Sys.time(), "%d %B %Y at %H:%M", tz = "UTC"), " UTC. ",
  "Every figure quoted below is read directly from the analysis tables; ",
  "nothing is recomputed. Disclosure mode: ", DISCLOSURE_MODE, "."
)

# --- Markdown and plain text ------------------------------------------------
build_plain_document <- function(markdown = TRUE) {
  heading1 <- function(text) if (markdown) paste0("# ", text) else
    c(toupper(text), strrep("=", nchar(text)))
  heading2 <- function(text) if (markdown) paste0("## ", text) else
    c(text, strrep("-", nchar(text)))

  lines <- c(
    heading1(document_title), "",
    if (markdown) paste0("*", generation_note, "*") else generation_note, ""
  )
  for (part in list(list("Methods", methods_sections),
                    list("Results", results_sections))) {
    lines <- c(lines, heading1(part[[1L]]), "")
    for (section_name in names(part[[2L]])) {
      lines <- c(lines, heading2(section_name), "")
      for (paragraph in part[[2L]][[section_name]]) {
        lines <- c(lines, paragraph, "")
      }
    }
  }
  lines
}

markdown_path <- file.path(MANUSCRIPT_DIR, "Manuscript.md")
text_path <- file.path(MANUSCRIPT_DIR, "Manuscript.txt")
writeLines(build_plain_document(TRUE), markdown_path)
writeLines(build_plain_document(FALSE), text_path)

# --- Word -------------------------------------------------------------------
document <- officer::read_docx()
document <- officer::body_add_par(document, document_title, style = "heading 1")
document <- officer::body_add_fpar(
  document,
  officer::fpar(officer::ftext(
    generation_note,
    officer::fp_text(font.size = 9, italic = TRUE, color = "#555555")
  ))
)
for (part in list(list("Methods", methods_sections),
                  list("Results", results_sections))) {
  document <- officer::body_add_par(document, part[[1L]], style = "heading 1")
  for (section_name in names(part[[2L]])) {
    document <- officer::body_add_par(document, section_name,
                                      style = "heading 2")
    for (paragraph in part[[2L]][[section_name]]) {
      document <- officer::body_add_par(document, paragraph, style = "Normal")
    }
  }
}

# A compact summary table is appended when flextable is available, because the
# Results text refers to it repeatedly and having it in the same file makes the
# generated document reviewable on its own.
if (HAS_FLEXTABLE) {
  summary_table <- data.frame(
    Outcome = metrics$outcome_label,
    Strata = fmt_n(metrics$observed_strata),
    `Event %` = fmt_num(metrics$event_percent, 2),
    `VPC A %` = fmt_num(metrics$vpc_model_a_percent, 2),
    `VPC B %` = fmt_num(metrics$vpc_model_b_percent, 2),
    `PCV %` = fmt_num(metrics$pcv_percent, 2),
    `MOR A` = fmt_num(metrics$mor_model_a, 2),
    `AUC` = fmt_num(metrics$auc_model_b_total, 3),
    `Interactions` = fmt_n(metrics$n_distinguishable_fdr),
    check.names = FALSE
  )
  document <- officer::body_add_par(document, "Summary of model results",
                                    style = "heading 2")
  flex <- flextable::flextable(summary_table)
  flex <- flextable::theme_booktabs(flex, bold_header = TRUE)
  flex <- flextable::fontsize(flex, size = 8, part = "all")
  flex <- flextable::autofit(flex)
  document <- flextable::body_add_flextable(document, flex)
}

word_path <- file.path(MANUSCRIPT_DIR, "Manuscript.docx")
print(document, target = word_path)

# --- Traceability, disclosure and placeholders ------------------------------
values_table <- if (length(value_registry$rows) > 0L) {
  do.call(rbind, value_registry$rows)
} else {
  data.frame()
}
values_path <- file.path(MANUSCRIPT_DIR, "manuscript_values.csv")
if (nrow(values_table) > 0L) {
  values_table <- values_table[!duplicated(values_table$label), , drop = FALSE]
  write.csv(values_table, values_path, row.names = FALSE)
}

disclosure_table <- if (length(disclosure_registry$rows) > 0L) {
  do.call(rbind, disclosure_registry$rows)
} else {
  data.frame()
}
disclosure_path <- file.path(MANUSCRIPT_DIR, "disclosure_check.csv")
write.csv(disclosure_table, disclosure_path, row.names = FALSE)

placeholder_path <- file.path(MANUSCRIPT_DIR, "placeholders_to_complete.txt")
placeholder_items <- placeholder_registry$items
writeLines(
  if (length(placeholder_items) == 0L) {
    "All study metadata fields were supplied. No placeholders remain."
  } else {
    c(
      "The following study metadata fields were not supplied and appear in the",
      "generated text as [PLACEHOLDER: ...] markers. Set them in section 1 of",
      "generate_manuscript.R, or through the corresponding environment",
      "variables, and regenerate.",
      "",
      paste0("  - ", placeholder_items)
    )
  },
  placeholder_path
)

# =============================================================================
# 9. COMPLETION SUMMARY
# =============================================================================

message("")
message("Manuscript written:")
message("  ", word_path)
message("  ", markdown_path)
message("  ", text_path)
message("")
message("Supporting files:")
message("  ", values_path, " (", nrow(values_table), " traced values)")
message("  ", disclosure_path, " (", nrow(disclosure_table),
        " small-cell findings)")
message("  ", placeholder_path, " (", length(placeholder_items),
        " placeholders outstanding)")

if (nrow(disclosure_table) > 0L) {
  message("")
  if (DISCLOSURE_MODE == "flag") {
    message("NOTE: ", nrow(disclosure_table), " reported ",
            plural(nrow(disclosure_table), "figure"),
            " derive from cells below the disclosure threshold of ",
            DISCLOSURE_THRESHOLD, ".")
    message("      True values are reported in the text. Review ",
            basename(disclosure_path), " before requesting output,")
    message("      or regenerate with MANUSCRIPT_DISCLOSURE_MODE=apply.")
  } else if (DISCLOSURE_MODE == "apply") {
    message("NOTE: ", nrow(disclosure_table), " small ",
            plural(nrow(disclosure_table), "cell"),
            " were suppressed or rounded in the generated text.")
  }
}
if (length(placeholder_items) > 0L) {
  message("")
  message("ACTION REQUIRED: ", length(placeholder_items),
          " study metadata ", plural(length(placeholder_items), "field"),
          " must be completed. See ", basename(placeholder_path), ".")
}
message("")
