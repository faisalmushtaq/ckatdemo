# =============================================================================
# Check (and optionally repair) the Stan toolchain before a Bayesian MAIHDA run
# =============================================================================
#
# PURPOSE
# -------
# Bayesian estimation needs brms, rstan and a working C++ toolchain with Boost
# headers. When any of that is missing the failure arrives from deep inside a
# Stan compile and is close to unreadable. This script checks each requirement
# in turn, reports exactly which one is unmet, and prints the specific command
# that fixes it on this machine.
#
# It is worth running once on the analysis machine before committing to an
# overnight run.
#
# USAGE
# -----
#   Rscript tools/check_stan_toolchain.R
#
# To attempt the safe repairs automatically:
#
#   Rscript tools/check_stan_toolchain.R --fix
#
# THE BOOST PROBLEM ON DEBIAN AND UBUNTU
# --------------------------------------
# The usual cause of "Boost not found; call install.packages('BH')" on a
# Debian or Ubuntu machine is that the distribution's r-cran-bh package is a
# shim. It installs an R package named BH containing metadata and help files
# but NO headers at all, and declares a dependency on the system libboost-dev
# instead. rstan looks for the headers inside the R package, does not find
# them, and stops.
#
# Two repairs work, and this script can apply the second without network
# access, which matters in a trusted research environment:
#
#   1. Install the real BH package from CRAN, which ships the Boost headers
#      inside the R package where rstan expects them.
#   2. Keep the distribution packages and bridge the gap, by creating the
#      include directory the shim omits and pointing it at the system Boost
#      headers that libboost-dev already provides.
#
# Repair 2 is not a bodge: libboost-dev genuinely supplies the same headers,
# just at a path rstan does not consult.
# =============================================================================

arguments <- commandArgs(trailingOnly = TRUE)
APPLY_FIXES <- "--fix" %in% arguments

results <- list()
record <- function(check, ok, detail, remedy = "") {
  results[[length(results) + 1L]] <<- list(
    check = check, ok = ok, detail = detail, remedy = remedy
  )
  symbol <- if (isTRUE(ok)) "PASS" else if (is.na(ok)) "WARN" else "FAIL"
  cat(sprintf("[%s] %-28s %s\n", symbol, check, detail))
  invisible(ok)
}

cat("Stan toolchain check for Bayesian MAIHDA\n")
cat(strrep("=", 68), "\n\n", sep = "")

# --- R version --------------------------------------------------------------
record("R version", TRUE, as.character(getRversion()))

# --- C++ compiler -----------------------------------------------------------
compiler <- Sys.which("g++")
if (!nzchar(compiler)) compiler <- Sys.which("clang++")
record(
  "C++ compiler", nzchar(compiler),
  if (nzchar(compiler)) compiler else "no g++ or clang++ on PATH",
  "Install build tools, e.g. apt-get install build-essential"
)

# --- make -------------------------------------------------------------------
record("make", nzchar(Sys.which("make")),
       if (nzchar(Sys.which("make"))) Sys.which("make") else "not found",
       "apt-get install make")

# --- R packages -------------------------------------------------------------
for (package in c("rstan", "brms", "RcppEigen", "StanHeaders", "Rcpp")) {
  present <- requireNamespace(package, quietly = TRUE)
  record(
    paste("package:", package), present,
    if (present) as.character(utils::packageVersion(package)) else "not installed",
    sprintf("install.packages(\"%s\")", package)
  )
}

# --- Boost headers, the usual culprit ---------------------------------------
boost_include <- system.file("include", package = "BH")
boost_headers_present <- nzchar(boost_include) &&
  dir.exists(file.path(boost_include, "boost"))

system_boost <- c("/usr/include/boost", "/usr/local/include/boost")
system_boost <- system_boost[dir.exists(system_boost)]

if (boost_headers_present) {
  record("Boost headers (BH)", TRUE,
         paste("found at", file.path(boost_include, "boost")))
} else {
  bh_installed <- requireNamespace("BH", quietly = TRUE) ||
    dir.exists(file.path(.libPaths(), "BH"))[1L]
  detail <- if (bh_installed) {
    paste0("BH is installed but ships no headers -- this is the ",
           "distribution shim, not the CRAN package")
  } else {
    "BH is not installed"
  }
  record("Boost headers (BH)", FALSE, detail)

  if (length(system_boost) > 0L) {
    cat("\n       System Boost IS present at ", system_boost[1L], ".\n",
        "       The headers exist; rstan simply does not look there.\n",
        sep = "")
  }
}

# --- Apply repairs ----------------------------------------------------------
if (!boost_headers_present) {
  cat("\n", strrep("-", 68), "\n", sep = "")
  cat("REPAIR OPTIONS\n\n")

  cat("Option 1 -- install the real BH package from CRAN (needs network):\n")
  cat("    install.packages(\"BH\")\n\n")

  # Locate the BH package where it is actually installed rather than assuming
  # the first library path, since the distribution shim commonly sits in a
  # different library from the user's default.
  installed_bh <- tryCatch(find.package("BH"), error = function(condition) NULL)
  bridge_target <- if (nzchar(boost_include)) {
    boost_include
  } else if (!is.null(installed_bh)) {
    file.path(installed_bh, "include")
  } else {
    file.path(.libPaths()[1L], "BH", "include")
  }
  if (length(system_boost) > 0L) {
    cat("Option 2 -- bridge the shim to system Boost (no network needed):\n")
    cat("    mkdir -p ", bridge_target, "\n", sep = "")
    cat("    ln -s ", system_boost[1L], " ", file.path(bridge_target, "boost"),
        "\n\n", sep = "")
  } else {
    cat("Option 2 unavailable: no system Boost found either.\n")
    cat("    apt-get install libboost-dev\n\n")
  }

  if (APPLY_FIXES && length(system_boost) > 0L) {
    cat("Applying option 2 (--fix was given) ...\n")
    created <- tryCatch({
      dir.create(bridge_target, recursive = TRUE, showWarnings = FALSE)
      link_path <- file.path(bridge_target, "boost")
      if (!dir.exists(link_path)) {
        file.symlink(system_boost[1L], link_path)
      }
      dir.exists(file.path(bridge_target, "boost"))
    }, error = function(condition) {
      cat("    Failed: ", conditionMessage(condition), "\n", sep = "")
      FALSE
    })
    if (isTRUE(created)) {
      cat("    Done. Boost headers now visible at ",
          file.path(bridge_target, "boost"), "\n", sep = "")
      boost_headers_present <- TRUE
    } else {
      cat("    Could not create the link. You may need write access to ",
          bridge_target, ", or run the command above with sudo.\n", sep = "")
    }
  } else if (!APPLY_FIXES) {
    cat("Re-run with --fix to apply option 2 automatically.\n")
  }
}

# --- End-to-end compile test ------------------------------------------------
# The only check that really settles it: compile and sample from a trivial
# model. Everything above can pass and this can still fail.
cat("\n", strrep("-", 68), "\n", sep = "")
if (boost_headers_present && requireNamespace("brms", quietly = TRUE)) {
  cat("Compiling a minimal Stan model to confirm the toolchain works.\n")
  cat("This takes a minute or two on a slow processor.\n\n")
  compile_test <- tryCatch({
    test_data <- data.frame(
      events = c(3L, 7L, 5L, 9L), n = c(10L, 20L, 15L, 25L),
      stratum = factor(c("a", "b", "c", "d"))
    )
    suppressMessages(suppressWarnings(
      brms::brm(
        brms::brmsformula(events | trials(n) ~ 1 + (1 | stratum),
                          family = binomial("logit")),
        data = test_data, chains = 1, iter = 200, warmup = 100,
        refresh = 0, seed = 1
      )
    ))
    TRUE
  }, error = function(condition) {
    cat("FAILED: ", conditionMessage(condition), "\n", sep = "")
    FALSE
  })

  if (isTRUE(compile_test)) {
    cat("[PASS] Stan compiled and sampled successfully.\n\n")
    cat("The toolchain is ready. Run the analysis with:\n")
    cat("    Rscript R/run_maihda_analysis.R\n")
  } else {
    cat("\n[FAIL] Stan could not compile.\n\n")
    cat("Fall back to maximum likelihood, which needs no C++ toolchain and\n")
    cat("produces the same outputs:\n")
    cat("    MAIHDA_ENGINE=mle Rscript R/run_maihda_analysis.R\n")
  }
} else {
  cat("Skipping the compile test: prerequisites above are unmet.\n\n")
  cat("If they cannot be resolved, maximum likelihood needs no C++ toolchain\n")
  cat("and produces the same outputs:\n")
  cat("    MAIHDA_ENGINE=mle Rscript R/run_maihda_analysis.R\n")
}

failed <- sum(vapply(results, function(r) isFALSE(r$ok), logical(1)))
cat("\n", strrep("=", 68), "\n", sep = "")
cat(sprintf("%d checks run, %d unmet.\n", length(results), failed))
