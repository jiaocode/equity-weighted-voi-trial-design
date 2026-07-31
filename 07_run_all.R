# ==============================================================================
# 07_run_all.R
# Clean, reproducible driver for the complete main and appendix workflow.
# Every substantive script runs in a separate --vanilla R session. Child stdout
# and stderr are retained and printed when a step fails.
# ==============================================================================

script_directory <- function() {
  arguments <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("^--file=", arguments, value = TRUE)

  if (length(file_argument) == 1L) {
    # Rscript encodes spaces in --file= as "~+~"; decode before resolving.
    dirname(normalizePath(
      gsub("~+~", " ", sub("^--file=", "", file_argument), fixed = TRUE)
    ))
  } else {
    normalizePath(getwd())
  }
}

ROOT <- script_directory()
source(file.path(ROOT, "00_analysis_parameters.R"), local = TRUE)

ANALYSIS_SIGNATURE <- compute_analysis_signature(ROOT)
RUNTIME_SIGNATURE <- compute_runtime_signature()

key_relative_paths <- c(
  file.path(OUTPUT$main_results_directory, "main_wevsi_summary.csv"),
  file.path(OUTPUT$main_results_directory, "main_wenbs_summary.csv"),
  file.path(OUTPUT$main_results_directory, "optimal_allocation_summary.csv"),
  file.path(OUTPUT$appendix_results_directory, "appendix_wevsi_summary.csv"),
  file.path(OUTPUT$appendix_results_directory, "appendix_wenbs_summary.csv"),
  file.path(OUTPUT$appendix_results_directory, "appendix_optimal_allocation_summary.csv")
)

# Preserve key outputs from a directly comparable prior run before cleaning.
comparison_directory <- tempfile("prior_comparable_outputs_")
dir.create(comparison_directory, recursive = TRUE, showWarnings = FALSE)
comparison_available <- FALSE
prior_manifest_path <- file.path(ROOT, OUTPUT$run_manifest)

if (file.exists(prior_manifest_path)) {
  prior_manifest <- tryCatch(
    utils::read.csv(prior_manifest_path, stringsAsFactors = FALSE),
    error = function(e) NULL
  )

  if (!is.null(prior_manifest) && nrow(prior_manifest) == 1L &&
      identical(as.character(prior_manifest$analysis_signature), ANALYSIS_SIGNATURE) &&
      identical(as.character(prior_manifest$runtime_signature), RUNTIME_SIGNATURE)) {
    prior_paths <- file.path(ROOT, key_relative_paths)
    if (all(file.exists(prior_paths))) {
      copied <- file.copy(
        from = prior_paths,
        to = file.path(comparison_directory, basename(prior_paths)),
        overwrite = TRUE
      )
      comparison_available <- all(copied)
    }
  }
}

# Remove every old result and figure so stale files cannot be reused.
unlink(file.path(ROOT, "results"), recursive = TRUE, force = TRUE)
unlink(file.path(ROOT, "figures"), recursive = TRUE, force = TRUE)

dir.create(file.path(ROOT, OUTPUT$main_results_directory), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(ROOT, OUTPUT$appendix_results_directory), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(ROOT, OUTPUT$main_figures_directory), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(ROOT, OUTPUT$appendix_figures_directory), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(ROOT, OUTPUT$logs_directory), recursive = TRUE, showWarnings = FALSE)

rscript <- file.path(R.home("bin"), "Rscript")
if (!file.exists(rscript)) {
  rscript <- Sys.which("Rscript")
}
if (!nzchar(rscript) || !file.exists(rscript)) {
  stop("Rscript could not be located.", call. = FALSE)
}

read_log <- function(path) {
  if (!file.exists(path)) {
    return("[log file was not created]")
  }
  lines <- readLines(path, warn = FALSE)
  if (length(lines) == 0L) {
    return("[log file is empty]")
  }
  paste(lines, collapse = "\n")
}

run_child_script <- function(script_name) {
  script_path <- file.path(ROOT, script_name)
  if (!file.exists(script_path)) {
    stop("Missing workflow script: ", script_name, call. = FALSE)
  }

  stem <- sub("\\.R$", "", script_name)
  stdout_path <- file.path(ROOT, OUTPUT$logs_directory, paste0(stem, "_stdout.log"))
  stderr_path <- file.path(ROOT, OUTPUT$logs_directory, paste0(stem, "_stderr.log"))

  message("Running ", script_name, " ...")
  status <- system2(
    command = rscript,
    args = c("--vanilla", shQuote(script_path)),
    stdout = stdout_path,
    stderr = stderr_path,
    wait = TRUE,
    env = c(
      "OMP_NUM_THREADS=1",
      "OPENBLAS_NUM_THREADS=1",
      "MKL_NUM_THREADS=1",
      "VECLIB_MAXIMUM_THREADS=1",
      "BLIS_NUM_THREADS=1"
    )
  )

  if (!identical(as.integer(status), 0L)) {
    stdout_text <- read_log(stdout_path)
    stderr_text <- read_log(stderr_path)
    stop(
      script_name, " failed with exit status ", status, ".\n\n",
      "----- STDOUT -----\n", stdout_text, "\n\n",
      "----- STDERR -----\n", stderr_text, "\n\n",
      "The complete logs remain in ", file.path(ROOT, OUTPUT$logs_directory), ".",
      call. = FALSE
    )
  }

  message("Completed ", script_name)
  invisible(TRUE)
}

workflow_scripts <- c(
  "03_run_main_analysis.R",
  "05_run_appendix_analysis.R",
  "04_make_main_figures.R",
  "06_make_appendix_figures.R",
  "08_run_validation_tests.R"
)

for (script_name in workflow_scripts) {
  run_child_script(script_name)
}

runtime <- runtime_information()
manifest <- data.frame(
  analysis_version = REPRODUCIBILITY$analysis_version,
  analysis_signature = ANALYSIS_SIGNATURE,
  runtime_signature = RUNTIME_SIGNATURE,
  r_version = runtime$r_version,
  platform = runtime$platform,
  rng_kind = runtime$rng_kind,
  normal_kind = runtime$normal_kind,
  sample_kind = runtime$sample_kind,
  simulation_nsim = SIMULATION$nsim,
  simulation_n_batches = SIMULATION$n_batches,
  appendix_simulation_nsim = APPENDIX_SIMULATION$nsim,
  appendix_simulation_n_batches = APPENDIX_SIMULATION$n_batches,
  burden_A_qaly = TRIAL$burden_A_qaly,
  burden_B_qaly = TRIAL$burden_B_qaly,
  cost_scenarios_usd = paste(TRIAL$recruitment_cost_increment_A_usd, collapse = ";"),
  stringsAsFactors = FALSE
)
utils::write.csv(
  manifest,
  file.path(ROOT, OUTPUT$run_manifest),
  row.names = FALSE,
  na = ""
)
writeLines(
  capture.output(utils::sessionInfo()),
  file.path(ROOT, "results", "sessionInfo.txt")
)

new_paths <- file.path(ROOT, key_relative_paths)
if (!all(file.exists(new_paths))) {
  stop("The workflow finished without creating all key output files.", call. = FALSE)
}
new_hashes <- unname(tools::md5sum(new_paths))

if (comparison_available) {
  old_paths <- file.path(comparison_directory, basename(new_paths))
  old_hashes <- unname(tools::md5sum(old_paths))
  identical_output <- new_hashes == old_hashes
  comparison_status <- ifelse(identical_output, "identical", "different")
} else {
  old_hashes <- rep(NA_character_, length(new_paths))
  identical_output <- rep(NA, length(new_paths))
  comparison_status <- rep("no directly comparable prior run", length(new_paths))
}

reproducibility_check <- data.frame(
  file = key_relative_paths,
  prior_md5 = old_hashes,
  current_md5 = new_hashes,
  byte_identical = identical_output,
  status = comparison_status,
  stringsAsFactors = FALSE
)
utils::write.csv(
  reproducibility_check,
  file.path(ROOT, OUTPUT$reproducibility_check),
  row.names = FALSE,
  na = ""
)

if (comparison_available && any(!identical_output)) {
  stop(
    "A directly comparable prior run was found, but at least one key output ",
    "was not byte-for-byte identical. See ",
    file.path(ROOT, OUTPUT$reproducibility_check), ".",
    call. = FALSE
  )
}

message("Complete workflow finished successfully.")
message("Analysis signature: ", ANALYSIS_SIGNATURE)
message("Results: ", file.path(ROOT, "results"))
message("Figures: ", file.path(ROOT, "figures"))
