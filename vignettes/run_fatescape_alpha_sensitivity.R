#!/usr/bin/env Rscript

# -------------------------
# 0. USER SETTINGS
# -------------------------
BASE_RUNNER <- "FateScape/vignettes/run_fatescape_batch.R"

METHOD_INPUT_ROOT <- "FateScape/data/Sensitivity_analysis"
METHOD_OUTPUT_ROOT <- "FateScape/results/Sensitivity_analysis_alpha"

# Optional: if FateScape is not installed as an R package, set package root here.
# Example: FATESCAPE_LOCAL_PATH <- "D:/projects/lineage tracing/FateScape"

CASE_INDICES <- "all"

# Batch behavior passed into the patched base runner.
SKIP_FINISHED <- TRUE
FORCE_RERUN <- FALSE
STOP_ON_ERROR <- FALSE

# Keep the same core FateScape settings except alpha.
# These will overwrite the same variables in BASE_RUNNER if they exist.
MAX_ITER <- 100
REPEAT_TIME <- 10
LABEL_SITE_FRACTION <- 0.70
DROP_DUP_ALPHA <- 1.5
DROP_DUP_BETA  <- 1.5

# Alpha grid.
# alpha1 = LAMBDA_STATE, alpha2 = LAMBDA_BARCODE.
# You can freely add/remove rows.
ALPHA_GRID <- data.frame(
  alpha_setting = c(
    "alpha1_0p00_alpha2_1p00",
    "alpha1_0p05_alpha2_0p95",
    "alpha1_0p10_alpha2_0p90_default",
    "alpha1_0p20_alpha2_0p80",
    "alpha1_0p50_alpha2_0p50",
    "alpha1_0p80_alpha2_0p20",
    "alpha1_1p00_alpha2_0p00"
  ),
  alpha1_state = c(0.00, 0.05, 0.10, 0.20, 0.50, 0.80, 1.00),
  alpha2_barcode = c(1.00, 0.95, 0.90, 0.80, 0.50, 0.20, 0.00),
  stringsAsFactors = FALSE
)

# If TRUE, enforce alpha1 + alpha2 = 1.
CHECK_ALPHA_SUM <- TRUE
ALPHA_SUM_TOL <- 1e-8

# If TRUE, create a basic sensitivity plot if ggplot2 is installed.
MAKE_PLOTS <- TRUE

# -------------------------
# 1. HELPERS
# -------------------------

safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

safe_write_csv <- function(x, file, row.names = FALSE) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, file = file, row.names = row.names)
}

quote_r_string <- function(x) {
  paste0("\"", gsub("\\\\", "/", x), "\"")
}

parse_case_indices <- function(case_indices, case_ids = NULL) {
  n <- if (is.null(case_ids)) NA_integer_ else length(case_ids)

  if (length(case_indices) == 1 && identical(tolower(as.character(case_indices)), "all")) {
    if (is.na(n)) return("all")
    return(seq_len(n))
  }

  if (is.numeric(case_indices)) {
    idx <- as.integer(case_indices)
    if (is.na(n)) return(idx)
    return(idx[idx >= 1 & idx <= n])
  }

  ci <- as.character(case_indices)

  if (length(ci) == 1 && grepl("^\\d+$", ci)) {
    idx <- as.integer(ci)
    if (is.na(n)) return(idx)
    return(idx[idx >= 1 & idx <= n])
  }

  if (length(ci) == 1 && grepl("^\\d+\\s*:\\s*\\d+$", ci)) {
    parts <- strsplit(ci, ":", fixed = TRUE)[[1]]
    idx <- seq.int(as.integer(trimws(parts[1])), as.integer(trimws(parts[2])))
    if (is.na(n)) return(idx)
    return(idx[idx >= 1 & idx <= n])
  }

  if (length(ci) == 1 && grepl(",", ci, fixed = TRUE)) {
    parts <- trimws(strsplit(ci, ",", fixed = TRUE)[[1]])
    if (all(grepl("^\\d+$", parts))) {
      idx <- as.integer(parts)
      if (is.na(n)) return(idx)
      return(idx[idx >= 1 & idx <= n])
    }
  }

  if (!is.null(case_ids)) {
    idx <- match(ci, case_ids)
    return(idx[!is.na(idx)])
  }

  ci
}

replace_assignment <- function(txt, var, value_code, warn_if_missing = TRUE) {
  pat <- paste0("(?m)^", var, "\\s*<-\\s*.*$")
  repl <- paste0(var, " <- ", value_code)
  txt2 <- gsub(pat, repl, txt, perl = TRUE)
  if (identical(txt, txt2) && warn_if_missing) {
    warning("Did not find assignment for ", var, " in BASE_RUNNER.")
  }
  txt2
}

append_alpha_metadata_to_runner <- function(txt, alpha_setting, alpha1_state, alpha2_barcode) {
  # Adds variables used only for metadata. Does not alter algorithm.
  meta <- paste0(
    "\n\n# ---- Injected by alpha-sensitivity wrapper ----\n",
    "ALPHA_SENSITIVITY_MODE <- TRUE\n",
    "ALPHA_SETTING <- ", quote_r_string(alpha_setting), "\n",
    "ALPHA1_STATE <- ", sprintf("%.12g", alpha1_state), "\n",
    "ALPHA2_BARCODE <- ", sprintf("%.12g", alpha2_barcode), "\n"
  )
  paste0(txt, meta)
}

write_patched_runner <- function(base_runner, out_file,
                                 method_input_root,
                                 method_output_root,
                                 case_indices,
                                 alpha_setting,
                                 alpha1_state,
                                 alpha2_barcode) {
  txt <- readLines(base_runner, warn = FALSE, encoding = "UTF-8")
  txt <- paste(txt, collapse = "\n")

  # Paths and batch selection.
  txt <- replace_assignment(txt, "METHOD_INPUT_ROOT", quote_r_string(method_input_root))
  txt <- replace_assignment(txt, "METHOD_OUTPUT_ROOT", quote_r_string(method_output_root))
  txt <- replace_assignment(txt, "CASE_INDICES", quote_r_string(as.character(case_indices)))
  txt <- replace_assignment(txt, "SKIP_FINISHED", ifelse(SKIP_FINISHED, "TRUE", "FALSE"))
  txt <- replace_assignment(txt, "FORCE_RERUN", ifelse(FORCE_RERUN, "TRUE", "FALSE"))
  txt <- replace_assignment(txt, "STOP_ON_ERROR", ifelse(STOP_ON_ERROR, "TRUE", "FALSE"))

  # Keep non-alpha settings fixed.
  txt <- replace_assignment(txt, "LABEL_SITE_FRACTION", as.character(LABEL_SITE_FRACTION), warn_if_missing = FALSE)
  txt <- replace_assignment(txt, "MAX_ITER", as.character(MAX_ITER), warn_if_missing = FALSE)
  txt <- replace_assignment(txt, "REPEAT_TIME", as.character(REPEAT_TIME), warn_if_missing = FALSE)
  txt <- replace_assignment(txt, "DROP_DUP_ALPHA", as.character(DROP_DUP_ALPHA), warn_if_missing = FALSE)
  txt <- replace_assignment(txt, "DROP_DUP_BETA", as.character(DROP_DUP_BETA), warn_if_missing = FALSE)

  # Alpha settings: alpha1 = state, alpha2 = barcode.
  txt <- replace_assignment(txt, "LAMBDA_STATE", sprintf("%.12g", alpha1_state))
  txt <- replace_assignment(txt, "LAMBDA_BARCODE", sprintf("%.12g", alpha2_barcode))

  txt <- append_alpha_metadata_to_runner(txt, alpha_setting, alpha1_state, alpha2_barcode)

  writeLines(txt, out_file, useBytes = TRUE)
  out_file
}

read_if_exists <- function(file) {
  if (!file.exists(file)) return(NULL)
  read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
}

standardize_metric_columns <- function(df) {
  if (is.null(df)) return(NULL)
  if (!"Nye" %in% names(df)) df$Nye <- NA_real_
  if (!"TreeDistance" %in% names(df)) df$TreeDistance <- NA_real_
  if (!"RF_normalized" %in% names(df)) df$RF_normalized <- NA_real_
  if (!"RF_raw" %in% names(df)) df$RF_raw <- NA_real_
  df
}

plot_sensitivity <- function(metrics_all, out_dir) {
  if (!MAKE_PLOTS) return(invisible(NULL))
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Package ggplot2 is not available; skipping alpha-sensitivity plots.")
    return(invisible(NULL))
  }

  metric_cols <- intersect(c("Nye", "TreeDistance", "RF_normalized", "RF_raw"), names(metrics_all))
  if (length(metric_cols) == 0) {
    message("No metric columns available for plotting.")
    return(invisible(NULL))
  }

  metrics_long <- data.frame()
  for (m in metric_cols) {
    tmp <- metrics_all[, intersect(c("case_id", "alpha_setting", "alpha1_state", "alpha2_barcode"), names(metrics_all)), drop = FALSE]
    tmp$metric <- m
    tmp$value <- metrics_all[[m]]
    metrics_long <- rbind(metrics_long, tmp)
  }

  suppressPackageStartupMessages(library(ggplot2))

  p <- ggplot(metrics_long, aes(x = alpha_setting, y = value)) +
    geom_boxplot(outlier.size = 0.6) +
    geom_jitter(width = 0.15, size = 0.7, alpha = 0.5) +
    facet_wrap(~ metric, scales = "free_y") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
    labs(
      x = "Alpha setting",
      y = "Metric value",
      title = "FateScape sensitivity to alpha1/state and alpha2/barcode weights"
    )

  ggplot2::ggsave(file.path(out_dir, "alpha_sensitivity_metrics_boxplot.pdf"), p, width = 10, height = 5)
  ggplot2::ggsave(file.path(out_dir, "alpha_sensitivity_metrics_boxplot.png"), p, width = 10, height = 5, dpi = 300)
}

# -------------------------
# 2. VALIDATE SETTINGS
# -------------------------

if (!file.exists(BASE_RUNNER)) stop("BASE_RUNNER not found: ", BASE_RUNNER)
if (!dir.exists(METHOD_INPUT_ROOT)) stop("METHOD_INPUT_ROOT not found: ", METHOD_INPUT_ROOT)

if (CHECK_ALPHA_SUM) {
  bad <- abs(ALPHA_GRID$alpha1_state + ALPHA_GRID$alpha2_barcode - 1) > ALPHA_SUM_TOL
  if (any(bad)) {
    stop("Some alpha rows do not satisfy alpha1_state + alpha2_barcode = 1: ",
         paste(ALPHA_GRID$alpha_setting[bad], collapse = ", "))
  }
}

SENS_OUTPUT_ROOT <- safe_dir_create(SENS_OUTPUT_ROOT)
RUNNER_DIR <- safe_dir_create(file.path(SENS_OUTPUT_ROOT, "_patched_runners"))
SUMMARY_DIR <- safe_dir_create(file.path(SENS_OUTPUT_ROOT, "_batch_summary"))

safe_write_csv(ALPHA_GRID, file.path(SENS_OUTPUT_ROOT, "alpha_sensitivity_grid.csv"), row.names = FALSE)

message("Base runner       : ", BASE_RUNNER)
message("Method input root : ", METHOD_INPUT_ROOT)
message("Output root       : ", SENS_OUTPUT_ROOT)
message("Case indices      : ", CASE_INDICES)
message("Number of alpha settings: ", nrow(ALPHA_GRID))

# -------------------------
# 3. RUN SENSITIVITY GRID
# -------------------------

all_metric_records <- list()
all_manifest_records <- list()
all_error_records <- list()
run_plan_records <- list()

for (ii in seq_len(nrow(ALPHA_GRID))) {
  alpha_setting <- ALPHA_GRID$alpha_setting[ii]
  alpha1_state <- ALPHA_GRID$alpha1_state[ii]
  alpha2_barcode <- ALPHA_GRID$alpha2_barcode[ii]

  out_root <- file.path(SENS_OUTPUT_ROOT, alpha_setting)
  safe_dir_create(out_root)

  message("\n============================================================")
  message("[FateScape alpha sensitivity] ", alpha_setting)
  message("alpha1/state   = ", alpha1_state)
  message("alpha2/barcode = ", alpha2_barcode)
  message("Output         = ", out_root)

  patched_runner <- file.path(RUNNER_DIR, paste0("run_fatescape_", alpha_setting, ".R"))

  write_patched_runner(
    base_runner = BASE_RUNNER,
    out_file = patched_runner,
    method_input_root = METHOD_INPUT_ROOT,
    method_output_root = out_root,
    case_indices = CASE_INDICES,
    alpha_setting = alpha_setting,
    alpha1_state = alpha1_state,
    alpha2_barcode = alpha2_barcode
  )

  run_plan_records[[length(run_plan_records) + 1L]] <- data.frame(
    alpha_setting = alpha_setting,
    alpha1_state = alpha1_state,
    alpha2_barcode = alpha2_barcode,
    method_input_root = METHOD_INPUT_ROOT,
    method_output_root = out_root,
    patched_runner = patched_runner,
    stringsAsFactors = FALSE
  )

  # Execute the patched v4 runner in an isolated environment.
  run_env <- new.env(parent = globalenv())
  result <- tryCatch(
    {
      source(patched_runner, local = run_env)
      "success"
    },
    error = function(e) {
      msg <- conditionMessage(e)
      message("[ERROR] Alpha setting failed: ", alpha_setting, " | ", msg)
      if (STOP_ON_ERROR) stop(e)
      msg
    }
  )

  # Read per-setting batch outputs.
  metrics_file <- file.path(out_root, "_batch_summary", "FateScape_batch_metrics_all.csv")
  manifest_file <- file.path(out_root, "_batch_summary", "FateScape_batch_manifest.csv")
  errors_file <- file.path(out_root, "_batch_summary", "FateScape_batch_errors.csv")

  metrics <- standardize_metric_columns(read_if_exists(metrics_file))
  manifest <- read_if_exists(manifest_file)
  errors <- read_if_exists(errors_file)

  if (!is.null(metrics)) {
    metrics$alpha_setting <- alpha_setting
    metrics$alpha1_state <- alpha1_state
    metrics$alpha2_barcode <- alpha2_barcode
    metrics$alpha_runner <- patched_runner
    all_metric_records[[length(all_metric_records) + 1L]] <- metrics
  }

  if (!is.null(manifest)) {
    manifest$alpha_setting <- alpha_setting
    manifest$alpha1_state <- alpha1_state
    manifest$alpha2_barcode <- alpha2_barcode
    manifest$alpha_runner <- patched_runner
    all_manifest_records[[length(all_manifest_records) + 1L]] <- manifest
  }

  if (!is.null(errors)) {
    errors$alpha_setting <- alpha_setting
    errors$alpha1_state <- alpha1_state
    errors$alpha2_barcode <- alpha2_barcode
    errors$alpha_runner <- patched_runner
    all_error_records[[length(all_error_records) + 1L]] <- errors
  }

  if (!identical(result, "success")) {
    all_error_records[[length(all_error_records) + 1L]] <- data.frame(
      alpha_setting = alpha_setting,
      alpha1_state = alpha1_state,
      alpha2_barcode = alpha2_barcode,
      status = "runner_failed",
      error = result,
      alpha_runner = patched_runner,
      stringsAsFactors = FALSE
    )
  }
}

# -------------------------
# 4. WRITE COMBINED SUMMARIES
# -------------------------

run_plan <- if (length(run_plan_records) > 0) do.call(rbind, run_plan_records) else data.frame()
safe_write_csv(run_plan, file.path(SUMMARY_DIR, "alpha_sensitivity_run_plan.csv"), row.names = FALSE)

if (length(all_metric_records) > 0) {
  metrics_all <- do.call(rbind, all_metric_records)
  safe_write_csv(metrics_all, file.path(SENS_OUTPUT_ROOT, "alpha_sensitivity_metrics_all.csv"), row.names = FALSE)
  safe_write_csv(metrics_all, file.path(SUMMARY_DIR, "alpha_sensitivity_metrics_all.csv"), row.names = FALSE)
  plot_sensitivity(metrics_all, SENS_OUTPUT_ROOT)
  message("Combined metrics written to: ", file.path(SENS_OUTPUT_ROOT, "alpha_sensitivity_metrics_all.csv"))
} else {
  message("No metrics were collected.")
}

if (length(all_manifest_records) > 0) {
  manifest_all <- do.call(rbind, all_manifest_records)
  safe_write_csv(manifest_all, file.path(SENS_OUTPUT_ROOT, "alpha_sensitivity_manifest_all.csv"), row.names = FALSE)
  safe_write_csv(manifest_all, file.path(SUMMARY_DIR, "alpha_sensitivity_manifest_all.csv"), row.names = FALSE)
  message("Combined manifest written to: ", file.path(SENS_OUTPUT_ROOT, "alpha_sensitivity_manifest_all.csv"))
  print(table(manifest_all$alpha_setting, manifest_all$status, useNA = "ifany"))
} else {
  message("No manifests were collected.")
}

if (length(all_error_records) > 0) {
  errors_all <- do.call(rbind, all_error_records)
  safe_write_csv(errors_all, file.path(SENS_OUTPUT_ROOT, "alpha_sensitivity_errors_all.csv"), row.names = FALSE)
  safe_write_csv(errors_all, file.path(SUMMARY_DIR, "alpha_sensitivity_errors_all.csv"), row.names = FALSE)
  message("Combined errors written to: ", file.path(SENS_OUTPUT_ROOT, "alpha_sensitivity_errors_all.csv"))
}

message("\nFateScape alpha-sensitivity completed.")
message("Output root: ", SENS_OUTPUT_ROOT)
