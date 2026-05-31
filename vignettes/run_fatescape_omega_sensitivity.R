#!/usr/bin/env Rscript

# -------------------------
# 0. USER SETTINGS
# -------------------------
BASE_RUNNER <- "FateScape/vignettes/run_fatescape_batch.R"

METHOD_INPUT_ROOT <- "FateScape/data/Sensitivity_analysis"
METHOD_OUTPUT_ROOT <- "FateScape/results/Sensitivity_analysis_omega"


CASE_INDICES <- "all"


# Batch behavior passed into each patched base runner.
SKIP_FINISHED <- TRUE
FORCE_RERUN <- FALSE
STOP_ON_ERROR <- FALSE

# Keep the same core FateScape settings except omega.
# These will overwrite the same variables in BASE_RUNNER if they exist.
USE_PREIMPUTED_BARCODE <- FALSE
LABEL_SITE_FRACTION <- 0.70
LAMBDA_STATE <- 0.10
LAMBDA_BARCODE <- 0.90
MAX_ITER <- 100
REPEAT_TIME <- 10
DROP_DUP_ALPHA <- 1.5

# Omega grid.
# omega = penalty for conflicting nonzero mutations at the same barcode site.
# omega = 1 is Hamming-like; omega = 1.5 matches the current default beta.
OMEGA_GRID <- data.frame(
  omega_setting = c(
    "omega_0p50_weak_conflict",
    "omega_1p00_hamming_like",
    "omega_1p50_default",
    "omega_2p00_strong_conflict",
    "omega_3p00_very_strong_conflict"
  ),
  omega = c(0.5, 1.0, 1.5, 2.0, 3.0),
  stringsAsFactors = FALSE
)

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

replace_assignment <- function(txt, var, value_code, warn_if_missing = TRUE) {
  pat <- paste0("(?m)^", var, "\\s*<-\\s*.*$")
  repl <- paste0(var, " <- ", value_code)
  txt2 <- gsub(pat, repl, txt, perl = TRUE)
  if (identical(txt, txt2) && warn_if_missing) {
    warning("Did not find assignment for ", var, " in BASE_RUNNER.")
  }
  txt2
}

insert_after_set_seed <- function(txt, insert_code) {
  lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  idx <- grep("^set\\.seed\\(GLOBAL_SEED\\)", lines)
  if (length(idx) == 0) {
    stop("Cannot find set.seed(GLOBAL_SEED) in BASE_RUNNER; cannot inject omega patch safely.")
  }
  idx <- idx[1]
  insert_lines <- strsplit(insert_code, "\n", fixed = TRUE)[[1]]
  paste(c(lines[seq_len(idx)], insert_lines, lines[(idx + 1):length(lines)]), collapse = "\n")
}

make_omega_patch_code <- function(omega_setting, omega_value) {
  paste0('

# ---- Injected by omega-sensitivity wrapper ----
OMEGA_SENSITIVITY_MODE <- TRUE
OMEGA_SETTING <- ', quote_r_string(omega_setting), '
OMEGA_CONFLICT_PENALTY <- ', sprintf("%.12g", omega_value), '

patch_fatescape_function <- function(fname, fun) {
  # Make patched function visible in the runner environment.
  assign(fname, fun, envir = .GlobalEnv)

  # Patch the FateScape namespace if the function exists there.
  if ("FateScape" %in% loadedNamespaces()) {
    ns <- asNamespace("FateScape")
    if (exists(fname, envir = ns, inherits = FALSE)) {
      was_locked <- bindingIsLocked(fname, ns)
      if (was_locked) unlockBinding(fname, ns)
      assign(fname, fun, envir = ns)
      if (was_locked) lockBinding(fname, ns)
    }
  }

  invisible(TRUE)
}

is_missing_barcode_entry_omega <- function(x) {
  if (length(x) == 0) return(TRUE)
  if (is.na(x)) return(TRUE)
  as.character(x) %in% c("", "NA", "NaN", "nan", "NULL", "null", "None", "none", "?", "-", "-1")
}

barcode_entry_numeric_omega <- function(x) {
  if (is_missing_barcode_entry_omega(x)) return(NA_real_)
  suppressWarnings(as.numeric(as.character(x)))
}

omega_pair_distance <- function(a, b, missing_penalty = 1) {
  n <- min(length(a), length(b))
  if (n == 0) return(0)
  d <- 0

  for (ii in seq_len(n)) {
    am <- is_missing_barcode_entry_omega(a[ii])
    bm <- is_missing_barcode_entry_omega(b[ii])

    if (am && bm) next
    if (xor(am, bm)) {
      d <- d + missing_penalty
      next
    }

    av <- barcode_entry_numeric_omega(a[ii])
    bv <- barcode_entry_numeric_omega(b[ii])

    if (is.na(av) || is.na(bv)) {
      if (!identical(as.character(a[ii]), as.character(b[ii]))) {
        d <- d + missing_penalty
      }
      next
    }

    if (av == bv) {
      next
    } else if (av == 0 || bv == 0) {
      d <- d + 1
    } else {
      d <- d + OMEGA_CONFLICT_PENALTY
    }
  }

  if (!is.finite(d) || is.na(d)) d <- Inf
  d
}

hamming_distance_omega <- function(barcodes) {
  barcodes <- as.data.frame(barcodes, check.names = FALSE)
  n <- nrow(barcodes)
  if (n < 2) stop("hamming_distance_omega requires at least two barcodes.")

  rn <- rownames(barcodes)
  if (is.null(rn) || any(is.na(rn)) || any(!nzchar(rn))) {
    rn <- paste0("cell_", seq_len(n))
  }

  d <- matrix(0, n, n, dimnames = list(rn, rn))

  if (n >= 2) {
    for (ii in seq_len(n - 1L)) {
      for (jj in (ii + 1L):n) {
        dij <- omega_pair_distance(
          barcodes[ii, , drop = TRUE],
          barcodes[jj, , drop = TRUE],
          missing_penalty = 1
        )
        d[ii, jj] <- dij
        d[jj, ii] <- dij
      }
    }
  }

  stats::as.dist(d)
}

improved_hamming_distance_omega <- function(barcode1, barcode2, alpha = 0, beta = OMEGA_CONFLICT_PENALTY, N_char = NULL) {
  if (is.null(N_char)) N_char <- min(length(barcode1), length(barcode2))
  if (N_char <= 0) return(0)

  barcode1 <- barcode1[seq_len(N_char)]
  barcode2 <- barcode2[seq_len(N_char)]
  d <- 0

  for (ii in seq_len(N_char)) {
    am <- is_missing_barcode_entry_omega(barcode1[ii])
    bm <- is_missing_barcode_entry_omega(barcode2[ii])

    if (am && bm) next
    if (xor(am, bm)) {
      d <- d + beta
      next
    }

    av <- barcode_entry_numeric_omega(barcode1[ii])
    bv <- barcode_entry_numeric_omega(barcode2[ii])

    if (is.na(av) || is.na(bv)) {
      s1 <- as.character(barcode1[ii])
      s2 <- as.character(barcode2[ii])
      if (identical(s1, s2) && !(s1 %in% c("0", "0.0"))) {
        d <- d - alpha
      } else if (!identical(s1, s2)) {
        d <- d + beta
      }
      next
    }

    if (av == bv && av != 0) {
      d <- d - alpha
    } else if (av != bv && (av == 0 || bv == 0)) {
      d <- d + 1
    } else if (av != bv && av != 0 && bv != 0) {
      d <- d + beta
    }
  }

  if (!is.finite(d) || is.na(d)) d <- Inf
  d
}

# Patch the two barcode-distance functions that influence initialization and
# duplicate-tip resolution. Other FateScape workflow calls are unchanged.
patch_fatescape_function("hamming_distance", hamming_distance_omega)
patch_fatescape_function("improved_hamming_distance", improved_hamming_distance_omega)

message("[Omega sensitivity] Patched barcode distances with omega=", OMEGA_CONFLICT_PENALTY)
')
}

write_patched_runner <- function(base_runner, out_file,
                                 method_input_root,
                                 method_output_root,
                                 case_indices,
                                 omega_setting,
                                 omega_value) {
  txt <- readLines(base_runner, warn = FALSE, encoding = "UTF-8")
  txt <- paste(txt, collapse = "\n")

  # Paths and batch selection.
  txt <- replace_assignment(txt, "METHOD_INPUT_ROOT", quote_r_string(method_input_root))
  txt <- replace_assignment(txt, "METHOD_OUTPUT_ROOT", quote_r_string(method_output_root))
  txt <- replace_assignment(txt, "CASE_INDICES", quote_r_string(as.character(case_indices)))
  txt <- replace_assignment(txt, "SKIP_FINISHED", ifelse(SKIP_FINISHED, "TRUE", "FALSE"))
  txt <- replace_assignment(txt, "FORCE_RERUN", ifelse(FORCE_RERUN, "TRUE", "FALSE"))
  txt <- replace_assignment(txt, "STOP_ON_ERROR", ifelse(STOP_ON_ERROR, "TRUE", "FALSE"))

  # Keep non-omega settings fixed.
  txt <- replace_assignment(txt, "USE_PREIMPUTED_BARCODE", ifelse(USE_PREIMPUTED_BARCODE, "TRUE", "FALSE"), warn_if_missing = FALSE)
  txt <- replace_assignment(txt, "LABEL_SITE_FRACTION", as.character(LABEL_SITE_FRACTION), warn_if_missing = FALSE)
  txt <- replace_assignment(txt, "LAMBDA_STATE", sprintf("%.12g", LAMBDA_STATE), warn_if_missing = FALSE)
  txt <- replace_assignment(txt, "LAMBDA_BARCODE", sprintf("%.12g", LAMBDA_BARCODE), warn_if_missing = FALSE)
  txt <- replace_assignment(txt, "MAX_ITER", as.character(MAX_ITER), warn_if_missing = FALSE)
  txt <- replace_assignment(txt, "REPEAT_TIME", as.character(REPEAT_TIME), warn_if_missing = FALSE)
  txt <- replace_assignment(txt, "DROP_DUP_ALPHA", as.character(DROP_DUP_ALPHA), warn_if_missing = FALSE)

  # In the original FateScape code, beta is the nonzero-conflict penalty used by
  # improved_hamming_distance in duplicate-tip resolution. Here beta is omega.
  txt <- replace_assignment(txt, "DROP_DUP_BETA", sprintf("%.12g", omega_value), warn_if_missing = FALSE)

  # Inject omega-specific distance functions after FateScape has been loaded.
  txt <- insert_after_set_seed(txt, make_omega_patch_code(omega_setting, omega_value))

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
  if (!"eval_error" %in% names(df)) df$eval_error <- NA_character_
  df
}

plot_sensitivity <- function(metrics_all, out_dir) {
  if (!MAKE_PLOTS) return(invisible(NULL))
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Package ggplot2 is not available; skipping omega-sensitivity plots.")
    return(invisible(NULL))
  }

  metric_cols <- intersect(c("Nye", "TreeDistance", "RF_normalized", "RF_raw"), names(metrics_all))
  if (length(metric_cols) == 0) {
    message("No metric columns available for plotting.")
    return(invisible(NULL))
  }

  metrics_long <- data.frame()
  for (m in metric_cols) {
    id_cols <- intersect(c("case_id", "omega_setting", "omega"), names(metrics_all))
    tmp <- metrics_all[, id_cols, drop = FALSE]
    tmp$metric <- m
    tmp$value <- metrics_all[[m]]
    metrics_long <- rbind(metrics_long, tmp)
  }

  suppressPackageStartupMessages(library(ggplot2))

  metrics_long$omega_setting <- factor(metrics_long$omega_setting, levels = OMEGA_GRID$omega_setting)

  p <- ggplot(metrics_long, aes(x = omega_setting, y = value)) +
    geom_boxplot(outlier.size = 0.6) +
    geom_jitter(width = 0.15, size = 0.7, alpha = 0.5) +
    facet_wrap(~ metric, scales = "free_y") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
    labs(
      x = "Omega setting",
      y = "Metric value",
      title = "FateScape sensitivity to omega conflict penalty"
    )

  ggplot2::ggsave(file.path(out_dir, "omega_sensitivity_metrics_boxplot.pdf"), p, width = 10, height = 5)
  ggplot2::ggsave(file.path(out_dir, "omega_sensitivity_metrics_boxplot.png"), p, width = 10, height = 5, dpi = 300)
}

# -------------------------
# 2. VALIDATE SETTINGS
# -------------------------

if (!file.exists(BASE_RUNNER)) stop("BASE_RUNNER not found: ", BASE_RUNNER)
if (!dir.exists(METHOD_INPUT_ROOT)) stop("METHOD_INPUT_ROOT not found: ", METHOD_INPUT_ROOT)
if (any(!is.finite(OMEGA_GRID$omega) | OMEGA_GRID$omega <= 0)) {
  stop("OMEGA_GRID$omega must contain positive finite values.")
}

SENS_OUTPUT_ROOT <- safe_dir_create(SENS_OUTPUT_ROOT)
RUNNER_DIR <- safe_dir_create(file.path(SENS_OUTPUT_ROOT, "_patched_runners"))
SUMMARY_DIR <- safe_dir_create(file.path(SENS_OUTPUT_ROOT, "_batch_summary"))

safe_write_csv(OMEGA_GRID, file.path(SENS_OUTPUT_ROOT, "omega_sensitivity_grid.csv"), row.names = FALSE)

message("Base runner       : ", BASE_RUNNER)
message("Method input root : ", METHOD_INPUT_ROOT)
message("Output root       : ", SENS_OUTPUT_ROOT)
message("Case indices      : ", CASE_INDICES)
message("Number of omega settings: ", nrow(OMEGA_GRID))

# -------------------------
# 3. RUN SENSITIVITY GRID
# -------------------------

all_metric_records <- list()
all_manifest_records <- list()
all_error_records <- list()
run_plan_records <- list()

for (ii in seq_len(nrow(OMEGA_GRID))) {
  omega_setting <- OMEGA_GRID$omega_setting[ii]
  omega_value <- OMEGA_GRID$omega[ii]

  out_root <- file.path(SENS_OUTPUT_ROOT, omega_setting)
  safe_dir_create(out_root)

  message("\n============================================================")
  message("[FateScape omega sensitivity] ", omega_setting)
  message("omega = ", omega_value)
  message("Output = ", out_root)

  patched_runner <- file.path(RUNNER_DIR, paste0("run_fatescape_", omega_setting, ".R"))

  write_patched_runner(
    base_runner = BASE_RUNNER,
    out_file = patched_runner,
    method_input_root = METHOD_INPUT_ROOT,
    method_output_root = out_root,
    case_indices = CASE_INDICES,
    omega_setting = omega_setting,
    omega_value = omega_value
  )

  run_plan_records[[length(run_plan_records) + 1L]] <- data.frame(
    omega_setting = omega_setting,
    omega = omega_value,
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
      message("[ERROR] Omega setting failed: ", omega_setting, " | ", msg)
      if (STOP_ON_ERROR) stop(e)
      msg
    }
  )

  # Read per-setting batch outputs from the v4 base runner.
  metrics_file <- file.path(out_root, "_batch_summary", "FateScape_batch_metrics_all.csv")
  manifest_file <- file.path(out_root, "_batch_summary", "FateScape_batch_manifest.csv")
  errors_file <- file.path(out_root, "_batch_summary", "FateScape_batch_errors.csv")

  metrics <- standardize_metric_columns(read_if_exists(metrics_file))
  manifest <- read_if_exists(manifest_file)
  errors <- read_if_exists(errors_file)

  if (!is.null(metrics)) {
    metrics$omega_setting <- omega_setting
    metrics$omega <- omega_value
    metrics$omega_runner <- patched_runner
    all_metric_records[[length(all_metric_records) + 1L]] <- metrics
  }

  if (!is.null(manifest)) {
    manifest$omega_setting <- omega_setting
    manifest$omega <- omega_value
    manifest$omega_runner <- patched_runner
    all_manifest_records[[length(all_manifest_records) + 1L]] <- manifest
  }

  if (!is.null(errors)) {
    errors$omega_setting <- omega_setting
    errors$omega <- omega_value
    errors$omega_runner <- patched_runner
    all_error_records[[length(all_error_records) + 1L]] <- errors
  }

  if (!identical(result, "success")) {
    all_error_records[[length(all_error_records) + 1L]] <- data.frame(
      omega_setting = omega_setting,
      omega = omega_value,
      status = "runner_failed",
      error = result,
      omega_runner = patched_runner,
      stringsAsFactors = FALSE
    )
  }
}

# -------------------------
# 4. WRITE COMBINED SUMMARIES
# -------------------------

run_plan <- if (length(run_plan_records) > 0) do.call(rbind, run_plan_records) else data.frame()
safe_write_csv(run_plan, file.path(SUMMARY_DIR, "omega_sensitivity_run_plan.csv"), row.names = FALSE)

if (length(all_metric_records) > 0) {
  metrics_all <- do.call(rbind, all_metric_records)
  safe_write_csv(metrics_all, file.path(SENS_OUTPUT_ROOT, "omega_sensitivity_metrics_all.csv"), row.names = FALSE)
  safe_write_csv(metrics_all, file.path(SUMMARY_DIR, "omega_sensitivity_metrics_all.csv"), row.names = FALSE)
  plot_sensitivity(metrics_all, SENS_OUTPUT_ROOT)
  message("Combined metrics written to: ", file.path(SENS_OUTPUT_ROOT, "omega_sensitivity_metrics_all.csv"))
} else {
  message("No metrics were collected.")
}

if (length(all_manifest_records) > 0) {
  manifest_all <- do.call(rbind, all_manifest_records)
  safe_write_csv(manifest_all, file.path(SENS_OUTPUT_ROOT, "omega_sensitivity_manifest_all.csv"), row.names = FALSE)
  safe_write_csv(manifest_all, file.path(SUMMARY_DIR, "omega_sensitivity_manifest_all.csv"), row.names = FALSE)
  message("Combined manifest written to: ", file.path(SENS_OUTPUT_ROOT, "omega_sensitivity_manifest_all.csv"))
  print(table(manifest_all$omega_setting, manifest_all$status, useNA = "ifany"))
} else {
  message("No manifests were collected.")
}

if (length(all_error_records) > 0) {
  errors_all <- do.call(rbind, all_error_records)
  safe_write_csv(errors_all, file.path(SENS_OUTPUT_ROOT, "omega_sensitivity_errors_all.csv"), row.names = FALSE)
  safe_write_csv(errors_all, file.path(SUMMARY_DIR, "omega_sensitivity_errors_all.csv"), row.names = FALSE)
  message("Combined errors written to: ", file.path(SENS_OUTPUT_ROOT, "omega_sensitivity_errors_all.csv"))
}

message("\nFateScape omega-sensitivity completed.")
message("Output root: ", SENS_OUTPUT_ROOT)
