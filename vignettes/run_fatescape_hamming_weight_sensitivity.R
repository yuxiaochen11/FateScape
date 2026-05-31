#!/usr/bin/env Rscript

# -------------------------
# 0. USER SETTINGS
# -------------------------

BASE_RUNNER <- "FateScape/vignettes/run_fatescape_batch.R"

METHOD_INPUT_ROOT <- "FateScape/data/Sensitivity_analysis"
METHOD_OUTPUT_ROOT <- "FateScape/results/Sensitivity_analysis_hammingweight"

# Optional: if FateScape is not installed as an R package, set package root here.
# Example: FATESCAPE_LOCAL_PATH <- "D:/projects/lineage tracing/FateScape"

CASE_INDICES <- "all"


# Batch behavior passed into each patched base runner.
SKIP_FINISHED <- TRUE
FORCE_RERUN <- FALSE
STOP_ON_ERROR <- FALSE

# Keep ordinary FateScape parameters fixed.
USE_PREIMPUTED_BARCODE <- FALSE
LABEL_SITE_FRACTION <- 0.70
LAMBDA_STATE <- 0.10
LAMBDA_BARCODE <- 0.90
MAX_ITER <- 100
REPEAT_TIME <- 10
DROP_DUP_ALPHA <- 1.5
DROP_DUP_BETA <- 1.5

# If TRUE, create a basic sensitivity plot if ggplot2 is installed.
MAKE_PLOTS <- TRUE

# Hamming-weight sensitivity grid.
# k denotes the Hamming distance class between two child barcodes.
HAMMING_WEIGHT_GRID <- data.frame(
  weight_setting = c(
    "weight_uniform",
    "weight_inverse_k_strict",
    "weight_inverse_sqrt_kplus1",
    "weight_inverse_kplus1_default",
    "weight_inverse_square_kplus1",
    "weight_exponential_mild",
    "weight_exponential_strong"
  ),
  mode = c(
    "uniform",
    "inverse_k",
    "inverse_power",
    "inverse_kplus1",
    "inverse_power",
    "exponential",
    "exponential"
  ),
  gamma = c(0, 1, 0.5, 1, 2, NA, NA),
  offset = c(1, 0, 1, 1, 1, NA, NA),
  tau = c(NA, NA, NA, NA, NA, 0.5, 1.0),
  label = c(
    "w(k)=1",
    "w(k)=1/k; w(0)=1",
    "w(k)=1/sqrt(k+1)",
    "w(k)=1/(k+1), default",
    "w(k)=1/(k+1)^2",
    "w(k)=exp(-0.5k)",
    "w(k)=exp(-k)"
  ),
  stringsAsFactors = FALSE
)

# For the proposal-sampling term, large Hamming distances can be treated as
# low-support rather than rewarded. With N_char=12 this threshold keeps all
# observed distances in the decayed-reward regime, matching the original code
# where diff_sites < 15 was almost always TRUE.
HAMMING_WEIGHT_MAX_REWARD_DISTANCE <- 15
HAMMING_WEIGHT_SCALE <- 5

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
    stop("Cannot find set.seed(GLOBAL_SEED) in BASE_RUNNER; cannot inject Hamming-weight patch safely.")
  }
  idx <- idx[1]
  insert_lines <- strsplit(insert_code, "\n", fixed = TRUE)[[1]]
  paste(c(lines[seq_len(idx)], insert_lines, lines[(idx + 1):length(lines)]), collapse = "\n")
}

make_hamming_weight_patch_code <- function(weight_setting, mode, gamma, offset, tau) {
  gamma_code <- ifelse(is.na(gamma), "NA_real_", sprintf("%.12g", gamma))
  offset_code <- ifelse(is.na(offset), "NA_real_", sprintf("%.12g", offset))
  tau_code <- ifelse(is.na(tau), "NA_real_", sprintf("%.12g", tau))

  paste0('

# ---- Injected by Hamming-weight sensitivity wrapper ----
HAMMING_WEIGHT_SENSITIVITY_MODE <- TRUE
HAMMING_WEIGHT_SETTING <- ', quote_r_string(weight_setting), '
HAMMING_WEIGHT_MODE <- ', quote_r_string(mode), '
HAMMING_WEIGHT_GAMMA <- ', gamma_code, '
HAMMING_WEIGHT_OFFSET <- ', offset_code, '
HAMMING_WEIGHT_TAU <- ', tau_code, '
HAMMING_WEIGHT_MAX_REWARD_DISTANCE <- ', sprintf("%.12g", HAMMING_WEIGHT_MAX_REWARD_DISTANCE), '
HAMMING_WEIGHT_SCALE <- ', sprintf("%.12g", HAMMING_WEIGHT_SCALE), '

patch_fatescape_function <- function(fname, fun) {
  assign(fname, fun, envir = .GlobalEnv)

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

is_missing_barcode_entry_hw <- function(x) {
  if (length(x) == 0) return(TRUE)
  if (is.na(x)) return(TRUE)
  as.character(x) %in% c("", "NA", "NaN", "nan", "NULL", "null", "None", "none", "?", "-", "-1")
}

barcode_entry_equal_hw <- function(a, b) {
  if (is_missing_barcode_entry_hw(a) && is_missing_barcode_entry_hw(b)) return(TRUE)
  if (is_missing_barcode_entry_hw(a) || is_missing_barcode_entry_hw(b)) return(FALSE)

  av <- suppressWarnings(as.numeric(as.character(a)))
  bv <- suppressWarnings(as.numeric(as.character(b)))
  if (!is.na(av) && !is.na(bv)) return(av == bv)

  identical(as.character(a), as.character(b))
}

count_hamming_diff_sites_hw <- function(a, b) {
  n <- min(length(a), length(b))
  if (n == 0) return(0)
  d <- 0

  for (ii in seq_len(n)) {
    am <- is_missing_barcode_entry_hw(a[ii])
    bm <- is_missing_barcode_entry_hw(b[ii])

    # two missing entries contain no information and are not counted as a difference
    if (am && bm) next

    # missing vs observed is treated as a one-site difference for this distance class
    if (xor(am, bm)) {
      d <- d + 1
      next
    }

    if (!barcode_entry_equal_hw(a[ii], b[ii])) d <- d + 1
  }

  if (!is.finite(d) || is.na(d)) d <- Inf
  d
}

hamming_weight_value <- function(k) {
  k <- as.numeric(k)
  out <- rep(NA_real_, length(k))
  valid <- is.finite(k) & !is.na(k)

  if (!any(valid)) return(rep(0, length(k)))

  if (HAMMING_WEIGHT_MODE == "uniform") {
    out[valid] <- 1
  } else if (HAMMING_WEIGHT_MODE == "inverse_k") {
    out[valid] <- ifelse(k[valid] <= 0, 1, 1 / k[valid])
  } else if (HAMMING_WEIGHT_MODE == "inverse_kplus1") {
    out[valid] <- 1 / (k[valid] + 1)
  } else if (HAMMING_WEIGHT_MODE == "inverse_power") {
    offset <- HAMMING_WEIGHT_OFFSET
    gamma <- HAMMING_WEIGHT_GAMMA
    if (!is.finite(offset) || is.na(offset)) offset <- 1
    if (!is.finite(gamma) || is.na(gamma)) gamma <- 1
    out[valid] <- 1 / ((k[valid] + offset)^gamma)
  } else if (HAMMING_WEIGHT_MODE == "exponential") {
    tau <- HAMMING_WEIGHT_TAU
    if (!is.finite(tau) || is.na(tau)) tau <- 1
    out[valid] <- exp(-tau * k[valid])
  } else {
    stop("Unknown HAMMING_WEIGHT_MODE: ", HAMMING_WEIGHT_MODE)
  }

  out[!is.finite(out) | is.na(out)] <- 0
  out
}

node_to_barcode_row_hw <- function(node, tree, barcodes) {
  ntip <- length(tree$tip.label)
  if (length(node) == 0 || is.na(node)) return(NA_integer_)

  if (node > ntip) {
    return(as.integer(node))
  }

  cell_id <- tree$tip.label[node]
  idx <- match(cell_id, rownames(barcodes))
  if (is.na(idx)) return(NA_integer_)
  as.integer(idx)
}

barcode_similarity_likelihood_hweight <- function(tree, N_char, edges, barcodes) {
  barcode_dist <- data.frame(diff_sites = 0:N_char, prob = rep(0, N_char + 1))

  if (is.null(edges) || nrow(edges) == 0) return(1)
  parent_nodes <- unique(edges[, 1])
  parent_nodes <- parent_nodes[!is.na(parent_nodes)]

  n_valid <- 0
  for (par in parent_nodes) {
    children_nodes <- edges[edges[, 1] == par, 2]
    children_nodes <- children_nodes[!is.na(children_nodes)]
    if (length(children_nodes) < 2) next

    child_id_1 <- node_to_barcode_row_hw(children_nodes[1], tree, barcodes)
    child_id_2 <- node_to_barcode_row_hw(children_nodes[2], tree, barcodes)
    if (is.na(child_id_1) || is.na(child_id_2)) next
    if (child_id_1 < 1 || child_id_2 < 1 || child_id_1 > nrow(barcodes) || child_id_2 > nrow(barcodes)) next

    child_barcodes_1 <- barcodes[child_id_1, , drop = TRUE]
    child_barcodes_2 <- barcodes[child_id_2, , drop = TRUE]

    diff_sites_num <- count_hamming_diff_sites_hw(child_barcodes_1, child_barcodes_2)
    if (!is.finite(diff_sites_num) || is.na(diff_sites_num)) next
    diff_sites_num <- max(0, min(N_char, as.integer(round(diff_sites_num))))

    barcode_dist[barcode_dist$diff_sites == diff_sites_num, "prob"] <-
      barcode_dist[barcode_dist$diff_sites == diff_sites_num, "prob"] + 1
    n_valid <- n_valid + 1
  }

  if (n_valid == 0 || sum(barcode_dist$prob) == 0) return(1)

  barcode_dist$prob <- barcode_dist$prob / sum(barcode_dist$prob)
  w <- hamming_weight_value(barcode_dist$diff_sites)
  score <- sum(w * barcode_dist$prob)

  if (!is.finite(score) || is.na(score)) score <- 1
  score
}

nodes_sampling_hweight <- function(tree, edges, cell_state_labels, barcodes, state_lineages) {
  ntip <- length(tree$tip.label)
  total_nodes <- max(c(edges, ntip), na.rm = TRUE)
  if (!is.finite(total_nodes) || total_nodes < 1) total_nodes <- 2 * ntip - 1

  renew_prob <- rep(0, total_nodes)

  if (is.null(edges) || nrow(edges) == 0) {
    renew_prob[seq_len(total_nodes)] <- 1 / total_nodes
    return(renew_prob)
  }

  parent_nodes <- unique(edges[, 1])
  parent_nodes <- parent_nodes[!is.na(parent_nodes)]

  for (par in parent_nodes) {
    children_nodes <- edges[edges[, 1] == par, 2]
    children_nodes <- children_nodes[!is.na(children_nodes)]
    if (length(children_nodes) < 2) next

    child_id_1 <- node_to_barcode_row_hw(children_nodes[1], tree, barcodes)
    child_id_2 <- node_to_barcode_row_hw(children_nodes[2], tree, barcodes)
    if (is.na(child_id_1) || is.na(child_id_2)) next
    if (child_id_1 < 1 || child_id_2 < 1 || child_id_1 > nrow(barcodes) || child_id_2 > nrow(barcodes)) next

    child_barcodes_1 <- barcodes[child_id_1, , drop = TRUE]
    child_barcodes_2 <- barcodes[child_id_2, , drop = TRUE]

    parent_state <- if (par <= length(cell_state_labels)) cell_state_labels[par] else NA
    child_state_1 <- if (child_id_1 <= length(cell_state_labels)) cell_state_labels[child_id_1] else NA
    child_state_2 <- if (child_id_2 <= length(cell_state_labels)) cell_state_labels[child_id_2] else NA

    if (is.na(parent_state) || identical(parent_state, 0) || identical(as.character(parent_state), "0")) {
      renew_prob[children_nodes[1]] <- renew_prob[children_nodes[1]] + 10
      renew_prob[children_nodes[2]] <- renew_prob[children_nodes[2]] + 10
    } else {
      for (lineage in state_lineages) {
        if ((parent_state %in% lineage) && (child_state_1 %in% lineage)) {
          state_dist_1 <- match(child_state_1, lineage) - match(parent_state, lineage)
          if (!is.na(state_dist_1) && state_dist_1 > 0 && state_dist_1 < 4) {
            renew_prob[children_nodes[1]] <- renew_prob[children_nodes[1]] - (5 / state_dist_1)
          } else {
            renew_prob[children_nodes[1]] <- renew_prob[children_nodes[1]] + 5
          }
        }

        if ((parent_state %in% lineage) && (child_state_2 %in% lineage)) {
          state_dist_2 <- match(child_state_2, lineage) - match(parent_state, lineage)
          if (!is.na(state_dist_2) && state_dist_2 > 0 && state_dist_2 < 4) {
            renew_prob[children_nodes[2]] <- renew_prob[children_nodes[2]] - (5 / state_dist_2)
          } else {
            renew_prob[children_nodes[2]] <- renew_prob[children_nodes[2]] + 5
          }
        }
      }
    }

    diff_sites <- count_hamming_diff_sites_hw(child_barcodes_1, child_barcodes_2)
    if (!is.finite(diff_sites) || is.na(diff_sites)) next

    if (diff_sites < HAMMING_WEIGHT_MAX_REWARD_DISTANCE) {
      reward <- HAMMING_WEIGHT_SCALE * hamming_weight_value(diff_sites)
      renew_prob[children_nodes[1]] <- renew_prob[children_nodes[1]] - reward
      renew_prob[children_nodes[2]] <- renew_prob[children_nodes[2]] - reward
    } else {
      renew_prob[children_nodes[1]] <- renew_prob[children_nodes[1]] + 5
      renew_prob[children_nodes[2]] <- renew_prob[children_nodes[2]] + 5
    }
  }

  renew_prob[!is.finite(renew_prob) | is.na(renew_prob)] <- 0
  renew_prob[renew_prob < 0] <- 0

  if (sum(renew_prob) <= 0) {
    valid_nodes <- unique(as.integer(edges[, 2]))
    valid_nodes <- valid_nodes[!is.na(valid_nodes) & valid_nodes >= 1 & valid_nodes <= length(renew_prob)]
    if (length(valid_nodes) == 0) valid_nodes <- seq_along(renew_prob)
    renew_prob[valid_nodes] <- 1 / length(valid_nodes)
  } else {
    renew_prob <- renew_prob / sum(renew_prob)
  }

  renew_prob
}

patch_fatescape_function("barcode_similarity_likelihood", barcode_similarity_likelihood_hweight)
patch_fatescape_function("nodes_sampling", nodes_sampling_hweight)

message("[Hamming-weight sensitivity] Patched barcode_similarity_likelihood() and nodes_sampling().")
message("[Hamming-weight sensitivity] setting=", HAMMING_WEIGHT_SETTING,
        "; mode=", HAMMING_WEIGHT_MODE,
        "; gamma=", HAMMING_WEIGHT_GAMMA,
        "; offset=", HAMMING_WEIGHT_OFFSET,
        "; tau=", HAMMING_WEIGHT_TAU)
')
}

write_patched_runner <- function(base_runner, out_file,
                                 method_input_root,
                                 method_output_root,
                                 case_indices,
                                 weight_setting, mode, gamma, offset, tau) {
  txt <- readLines(base_runner, warn = FALSE, encoding = "UTF-8")
  txt <- paste(txt, collapse = "\n")

  txt <- replace_assignment(txt, "METHOD_INPUT_ROOT", quote_r_string(method_input_root))
  txt <- replace_assignment(txt, "METHOD_OUTPUT_ROOT", quote_r_string(method_output_root))
  txt <- replace_assignment(txt, "CASE_INDICES", quote_r_string(as.character(case_indices)))
  txt <- replace_assignment(txt, "SKIP_FINISHED", ifelse(SKIP_FINISHED, "TRUE", "FALSE"))
  txt <- replace_assignment(txt, "FORCE_RERUN", ifelse(FORCE_RERUN, "TRUE", "FALSE"))
  txt <- replace_assignment(txt, "STOP_ON_ERROR", ifelse(STOP_ON_ERROR, "TRUE", "FALSE"))

  txt <- replace_assignment(txt, "USE_PREIMPUTED_BARCODE", ifelse(USE_PREIMPUTED_BARCODE, "TRUE", "FALSE"), warn_if_missing = FALSE)
  txt <- replace_assignment(txt, "LABEL_SITE_FRACTION", as.character(LABEL_SITE_FRACTION), warn_if_missing = FALSE)
  txt <- replace_assignment(txt, "LAMBDA_STATE", sprintf("%.12g", LAMBDA_STATE), warn_if_missing = FALSE)
  txt <- replace_assignment(txt, "LAMBDA_BARCODE", sprintf("%.12g", LAMBDA_BARCODE), warn_if_missing = FALSE)
  txt <- replace_assignment(txt, "MAX_ITER", as.character(MAX_ITER), warn_if_missing = FALSE)
  txt <- replace_assignment(txt, "REPEAT_TIME", as.character(REPEAT_TIME), warn_if_missing = FALSE)
  txt <- replace_assignment(txt, "DROP_DUP_ALPHA", as.character(DROP_DUP_ALPHA), warn_if_missing = FALSE)
  txt <- replace_assignment(txt, "DROP_DUP_BETA", as.character(DROP_DUP_BETA), warn_if_missing = FALSE)

  txt <- insert_after_set_seed(
    txt,
    make_hamming_weight_patch_code(
      weight_setting = weight_setting,
      mode = mode,
      gamma = gamma,
      offset = offset,
      tau = tau
    )
  )

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
    message("Package ggplot2 is not available; skipping Hamming-weight sensitivity plots.")
    return(invisible(NULL))
  }

  metric_cols <- intersect(c("Nye", "TreeDistance", "RF_normalized", "RF_raw"), names(metrics_all))
  if (length(metric_cols) == 0) {
    message("No metric columns available for plotting.")
    return(invisible(NULL))
  }

  metrics_long <- data.frame()
  for (m in metric_cols) {
    id_cols <- intersect(c("case_id", "weight_setting", "weight_label"), names(metrics_all))
    tmp <- metrics_all[, id_cols, drop = FALSE]
    tmp$metric <- m
    tmp$value <- metrics_all[[m]]
    metrics_long <- rbind(metrics_long, tmp)
  }

  suppressPackageStartupMessages(library(ggplot2))
  metrics_long$weight_setting <- factor(metrics_long$weight_setting, levels = HAMMING_WEIGHT_GRID$weight_setting)

  p <- ggplot(metrics_long, aes(x = weight_setting, y = value)) +
    geom_boxplot(outlier.size = 0.6) +
    geom_jitter(width = 0.15, size = 0.7, alpha = 0.5) +
    facet_wrap(~ metric, scales = "free_y") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
    labs(
      x = "Hamming-distance weight setting",
      y = "Metric value",
      title = "FateScape sensitivity to Hamming-distance class weighting"
    )

  ggplot2::ggsave(file.path(out_dir, "hamming_weight_sensitivity_metrics_boxplot.pdf"), p, width = 11, height = 5)
  ggplot2::ggsave(file.path(out_dir, "hamming_weight_sensitivity_metrics_boxplot.png"), p, width = 11, height = 5, dpi = 300)
}

# -------------------------
# 2. VALIDATE SETTINGS
# -------------------------

if (!file.exists(BASE_RUNNER)) stop("BASE_RUNNER not found: ", BASE_RUNNER)
if (!dir.exists(METHOD_INPUT_ROOT)) stop("METHOD_INPUT_ROOT not found: ", METHOD_INPUT_ROOT)

SENS_OUTPUT_ROOT <- safe_dir_create(SENS_OUTPUT_ROOT)
RUNNER_DIR <- safe_dir_create(file.path(SENS_OUTPUT_ROOT, "_patched_runners"))
SUMMARY_DIR <- safe_dir_create(file.path(SENS_OUTPUT_ROOT, "_batch_summary"))

safe_write_csv(HAMMING_WEIGHT_GRID, file.path(SENS_OUTPUT_ROOT, "hamming_weight_sensitivity_grid.csv"), row.names = FALSE)

message("Base runner       : ", BASE_RUNNER)
message("Method input root : ", METHOD_INPUT_ROOT)
message("Output root       : ", SENS_OUTPUT_ROOT)
message("Case indices      : ", CASE_INDICES)
message("Number of Hamming-weight settings: ", nrow(HAMMING_WEIGHT_GRID))

# -------------------------
# 3. RUN SENSITIVITY GRID
# -------------------------

all_metric_records <- list()
all_manifest_records <- list()
all_error_records <- list()
run_plan_records <- list()

for (ii in seq_len(nrow(HAMMING_WEIGHT_GRID))) {
  weight_setting <- HAMMING_WEIGHT_GRID$weight_setting[ii]
  mode <- HAMMING_WEIGHT_GRID$mode[ii]
  gamma <- HAMMING_WEIGHT_GRID$gamma[ii]
  offset <- HAMMING_WEIGHT_GRID$offset[ii]
  tau <- HAMMING_WEIGHT_GRID$tau[ii]
  weight_label <- HAMMING_WEIGHT_GRID$label[ii]

  out_root <- file.path(SENS_OUTPUT_ROOT, weight_setting)
  safe_dir_create(out_root)

  message("\n============================================================")
  message("[FateScape Hamming-weight sensitivity] ", weight_setting)
  message("label  = ", weight_label)
  message("mode   = ", mode, "; gamma=", gamma, "; offset=", offset, "; tau=", tau)
  message("Output = ", out_root)

  patched_runner <- file.path(RUNNER_DIR, paste0("run_fatescape_", weight_setting, ".R"))

  write_patched_runner(
    base_runner = BASE_RUNNER,
    out_file = patched_runner,
    method_input_root = METHOD_INPUT_ROOT,
    method_output_root = out_root,
    case_indices = CASE_INDICES,
    weight_setting = weight_setting,
    mode = mode,
    gamma = gamma,
    offset = offset,
    tau = tau
  )

  run_plan_records[[length(run_plan_records) + 1L]] <- data.frame(
    weight_setting = weight_setting,
    weight_label = weight_label,
    mode = mode,
    gamma = gamma,
    offset = offset,
    tau = tau,
    method_input_root = METHOD_INPUT_ROOT,
    method_output_root = out_root,
    patched_runner = patched_runner,
    stringsAsFactors = FALSE
  )

  run_env <- new.env(parent = globalenv())
  result <- tryCatch(
    {
      source(patched_runner, local = run_env)
      "success"
    },
    error = function(e) {
      msg <- conditionMessage(e)
      message("[ERROR] Hamming-weight setting failed: ", weight_setting, " | ", msg)
      if (STOP_ON_ERROR) stop(e)
      msg
    }
  )

  metrics_file <- file.path(out_root, "_batch_summary", "FateScape_batch_metrics_all.csv")
  manifest_file <- file.path(out_root, "_batch_summary", "FateScape_batch_manifest.csv")
  errors_file <- file.path(out_root, "_batch_summary", "FateScape_batch_errors.csv")

  metrics <- standardize_metric_columns(read_if_exists(metrics_file))
  manifest <- read_if_exists(manifest_file)
  errors <- read_if_exists(errors_file)

  if (!is.null(metrics)) {
    metrics$weight_setting <- weight_setting
    metrics$weight_label <- weight_label
    metrics$hamming_weight_mode <- mode
    metrics$hamming_weight_gamma <- gamma
    metrics$hamming_weight_offset <- offset
    metrics$hamming_weight_tau <- tau
    metrics$hamming_weight_runner <- patched_runner
    all_metric_records[[length(all_metric_records) + 1L]] <- metrics
  }

  if (!is.null(manifest)) {
    manifest$weight_setting <- weight_setting
    manifest$weight_label <- weight_label
    manifest$hamming_weight_mode <- mode
    manifest$hamming_weight_gamma <- gamma
    manifest$hamming_weight_offset <- offset
    manifest$hamming_weight_tau <- tau
    manifest$hamming_weight_runner <- patched_runner
    all_manifest_records[[length(all_manifest_records) + 1L]] <- manifest
  }

  if (!is.null(errors)) {
    errors$weight_setting <- weight_setting
    errors$weight_label <- weight_label
    errors$hamming_weight_mode <- mode
    errors$hamming_weight_gamma <- gamma
    errors$hamming_weight_offset <- offset
    errors$hamming_weight_tau <- tau
    errors$hamming_weight_runner <- patched_runner
    all_error_records[[length(all_error_records) + 1L]] <- errors
  }

  if (!identical(result, "success")) {
    all_error_records[[length(all_error_records) + 1L]] <- data.frame(
      weight_setting = weight_setting,
      weight_label = weight_label,
      hamming_weight_mode = mode,
      hamming_weight_gamma = gamma,
      hamming_weight_offset = offset,
      hamming_weight_tau = tau,
      status = "runner_failed",
      error = result,
      hamming_weight_runner = patched_runner,
      stringsAsFactors = FALSE
    )
  }
}

# -------------------------
# 4. WRITE COMBINED SUMMARIES
# -------------------------

run_plan <- if (length(run_plan_records) > 0) do.call(rbind, run_plan_records) else data.frame()
safe_write_csv(run_plan, file.path(SUMMARY_DIR, "hamming_weight_sensitivity_run_plan.csv"), row.names = FALSE)

if (length(all_metric_records) > 0) {
  metrics_all <- do.call(rbind, all_metric_records)
  safe_write_csv(metrics_all, file.path(SENS_OUTPUT_ROOT, "hamming_weight_sensitivity_metrics_all.csv"), row.names = FALSE)
  safe_write_csv(metrics_all, file.path(SUMMARY_DIR, "hamming_weight_sensitivity_metrics_all.csv"), row.names = FALSE)
  plot_sensitivity(metrics_all, SENS_OUTPUT_ROOT)
  message("Combined metrics written to: ", file.path(SENS_OUTPUT_ROOT, "hamming_weight_sensitivity_metrics_all.csv"))
} else {
  message("No metrics were collected.")
}

if (length(all_manifest_records) > 0) {
  manifest_all <- do.call(rbind, all_manifest_records)
  safe_write_csv(manifest_all, file.path(SENS_OUTPUT_ROOT, "hamming_weight_sensitivity_manifest_all.csv"), row.names = FALSE)
  safe_write_csv(manifest_all, file.path(SUMMARY_DIR, "hamming_weight_sensitivity_manifest_all.csv"), row.names = FALSE)
  message("Combined manifest written to: ", file.path(SENS_OUTPUT_ROOT, "hamming_weight_sensitivity_manifest_all.csv"))
  print(table(manifest_all$weight_setting, manifest_all$status, useNA = "ifany"))
} else {
  message("No manifests were collected.")
}

if (length(all_error_records) > 0) {
  errors_all <- do.call(rbind, all_error_records)
  safe_write_csv(errors_all, file.path(SENS_OUTPUT_ROOT, "hamming_weight_sensitivity_errors_all.csv"), row.names = FALSE)
  safe_write_csv(errors_all, file.path(SUMMARY_DIR, "hamming_weight_sensitivity_errors_all.csv"), row.names = FALSE)
  message("Combined errors written to: ", file.path(SENS_OUTPUT_ROOT, "hamming_weight_sensitivity_errors_all.csv"))
}

message("\nFateScape Hamming-weight sensitivity completed.")
message("Output root: ", SENS_OUTPUT_ROOT)
