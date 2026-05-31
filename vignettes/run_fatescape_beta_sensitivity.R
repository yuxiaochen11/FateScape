#!/usr/bin/env Rscript



# -------------------------
# 0. USER SETTINGS
# -------------------------
BASE_RUNNER <- "FateScape/vignettes/run_fatescape_batch.R"

METHOD_INPUT_ROOT <- "FateScape/data/Sensitivity_analysis"
METHOD_OUTPUT_ROOT <- "FateScape/results/Sensitivity_analysis_beta"

# Optional: if FateScape is not installed as an R package, set package root here.
# Example: FATESCAPE_LOCAL_PATH <- "D:/projects/lineage tracing/FateScape"

CASE_INDICES <- "all"


# Batch behavior passed into the patched base runner.
SKIP_FINISHED <- TRUE
FORCE_RERUN <- FALSE
STOP_ON_ERROR <- FALSE

# Keep the same core FateScape settings except beta proposal weights.
MAX_ITER <- 100
REPEAT_TIME <- 10
LABEL_SITE_FRACTION <- 0.70
LAMBDA_STATE <- 0.10
LAMBDA_BARCODE <- 0.90
DROP_DUP_ALPHA <- 1.5
DROP_DUP_BETA  <- 1.5

# Beta proposal-weight grid.
# beta1_state = proposal_state_weight
# beta2_barcode = proposal_barcode_weight
#
# The original nodes_sampling() defaults are approximately state_weight = 5 and
# barcode_weight = 5. The grid below tests barcode-only, state-only, balanced,
# and moderately dominant settings while keeping the scale close to the default.
BETA_GRID <- data.frame(
  beta_setting = c(
    "beta1_0_beta2_5_barcode_only",
    "beta1_5_beta2_0_state_only",
    "beta1_5_beta2_5_default",
    "beta1_5_beta2_10_barcode_dominant",
    "beta1_10_beta2_5_state_dominant",
    "beta1_2p5_beta2_7p5_barcode_moderate",
    "beta1_7p5_beta2_2p5_state_moderate"
  ),
  beta1_state = c(0.0, 5.0, 5.0, 5.0, 10.0, 2.5, 7.5),
  beta2_barcode = c(5.0, 0.0, 5.0, 10.0, 5.0, 7.5, 2.5),
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
    message("Package ggplot2 is not available; skipping beta-sensitivity plots.")
    return(invisible(NULL))
  }

  metric_cols <- intersect(c("Nye", "TreeDistance", "RF_normalized", "RF_raw"), names(metrics_all))
  if (length(metric_cols) == 0) {
    message("No metric columns available for plotting.")
    return(invisible(NULL))
  }

  metrics_long <- data.frame()
  for (m in metric_cols) {
    tmp <- metrics_all[, intersect(c("case_id", "beta_setting", "beta1_state", "beta2_barcode"), names(metrics_all)), drop = FALSE]
    tmp$metric <- m
    tmp$value <- metrics_all[[m]]
    metrics_long <- rbind(metrics_long, tmp)
  }

  suppressPackageStartupMessages(library(ggplot2))

  p <- ggplot(metrics_long, aes(x = beta_setting, y = value)) +
    geom_boxplot(outlier.size = 0.6) +
    geom_jitter(width = 0.15, size = 0.7, alpha = 0.5) +
    facet_wrap(~ metric, scales = "free_y") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
    labs(
      x = "Beta proposal-weight setting",
      y = "Metric value",
      title = "FateScape sensitivity to beta1/state and beta2/barcode proposal weights"
    )

  ggplot2::ggsave(file.path(out_dir, "beta_sensitivity_metrics_boxplot.pdf"), p, width = 11, height = 5)
  ggplot2::ggsave(file.path(out_dir, "beta_sensitivity_metrics_boxplot.png"), p, width = 11, height = 5, dpi = 300)
}

# -------------------------
# 2. INJECTED BETA FUNCTIONS
# -------------------------

beta_function_block <- function(beta1_state, beta2_barcode, beta_setting) {
  paste0(
'
# ---- Injected by beta-sensitivity wrapper ----
BETA_SENSITIVITY_MODE <- TRUE
BETA_SETTING <- ', quote_r_string(beta_setting), '
BETA1_STATE <- ', sprintf("%.12g", beta1_state), '
BETA2_BARCODE <- ', sprintf("%.12g", beta2_barcode), '

barcode_diff_count_beta_local <- function(x, y) {
  x <- as.vector(x)
  y <- as.vector(y)
  if (length(x) != length(y)) {
    stop("barcode_diff_count_beta_local(): x and y must have the same length.")
  }
  missing_tokens <- c("-", "?", "", "NA", "NaN", "nan", "NULL", "None", "-1")
  x_chr <- as.character(x)
  y_chr <- as.character(y)
  x_missing <- is.na(x) | x_chr %in% missing_tokens
  y_missing <- is.na(y) | y_chr %in% missing_tokens
  one_missing <- xor(x_missing, y_missing)
  both_obs <- !x_missing & !y_missing
  sum(one_missing | (both_obs & x_chr != y_chr), na.rm = TRUE)
}

nodes_sampling_beta_exposed <- function(tree, edges, cell_state_labels, barcodes, state_lineages,
                                        state_threshold = 4,
                                        barcode_threshold = 15,
                                        state_weight = BETA1_STATE,
                                        barcode_weight = BETA2_BARCODE,
                                        undefined_state_penalty = 10) {
  total_nodes <- 2 * length(tree$tip.label) - 1
  renew_prob <- rep(0, total_nodes)
  parent_nodes <- unique(edges[, 1])
  parent_nodes <- parent_nodes[!is.na(parent_nodes)]

  for (par in parent_nodes) {
    children_nodes <- edges[edges[, 1] == par, 2]
    children_nodes <- children_nodes[!is.na(children_nodes)]
    if (length(children_nodes) < 2) next
    children_nodes <- children_nodes[seq_len(2)]

    if (children_nodes[1] > length(tree$tip.label)) {
      child_id_1 <- children_nodes[1]
    } else {
      cell_id_1 <- tree$tip.label[children_nodes[1]]
      child_id_1 <- match(cell_id_1, rownames(barcodes))
    }

    if (children_nodes[2] > length(tree$tip.label)) {
      child_id_2 <- children_nodes[2]
    } else {
      cell_id_2 <- tree$tip.label[children_nodes[2]]
      child_id_2 <- match(cell_id_2, rownames(barcodes))
    }

    if (is.na(child_id_1) || is.na(child_id_2)) next
    if (child_id_1 < 1 || child_id_1 > nrow(barcodes)) next
    if (child_id_2 < 1 || child_id_2 > nrow(barcodes)) next

    child_barcodes_1 <- barcodes[child_id_1, , drop = TRUE]
    child_barcodes_2 <- barcodes[child_id_2, , drop = TRUE]

    parent_state  <- cell_state_labels[par]
    child_state_1 <- cell_state_labels[child_id_1]
    child_state_2 <- cell_state_labels[child_id_2]

    if (is.na(parent_state) || parent_state == 0) {
      renew_prob[children_nodes[1]] <- renew_prob[children_nodes[1]] + undefined_state_penalty
      renew_prob[children_nodes[2]] <- renew_prob[children_nodes[2]] + undefined_state_penalty
    } else {
      for (lineage in state_lineages) {
        if (!is.na(child_state_1) && (parent_state %in% lineage) && (child_state_1 %in% lineage)) {
          state_dist_1 <- match(child_state_1, lineage) - match(parent_state, lineage)
          if (!is.na(state_dist_1) && state_dist_1 > 0 && state_dist_1 < state_threshold) {
            renew_prob[children_nodes[1]] <- renew_prob[children_nodes[1]] - (state_weight / state_dist_1)
          } else {
            renew_prob[children_nodes[1]] <- renew_prob[children_nodes[1]] + state_weight
          }
        }
      }
      for (lineage in state_lineages) {
        if (!is.na(child_state_2) && (parent_state %in% lineage) && (child_state_2 %in% lineage)) {
          state_dist_2 <- match(child_state_2, lineage) - match(parent_state, lineage)
          if (!is.na(state_dist_2) && state_dist_2 > 0 && state_dist_2 < state_threshold) {
            renew_prob[children_nodes[2]] <- renew_prob[children_nodes[2]] - (state_weight / state_dist_2)
          } else {
            renew_prob[children_nodes[2]] <- renew_prob[children_nodes[2]] + state_weight
          }
        }
      }
    }

    diff_sites <- barcode_diff_count_beta_local(child_barcodes_1, child_barcodes_2)
    if (!is.na(diff_sites) && diff_sites < barcode_threshold) {
      renew_prob[children_nodes[1]] <- renew_prob[children_nodes[1]] - (barcode_weight / (diff_sites + 1))
      renew_prob[children_nodes[2]] <- renew_prob[children_nodes[2]] - (barcode_weight / (diff_sites + 1))
    } else {
      renew_prob[children_nodes[1]] <- renew_prob[children_nodes[1]] + barcode_weight
      renew_prob[children_nodes[2]] <- renew_prob[children_nodes[2]] + barcode_weight
    }
  }

  renew_prob[is.na(renew_prob)] <- 0
  renew_prob[renew_prob < 0] <- 0
  total <- sum(renew_prob)

  if (total <= 0 || is.na(total)) {
    renew_prob <- rep(0, total_nodes)
    root_node <- setdiff(edges[, 1], edges[, 2])
    root_node <- root_node[!is.na(root_node)]
    if (length(root_node) == 0) root_node <- unique(edges[, 1])[1]
    non_root <- setdiff(seq_len(total_nodes), root_node[1])
    if (length(non_root) > 0) {
      renew_prob[non_root] <- 1 / length(non_root)
    } else {
      renew_prob[] <- 1 / length(renew_prob)
    }
  } else {
    renew_prob <- renew_prob / total
  }

  renew_prob
}

subtree_refinement_beta_exposed <- function(Trees_initial, state_lineages, barcodes_lineages, N_char,
                                            state_labels_lineages, lambda1, lambda2,
                                            maxIter = 100, repeat_time = 10) {
  ptm_FateScape <- proc.time()
  bestsubtreescore <- list()
  bestsubtree <- list()

  for (j in seq_along(state_lineages)) {
    lineage_label <- paste0("L", j)
    tree <- Trees_initial[[lineage_label]]
    if (is.null(tree)) next
    tree$edge.length <- rep(1, nrow(tree$edge))

    score <- FateScape:::combined_likelihood(
      tree = tree,
      barcodes = barcodes_lineages[[lineage_label]],
      N_char = N_char,
      cell_state_labels = state_labels_lineages[[lineage_label]],
      state_lineages = state_lineages,
      state_score = NULL,
      barcode_score = NULL,
      lambda_1 = lambda1,
      lambda_2 = lambda2
    )
    maxscore <- score

    ances_res <- ancestor_inference(
      tree = tree,
      N_char = N_char,
      barcodes = barcodes_lineages[[lineage_label]],
      cell_state_labels = state_labels_lineages[[lineage_label]],
      state_lineages = state_lineages
    )
    barcodes_ancestral <- ances_res[[1]]
    cell_state_labels_ancestral <- ances_res[[2]]

    score_all <- c()
    iter_repeat <- repeat_time
    best_tree <- tree
    best_tree_list <- list()
    best_score_list <- c()

    for (i in seq_len(maxIter)) {
      renew_prob <- nodes_sampling_beta_exposed(
        tree = tree,
        edges = tree$edge,
        cell_state_labels = cell_state_labels_ancestral,
        barcodes = barcodes_ancestral,
        state_lineages = state_lineages,
        state_weight = BETA1_STATE,
        barcode_weight = BETA2_BARCODE
      )

      tree_new <- FateScape:::subtrees_swapping(
        tree = tree,
        edges = tree$edge,
        cell_state_labels = cell_state_labels_ancestral,
        barcodes = barcodes_ancestral,
        state_lineages = state_lineages,
        tree_renew_prob = renew_prob
      )

      score <- FateScape:::combined_likelihood(
        tree_new,
        barcodes_lineages[[lineage_label]],
        N_char = N_char,
        cell_state_labels = state_labels_lineages[[lineage_label]],
        state_lineages = state_lineages,
        state_score = NULL,
        barcode_score = NULL,
        lambda_1 = lambda1,
        lambda_2 = lambda2
      )

      ances_res <- ancestor_inference(
        tree_new,
        N_char = N_char,
        barcodes = barcodes_lineages[[lineage_label]],
        cell_state_labels = state_labels_lineages[[lineage_label]],
        state_lineages = state_lineages
      )
      barcodes_new <- ances_res[[1]]
      cell_state_labels_new <- ances_res[[2]]

      score_all <- c(score_all, maxscore)

      if (!is.na(score) && score > maxscore) {
        maxscore <- score
        best_tree <- tree_new
        tree <- tree_new
        barcodes_ancestral <- barcodes_new
        cell_state_labels_ancestral <- cell_state_labels_new
      }

      if (i %% iter_repeat == 0) {
        recent_scores <- score_all[(i - iter_repeat + 1):i]
        if (length(unique(recent_scores)) == 1) {
          best_tree_list[[length(best_tree_list) + 1]] <- best_tree
          best_score_list <- c(best_score_list, maxscore)

          tree <- Trees_initial[[lineage_label]]
          tree$edge.length <- rep(1, nrow(tree$edge))
          score <- FateScape:::combined_likelihood(
            tree,
            barcodes_lineages[[lineage_label]],
            N_char,
            state_labels_lineages[[lineage_label]],
            state_lineages,
            state_score = NULL,
            barcode_score = NULL,
            lambda_1 = lambda1,
            lambda_2 = lambda2
          )
          maxscore <- score
          ances_res <- ancestor_inference(
            tree,
            N_char,
            barcodes_lineages[[lineage_label]],
            state_labels_lineages[[lineage_label]],
            state_lineages
          )
          barcodes_ancestral <- ances_res[[1]]
          cell_state_labels_ancestral <- ances_res[[2]]
        }
      }
    }

    if (length(best_score_list) == 0) {
      best_score_list <- maxscore
      best_tree_list <- list(best_tree)
    }

    bestsubtreescore[[lineage_label]] <- max(best_score_list)
    best_index <- match(max(best_score_list), best_score_list)
    bestsubtree[[lineage_label]] <- best_tree_list[[best_index]]
  }

  total_time <- proc.time() - ptm_FateScape
  list(bestsubtree = bestsubtree, bestsubtreescore = bestsubtreescore, total_time = total_time)
}
'
  )
}

insert_beta_functions <- function(txt, beta1_state, beta2_barcode, beta_setting) {
  marker <- "set.seed(GLOBAL_SEED)"
  if (!grepl(marker, txt, fixed = TRUE)) {
    stop("Cannot find insertion marker in BASE_RUNNER: ", marker)
  }
  block <- beta_function_block(beta1_state, beta2_barcode, beta_setting)
  sub(marker, paste0(marker, "\n\n", block), txt, fixed = TRUE)
}

patch_subtree_refinement_call <- function(txt) {
  # Replace only the main tree-refinement call in the base runner.
  # This keeps all downstream FateScape steps unchanged.
  txt2 <- sub(
    "LiterateRefine\\s*<-\\s*subtree_refinement\\s*\\(",
    "LiterateRefine <- subtree_refinement_beta_exposed(",
    txt,
    perl = TRUE
  )
  if (identical(txt, txt2)) {
    stop("Could not patch the subtree_refinement() call in BASE_RUNNER.")
  }
  txt2
}

write_patched_runner <- function(base_runner, out_file,
                                 method_input_root,
                                 method_output_root,
                                 case_indices,
                                 beta_setting,
                                 beta1_state,
                                 beta2_barcode) {
  txt <- readLines(base_runner, warn = FALSE, encoding = "UTF-8")
  txt <- paste(txt, collapse = "\n")

  # Paths and batch selection.
  txt <- replace_assignment(txt, "METHOD_INPUT_ROOT", quote_r_string(method_input_root))
  txt <- replace_assignment(txt, "METHOD_OUTPUT_ROOT", quote_r_string(method_output_root))
  txt <- replace_assignment(txt, "CASE_INDICES", quote_r_string(as.character(case_indices)))
  txt <- replace_assignment(txt, "SKIP_FINISHED", ifelse(SKIP_FINISHED, "TRUE", "FALSE"))
  txt <- replace_assignment(txt, "FORCE_RERUN", ifelse(FORCE_RERUN, "TRUE", "FALSE"))
  txt <- replace_assignment(txt, "STOP_ON_ERROR", ifelse(STOP_ON_ERROR, "TRUE", "FALSE"))

  # Keep non-beta settings fixed.
  txt <- replace_assignment(txt, "LABEL_SITE_FRACTION", as.character(LABEL_SITE_FRACTION), warn_if_missing = FALSE)
  txt <- replace_assignment(txt, "MAX_ITER", as.character(MAX_ITER), warn_if_missing = FALSE)
  txt <- replace_assignment(txt, "REPEAT_TIME", as.character(REPEAT_TIME), warn_if_missing = FALSE)
  txt <- replace_assignment(txt, "LAMBDA_STATE", as.character(LAMBDA_STATE), warn_if_missing = FALSE)
  txt <- replace_assignment(txt, "LAMBDA_BARCODE", as.character(LAMBDA_BARCODE), warn_if_missing = FALSE)
  txt <- replace_assignment(txt, "DROP_DUP_ALPHA", as.character(DROP_DUP_ALPHA), warn_if_missing = FALSE)
  txt <- replace_assignment(txt, "DROP_DUP_BETA", as.character(DROP_DUP_BETA), warn_if_missing = FALSE)

  # Inject beta-specific local refinement and patch the refinement call.
  txt <- insert_beta_functions(txt, beta1_state, beta2_barcode, beta_setting)
  txt <- patch_subtree_refinement_call(txt)

  writeLines(txt, out_file, useBytes = TRUE)
  out_file
}

# -------------------------
# 3. VALIDATE SETTINGS
# -------------------------

if (!file.exists(BASE_RUNNER)) stop("BASE_RUNNER not found: ", BASE_RUNNER)
if (!dir.exists(METHOD_INPUT_ROOT)) stop("METHOD_INPUT_ROOT not found: ", METHOD_INPUT_ROOT)

SENS_OUTPUT_ROOT <- safe_dir_create(SENS_OUTPUT_ROOT)
RUNNER_DIR <- safe_dir_create(file.path(SENS_OUTPUT_ROOT, "_patched_runners"))
SUMMARY_DIR <- safe_dir_create(file.path(SENS_OUTPUT_ROOT, "_batch_summary"))

safe_write_csv(BETA_GRID, file.path(SENS_OUTPUT_ROOT, "beta_sensitivity_grid.csv"), row.names = FALSE)

message("Base runner       : ", BASE_RUNNER)
message("Method input root : ", METHOD_INPUT_ROOT)
message("Output root       : ", SENS_OUTPUT_ROOT)
message("Case indices      : ", CASE_INDICES)
message("Number of beta settings: ", nrow(BETA_GRID))

# -------------------------
# 4. RUN SENSITIVITY GRID
# -------------------------

all_metric_records <- list()
all_manifest_records <- list()
all_error_records <- list()
run_plan_records <- list()

for (ii in seq_len(nrow(BETA_GRID))) {
  beta_setting <- BETA_GRID$beta_setting[ii]
  beta1_state <- BETA_GRID$beta1_state[ii]
  beta2_barcode <- BETA_GRID$beta2_barcode[ii]

  out_root <- file.path(SENS_OUTPUT_ROOT, beta_setting)
  safe_dir_create(out_root)

  message("\n============================================================")
  message("[FateScape beta sensitivity] ", beta_setting)
  message("beta1/state   = ", beta1_state)
  message("beta2/barcode = ", beta2_barcode)
  message("Output        = ", out_root)

  patched_runner <- file.path(RUNNER_DIR, paste0("run_fatescape_", beta_setting, ".R"))

  write_patched_runner(
    base_runner = BASE_RUNNER,
    out_file = patched_runner,
    method_input_root = METHOD_INPUT_ROOT,
    method_output_root = out_root,
    case_indices = CASE_INDICES,
    beta_setting = beta_setting,
    beta1_state = beta1_state,
    beta2_barcode = beta2_barcode
  )

  run_plan_records[[length(run_plan_records) + 1L]] <- data.frame(
    beta_setting = beta_setting,
    beta1_state = beta1_state,
    beta2_barcode = beta2_barcode,
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
      message("[ERROR] Beta setting failed: ", beta_setting, " | ", msg)
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
    metrics$beta_setting <- beta_setting
    metrics$beta1_state <- beta1_state
    metrics$beta2_barcode <- beta2_barcode
    metrics$beta_runner <- patched_runner
    all_metric_records[[length(all_metric_records) + 1L]] <- metrics
  }

  if (!is.null(manifest)) {
    manifest$beta_setting <- beta_setting
    manifest$beta1_state <- beta1_state
    manifest$beta2_barcode <- beta2_barcode
    manifest$beta_runner <- patched_runner
    all_manifest_records[[length(all_manifest_records) + 1L]] <- manifest
  }

  if (!is.null(errors)) {
    errors$beta_setting <- beta_setting
    errors$beta1_state <- beta1_state
    errors$beta2_barcode <- beta2_barcode
    errors$beta_runner <- patched_runner
    all_error_records[[length(all_error_records) + 1L]] <- errors
  }

  if (!identical(result, "success")) {
    all_error_records[[length(all_error_records) + 1L]] <- data.frame(
      beta_setting = beta_setting,
      beta1_state = beta1_state,
      beta2_barcode = beta2_barcode,
      status = "runner_failed",
      error = result,
      beta_runner = patched_runner,
      stringsAsFactors = FALSE
    )
  }
}

# -------------------------
# 5. WRITE COMBINED SUMMARIES
# -------------------------

run_plan <- if (length(run_plan_records) > 0) do.call(rbind, run_plan_records) else data.frame()
safe_write_csv(run_plan, file.path(SUMMARY_DIR, "beta_sensitivity_run_plan.csv"), row.names = FALSE)

if (length(all_metric_records) > 0) {
  metrics_all <- do.call(rbind, all_metric_records)
  safe_write_csv(metrics_all, file.path(SENS_OUTPUT_ROOT, "beta_sensitivity_metrics_all.csv"), row.names = FALSE)
  safe_write_csv(metrics_all, file.path(SUMMARY_DIR, "beta_sensitivity_metrics_all.csv"), row.names = FALSE)
  plot_sensitivity(metrics_all, SENS_OUTPUT_ROOT)
  message("Combined metrics written to: ", file.path(SENS_OUTPUT_ROOT, "beta_sensitivity_metrics_all.csv"))
} else {
  message("No metrics were collected.")
}

if (length(all_manifest_records) > 0) {
  manifest_all <- do.call(rbind, all_manifest_records)
  safe_write_csv(manifest_all, file.path(SENS_OUTPUT_ROOT, "beta_sensitivity_manifest_all.csv"), row.names = FALSE)
  safe_write_csv(manifest_all, file.path(SUMMARY_DIR, "beta_sensitivity_manifest_all.csv"), row.names = FALSE)
  message("Combined manifest written to: ", file.path(SENS_OUTPUT_ROOT, "beta_sensitivity_manifest_all.csv"))
  print(table(manifest_all$beta_setting, manifest_all$status, useNA = "ifany"))
} else {
  message("No manifests were collected.")
}

if (length(all_error_records) > 0) {
  errors_all <- do.call(rbind, all_error_records)
  safe_write_csv(errors_all, file.path(SENS_OUTPUT_ROOT, "beta_sensitivity_errors_all.csv"), row.names = FALSE)
  safe_write_csv(errors_all, file.path(SUMMARY_DIR, "beta_sensitivity_errors_all.csv"), row.names = FALSE)
  message("Combined errors written to: ", file.path(SENS_OUTPUT_ROOT, "beta_sensitivity_errors_all.csv"))
}

message("\nFateScape beta-sensitivity completed.")
message("Output root: ", SENS_OUTPUT_ROOT)
