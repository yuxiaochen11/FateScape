# ============================================================
# Generate 100-cell small-scale benchmark datasets
# for TiDeTree-scale supplemental comparison
#
# IMPORTANT:
#   Run this script AFTER running the main data-generation Rmd
#   so that generate_one_case(), OUT_ROOT, safe_write_csv(), and
#   other helper functions are already available in the R session.
#
#   This v2 fixes the previous failure "更换参数长度为零" by using
#   the exact parameter names expected by generate_one_case():
#     n_cif, cif_step, N_ms, r_n
#   rather than n_CIF / missing cif_step.
# ============================================================

# -------------------------
# 0. Output location
# -------------------------
OUT_ROOT <- "D:/projects/lineage tracing/FateScape/Sensitivity_data/input"

# -------------------------
# 1. Required function checks
# -------------------------

generate_one_case <- function(params, out_root = OUT_ROOT) {
  params <- as.list(params)
  case_dir <- file.path(out_root, params$pool, params$case_id)
  case_manifest_file <- file.path(case_dir, "manifest.csv")

  if (dir.exists(case_dir)) {
    if (isTRUE(SKIP_EXISTING_CASES)) {
      message("\n[Skip existing] ", params$case_id)
      message("Existing output: ", normalizePath(case_dir, winslash = "/", mustWork = FALSE))

      if (file.exists(case_manifest_file)) {
        manifest <- read.csv(case_manifest_file, stringsAsFactors = FALSE)
        manifest$generation_status <- "skipped_existing"
        return(manifest)
      }

      return(data.frame(
        case_id = params$case_id,
        case_index = params$case_index,
        pool = params$pool,
        figure = params$figure,
        panel = params$panel,
        x_variable = params$x_variable,
        ncells = params$ncells,
        N_char = params$N_char,
        mu = params$mu,
        p_ds = params$p_ds,
        p_d = params$p_d,
        replicate = params$replicate,
        seed = params$seed,
        case_dir = normalizePath(case_dir, winslash = "/", mustWork = FALSE),
        generation_status = "skipped_existing_no_manifest",
        warning = "Case directory exists but manifest.csv was not found; no files were overwritten.",
        stringsAsFactors = FALSE
      ))
    }

    if (isTRUE(OVERWRITE_CASES)) {
      message("\n[Overwrite existing] ", params$case_id)
      unlink(case_dir, recursive = TRUE, force = TRUE)
    } else {
      stop(
        "Case directory already exists: ", case_dir,
        ". Set SKIP_EXISTING_CASES <- TRUE to skip it, or OVERWRITE_CASES <- TRUE to regenerate it."
      )
    }
  }

  dir.create(case_dir, recursive = TRUE, showWarnings = FALSE)

  message("\n[Generating] ", params$case_id)
  message("Output: ", normalizePath(case_dir, winslash = "/", mustWork = FALSE))

  # Fixed cell-state tree from the original Rmd.
  cell_state_tree_text <- "((((((t1:1)n9:1,(t2:1)n10:1)n8:1)n7:1,(((t3:1)n13:1)n12:1)n11:1)n6:1)n5:1)n4;"
  cell_state_tree <- ape::read.tree(text = cell_state_tree_text)
  state_lineages <- make_state_lineages(cell_state_tree)
  s_depth <- state_depths(cell_state_tree)
  terminal_states <- seq_along(cell_state_tree$tip.label)

  # Generate raw lineage barcodes, cell states, CIFs, and true tree using original TedSim logic.
  set.seed(params$seed)
  returnlist <- SIFGenerate(cell_state_tree, params$n_diff, step = params$cif_step)

  cifs <- simulate_cifs_retry(
    ncells = params$ncells,
    cell_state_tree = cell_state_tree,
    p_a = params$p_a,
    n_CIF = params$n_cif,
    n_diff = params$n_diff,
    step = params$cif_step,
    p_d = params$p_d,
    Sigma = params$Sigma,
    mu = params$mu,
    N_char = params$N_char,
    max_walk_candidates = unique(c(params$max_walk, MAX_WALK_CANDIDATES)),
    SIF_res = returnlist,
    unif_on = FALSE,
    n_retries = SIMULATE_CIFS_RETRIES,
    seed = params$seed
  )

  params$simulate_attempt <- attr(cifs, "simulate_attempt")
  params$max_walk_used <- attr(cifs, "max_walk_used")
  params$seed_used <- attr(cifs, "seed_used")
  params$simulator_backend <- "TedSim_SimulateCIFs"

  barcodes_all <- cifs[[7]]

  # CIF list for leaves. Original code assumes three terminal cell states.
  n_terminal <- length(cell_state_tree$tip.label)
  cif_leaves <- vector("list", n_terminal)
  for (i in seq_len(n_terminal)) {
    cif_leaves[[i]] <- cifs[[1]][[i]][seq_len(params$ncells), , drop = FALSE]
  }

  true_counts_res <- CIF2Truecounts(
    ngenes = params$ngenes,
    ncif = params$n_cif,
    ge_prob = params$ge_prob,
    ncells = params$ncells,
    cif_res = list(cif_leaves, cifs[[2]])
  )

  counts_true <- t(true_counts_res[[1]])
  states_leaves <- cifs[[2]][seq_len(params$ncells), , drop = FALSE]

  cell_ids <- paste0("cell_", states_leaves[, "cellID"])
  rownames(counts_true) <- cell_ids
  colnames(counts_true) <- paste0("Gene", seq_len(ncol(counts_true)))

  counts_expr <- counts_true
  expression_backend <- "TedSim_CIF2Truecounts"

  if (isTRUE(SIMULATE_OBSERVED_COUNTS) && exists("True2ObservedCounts", mode = "function", inherits = TRUE)) {
    obs_try <- tryCatch(
      True2ObservedCounts(counts_true),
      error = function(e) {
        warning("True2ObservedCounts failed: ", conditionMessage(e), ". Using true counts as expression.")
        NULL
      }
    )
    if (!is.null(obs_try)) {
      counts_expr <- obs_try
      expression_backend <- "TedSim_True2ObservedCounts"
    }
  }
  params$expression_backend <- expression_backend

  # Ground-truth tree.
  tree_ct <- cifs[[4]]
  tree_ct <- rename_tree_tips_to_cells(tree_ct, params$ncells)

  # Barcode matrices.
  barcodes_initial <- barcodes_all[seq_len(params$ncells), , drop = FALSE]
  barcode_true <- standardize_barcode_matrix(barcodes_initial, cell_ids = cell_ids, prefix = "S")

  qc <- barcode_qc(barcode_true, params)
  params <- c(params, qc)

  dropout_res <- add_controlled_stochastic_dropout(
    barcode_true,
    p_ds = params$p_ds,
    seed = params$seed + 777L
  )
  barcode_observed <- dropout_res$barcode_observed
  barcode_dropout_type <- dropout_res$dropout_type


  # barcode_imputed <- try_fatescape_imputation(
  #   barcode_observed = barcode_observed,
  #   N_char = params$N_char,
  #   ncells = params$ncells,
  #   r_n = params$r_n
  # )


  # Metadata.
  state_label <- as.integer(states_leaves[, "cluster"])
  state_type <- ifelse(state_label %in% terminal_states, "terminal", "internal_or_root")
  state_depth <- as.integer(s_depth[as.character(state_label)])

  cell_metadata <- data.frame(
    cell_id = cell_ids,
    state_label = state_label,
    state_type = state_type,
    state_depth = state_depth,
    ted_parent_state = as.integer(states_leaves[, "parent"]),
    ted_cluster = as.integer(states_leaves[, "cluster"]),
    ted_depth = as.integer(states_leaves[, "depth"]),
    ted_cellID = as.integer(states_leaves[, "cellID"]),
    stringsAsFactors = FALSE
  )

  represented_states <- sort(unique(state_label))
  params$n_sampled_states <- length(represented_states)
  params$n_sampled_terminal_states <- sum(represented_states %in% terminal_states)
  params$n_sampled_internal_states <- sum(!(represented_states %in% terminal_states))
  params$terminal_state_cell_fraction <- mean(state_type == "terminal")
  params$internal_state_cell_fraction <- mean(state_type == "internal_or_root")
  params$sampled_states <- paste(represented_states, collapse = ";")

  # Canonical standardized outputs only.
  # Method-specific inputs should be produced later by converter scripts from these files.
  expression_file <- file.path(case_dir, "expression.csv")
  true_counts_file <- file.path(case_dir, "true_counts.csv")
  barcode_true_file <- file.path(case_dir, "barcode_true.csv")
  barcode_observed_file <- file.path(case_dir, "barcode_observed.csv")
  # barcode_imputed_file <- file.path(case_dir, "barcode_imputed_fatescape.csv")
  barcode_dropout_type_file <- file.path(case_dir, "barcode_dropout_type.csv")
  metadata_file <- file.path(case_dir, "cell_metadata.csv")
  true_tree_file <- file.path(case_dir, "true_tree.nwk")
  state_tree_file <- file.path(case_dir, "state_tree.nwk")
  state_lineages_file <- file.path(case_dir, "state_lineages.rds")

  safe_write_csv(as.data.frame(counts_expr, check.names = FALSE), expression_file, row.names = TRUE)
  safe_write_csv(as.data.frame(counts_true, check.names = FALSE), true_counts_file, row.names = TRUE)
  safe_write_csv(as.data.frame(barcode_true, check.names = FALSE), barcode_true_file, row.names = TRUE)
  safe_write_csv(as.data.frame(barcode_observed, check.names = FALSE), barcode_observed_file, row.names = TRUE)
  # safe_write_csv(as.data.frame(barcode_imputed, check.names = FALSE), barcode_imputed_file, row.names = TRUE)
  safe_write_csv(as.data.frame(barcode_dropout_type, check.names = FALSE), barcode_dropout_type_file, row.names = TRUE)
  safe_write_csv(cell_metadata, metadata_file, row.names = FALSE)
  safe_write_tree(tree_ct, true_tree_file)
  safe_write_tree(cell_state_tree, state_tree_file)
  safe_save_rds(state_lineages, state_lineages_file)

  # RDS with original TedSim objects.
  safe_save_rds(
    list(
      cifs = cifs,
      returnlist = returnlist,
      true_counts_res = true_counts_res,
      state_lineages = state_lineages,
      params = params
    ),
    file.path(case_dir, "tedsim_objects.rds")
  )

  params_df <- as.data.frame(params, stringsAsFactors = FALSE)
  safe_write_csv(params_df, file.path(case_dir, "params.csv"), row.names = FALSE)

  manifest <- data.frame(
    case_id = params$case_id,
    case_index = params$case_index,
    pool = params$pool,
    figure = params$figure,
    panel = params$panel,
    x_variable = params$x_variable,
    ncells = params$ncells,
    N_char = params$N_char,
    mu = params$mu,
    p_ds = params$p_ds,
    p_d = params$p_d,
    replicate = params$replicate,
    seed = params$seed,
    seed_used = params$seed_used,
    max_walk_used = params$max_walk_used,
    simulate_attempt = params$simulate_attempt,
    simulator_backend = params$simulator_backend,
    expression_backend = params$expression_backend,
    case_dir = normalizePath(case_dir, winslash = "/", mustWork = TRUE),
    expression_file = normalizePath(expression_file, winslash = "/", mustWork = TRUE),
    true_counts_file = normalizePath(true_counts_file, winslash = "/", mustWork = TRUE),
    barcode_true_file = normalizePath(barcode_true_file, winslash = "/", mustWork = TRUE),
    barcode_observed_file = normalizePath(barcode_observed_file, winslash = "/", mustWork = TRUE),
    # barcode_imputed_file = normalizePath(barcode_imputed_file, winslash = "/", mustWork = TRUE),
    barcode_dropout_type_file = normalizePath(barcode_dropout_type_file, winslash = "/", mustWork = TRUE),
    metadata_file = normalizePath(metadata_file, winslash = "/", mustWork = TRUE),
    true_tree_file = normalizePath(true_tree_file, winslash = "/", mustWork = TRUE),
    state_tree_file = normalizePath(state_tree_file, winslash = "/", mustWork = TRUE),
    state_lineages_file = normalizePath(state_lineages_file, winslash = "/", mustWork = TRUE),
    params_file = normalizePath(file.path(case_dir, "params.csv"), winslash = "/", mustWork = TRUE),
    tedsim_objects_file = normalizePath(file.path(case_dir, "tedsim_objects.rds"), winslash = "/", mustWork = TRUE),
    n_mutated_entries_true = params$n_mutated_entries_true,
    n_missing_entries_true = params$n_missing_entries_true,
    mutation_fraction_true = params$mutation_fraction_true,
    missing_fraction_true = params$missing_fraction_true,
    n_sampled_states = params$n_sampled_states,
    n_sampled_terminal_states = params$n_sampled_terminal_states,
    n_sampled_internal_states = params$n_sampled_internal_states,
    terminal_state_cell_fraction = params$terminal_state_cell_fraction,
    internal_state_cell_fraction = params$internal_state_cell_fraction,
    sampled_states = params$sampled_states,
    generation_status = "generated",
    stringsAsFactors = FALSE
  )

  safe_write_csv(manifest, file.path(case_dir, "manifest.csv"), row.names = FALSE)

  message("[Done] ", params$case_id)
  manifest
}



if (!exists("generate_one_case", mode = "function", inherits = TRUE)) {
  stop(
    "generate_one_case() is not found. Please first run the main ",
    "data_generation_batch_original_results_only.Rmd in the same R session."
  )
}

if (!exists("make_case_id", mode = "function", inherits = TRUE)) {
  # Same convention as the main Rmd.
  p_tag <- function(x) gsub("\\\\.", "p", as.character(x))
  make_case_id <- function(pool, ncells, N_char, mu, p_ds, p_d, replicate) {
    paste0(
      pool,
      "_n", sprintf("%04d", ncells),
      "_char", N_char,
      "_mu", p_tag(mu),
      "_sd", p_tag(p_ds),
      "_hd", p_tag(p_d),
      "_rep", sprintf("%03d", replicate)
    )
  }
}



# Local grid builder. Does not overwrite the main Rmd's base_grid().
make_small_grid <- function(pool, figure, panel, x_variable,
                            ncells, N_char, mu, p_ds, p_d, replicates) {
  expand.grid(
    pool = pool,
    figure = figure,
    panel = panel,
    x_variable = x_variable,
    ncells = ncells,
    N_char = N_char,
    mu = mu,
    p_ds = p_ds,
    p_d = p_d,
    replicate = seq_len(replicates),
    stringsAsFactors = FALSE
  )
}

# -------------------------
# 2. Small-scale parameter settings
# -------------------------
N_REPLICATES_SMALL <- 10
INCLUDE_STRESS <- TRUE
BASE_SEED_SMALL <- 319

small_grid <- make_small_grid(
  pool = "1024_default",
  figure = "Supplemental",
  panel = "Sensitivity",
  x_variable = "small_scale_default",
  ncells = 1024,
  N_char = 12,
  mu = 0.2,
  p_ds = 0.2,
  p_d = 0.05,
  replicates = N_REPLICATES_SMALL
)

if (isTRUE(INCLUDE_STRESS)) {
  stress_grid <- make_small_grid(
    pool = "1024_dropout_stress",
    figure = "Supplemental",
    panel = "Sensitivity",
    x_variable = "small_scale_dropout_stress",
    ncells = 1024,
    N_char = 12,
    mu = 0.1,
    p_ds = 0.3,
    p_d = 0.05,
    replicates = N_REPLICATES_SMALL
  )
  small_grid <- rbind(small_grid, stress_grid)
}

# -------------------------
# 3. Add required columns expected by generate_one_case()
# -------------------------
small_grid$ngenes <- 500
small_grid$max_walk <- 3
small_grid$p_a <- 0.8

# IMPORTANT: generate_one_case() expects lowercase n_cif, not n_CIF.
small_grid$n_cif <- 30
small_grid$n_diff <- 20
small_grid$cif_step <- 1
small_grid$ge_prob <- 0.3
small_grid$Sigma <- 0.5
small_grid$N_ms <- 10
small_grid$r_n <- ceiling(small_grid$N_char * 0.7)
small_grid$seed <- BASE_SEED_SMALL + seq_len(nrow(small_grid)) * 1000L

small_grid$case_id <- mapply(
  make_case_id,
  pool = small_grid$pool,
  ncells = small_grid$ncells,
  N_char = small_grid$N_char,
  mu = small_grid$mu,
  p_ds = small_grid$p_ds,
  p_d = small_grid$p_d,
  replicate = small_grid$replicate,
  USE.NAMES = FALSE
)
small_grid <- small_grid[!duplicated(small_grid$case_id), , drop = FALSE]
small_grid$case_index <- seq_len(nrow(small_grid))

# -------------------------
# 4. Write planned cases
# -------------------------
dir.create(file.path(OUT_ROOT, "_manifests"), recursive = TRUE, showWarnings = FALSE)
planned_small_file <- file.path(OUT_ROOT, "_manifests", "planned_cases_1024.csv")
write.csv(small_grid, planned_small_file, row.names = FALSE, quote = FALSE)

message("sensitivity planned cases written to: ", planned_small_file)
message("Number of sensitivity cases: ", nrow(small_grid))
print(small_grid[, c(
  "case_index", "pool", "case_id", "ncells", "N_char",
  "mu", "p_ds", "p_d", "replicate", "seed", "n_cif", "cif_step"
)])

# -------------------------
# 5. Generation settings
# -------------------------
# This script should not overwrite completed cases. It removes only incomplete
# case folders that lack manifest.csv, which commonly remain after failed runs.
SKIP_EXISTING_CASES <- TRUE
OVERWRITE_CASES <- FALSE

# generate_one_case() uses SIMULATE_CIFS_RETRIES and MAX_WALK_CANDIDATES.
# Increase retries for small 100-cell cases because SimulateCIFs can fail for
# specific random seeds.
SIMULATE_CIFS_RETRIES <- 20
MAX_WALK_CANDIDATES <- unique(c(3, 4, 5))

# -------------------------
# 6. Remove incomplete folders from previous failed runs
# -------------------------
case_dirs <- file.path(OUT_ROOT, small_grid$pool, small_grid$case_id)
for (d in case_dirs) {
  manifest_file <- file.path(d, "manifest.csv")
  if (dir.exists(d) && !file.exists(manifest_file)) {
    message("[Remove incomplete case directory] ", d)
    unlink(d, recursive = TRUE, force = TRUE)
  }
}

# -------------------------
# 7. Generate cases
# -------------------------
result_records <- list()
error_records <- list()

for (i in seq_len(nrow(small_grid))) {
  params <- small_grid[i, , drop = FALSE]
  message("\n============================================================")
  message("[sensitivity generation] ", params$case_id)
  message("Output: ", file.path(OUT_ROOT, params$pool, params$case_id))

  res <- tryCatch(
    {
      out <- generate_one_case(params, OUT_ROOT)
      result_records[[length(result_records) + 1L]] <- data.frame(
        case_id = as.character(params$case_id[[1]]),
        pool = as.character(params$pool[[1]]),
        case_dir = normalizePath(
          file.path(OUT_ROOT, as.character(params$pool[[1]]), as.character(params$case_id[[1]])),
          winslash = "/",
          mustWork = FALSE
        ),
        status = "success",
        error = NA_character_,
        stringsAsFactors = FALSE
      )
      out
    },
    error = function(e) {
      msg <- conditionMessage(e)
      message("[ERROR] ", params$case_id, ": ", msg)
      error_records[[length(error_records) + 1L]] <- data.frame(
        case_id = as.character(params$case_id[[1]]),
        pool = as.character(params$pool[[1]]),
        case_dir = normalizePath(
          file.path(OUT_ROOT, as.character(params$pool[[1]]), as.character(params$case_id[[1]])),
          winslash = "/",
          mustWork = FALSE
        ),
        status = "failed",
        error = msg,
        stringsAsFactors = FALSE
      )
      NULL
    }
  )
}

# -------------------------
# 8. Summarize generation status
# -------------------------
status_df <- do.call(rbind, c(result_records, error_records))
if (is.null(status_df)) {
  status_df <- data.frame(case_id = character(), pool = character(), case_dir = character(), status = character(), error = character())
}

status_file <- file.path(OUT_ROOT, "_manifests", "sensitivity_generation_status.csv")
write.csv(status_df, status_file, row.names = FALSE, quote = FALSE)

message("\nsensitivity generation status written to: ", status_file)
print(table(status_df$status))

# -------------------------
# 9. Collect manifest files for small-scale cases
# -------------------------
manifest_files <- file.path(case_dirs, "manifest.csv")
manifest_files <- manifest_files[file.exists(manifest_files)]

if (length(manifest_files) > 0) {
  small_manifest <- do.call(
    rbind,
    lapply(manifest_files, function(f) {
      read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    })
  )

  small_manifest_file <- file.path(OUT_ROOT, "_manifests", "all_cases_manifest_sensitivity.csv")
  write.csv(small_manifest, small_manifest_file, row.names = FALSE, quote = FALSE)

  message("sensitivity all-case manifest written to: ", small_manifest_file)
  message("Number of successfully generated sensitivity cases: ", nrow(small_manifest))
  print(table(small_manifest$pool))

  missing_cases <- setdiff(small_grid$case_id, small_manifest$case_id)
} else {
  warning("No sensitivity case manifest.csv files were generated.")
  small_manifest <- data.frame()
  missing_cases <- small_grid$case_id
}

# -------------------------
# 10. Completeness check
# -------------------------
missing_file <- file.path(OUT_ROOT, "_manifests", "missing_cases_sensitivity.csv")
if (length(missing_cases) > 0) {
  missing_df <- small_grid[small_grid$case_id %in% missing_cases, , drop = FALSE]
  write.csv(missing_df, missing_file, row.names = FALSE, quote = FALSE)
} else {
  write.csv(data.frame(), missing_file, row.names = FALSE, quote = FALSE)
}

message("Planned sensitivity cases: ", nrow(small_grid))
message("Generated sensitivity cases: ", length(manifest_files))
message("Missing sensitivity cases: ", length(missing_cases))
if (length(missing_cases) > 0) {
  message("Missing cases written to: ", missing_file)
  print(missing_cases)
}
