#!/usr/bin/env Rscript


# -------------------------
# 0. USER SETTINGS
# -------------------------
METHOD_INPUT_ROOT <- "FateScape/data/Simulation/method_inputs"
METHOD_OUTPUT_ROOT <- "FateScape/results/Simulation/benchmark_result"

# Optional: if FateScape is not installed as an R package, set package root here.
# Example: FATESCAPE_LOCAL_PATH <- "D:/projects/lineage tracing/FateScape"

CASE_INDICES <- "all"

# Skip/rerun/error behavior.
SKIP_FINISHED <- TRUE
FORCE_RERUN <- FALSE
STOP_ON_ERROR <- FALSE

# Barcode and FateScape reconstruction parameters.
USE_PREIMPUTED_BARCODE <- FALSE
LABEL_SITE_FRACTION <- 0.70
LAMBDA_BARCODE <- 0.90
LAMBDA_STATE <- 0.10
MAX_ITER <- 100
REPEAT_TIME <- 10
DROP_DUP_ALPHA <- 1.5
DROP_DUP_BETA <- 1.5

# Monocle3 settings.
MONOCLE_NUM_DIM <- 20
MONOCLE_REDUCTION <- "UMAP"
MONOCLE_PREPROCESS_METHOD <- "PCA"
MONOCLE_K <- 20
MONOCLE_RESOLUTION <- 1e-3
MONOCLE_MIN_CLUSTERS <- 2
# Root selection for pseudotime: "pc1_min", "umap1_min", or "largest_low_pc1_cluster".
MONOCLE_ROOT_MODE <- "pc1_min"
# If Monocle3 graph / pseudotime fails, fall back to clusters ordered by PC1.
ALLOW_MONOCLE_FALLBACK <- TRUE

# Batch storage options.
SAVE_INTERMEDIATE_OBJECTS <- FALSE
SAVE_FAILED_INTERMEDIATE_OBJECTS <- FALSE

# Reproducibility.
GLOBAL_SEED <- 20260516

# -------------------------
# 1. LOAD PACKAGES
# -------------------------
message("Loading packages...")

if (!is.null(FATESCAPE_LOCAL_PATH) && dir.exists(FATESCAPE_LOCAL_PATH)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::load_all(FATESCAPE_LOCAL_PATH)
} else {
  library(FateScape)
}

suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
  library(TreeTools)
  library(TreeDist)
  library(monocle3)
  library(Matrix)
})

set.seed(GLOBAL_SEED)

# -------------------------
# 2. HELPERS
# -------------------------
safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

safe_write_csv <- function(x, file, row.names = FALSE) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, file = file, row.names = row.names)
}

read_table_with_rownames <- function(file) {
  if (!file.exists(file)) stop("File not found: ", file)
  x <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(x) < 2) stop("Input table has fewer than two columns: ", file)

  first_col <- names(x)[1]
  if (tolower(first_col) %in% c("cell_id", "cell", "cells", "id", "x")) {
    rn <- x[[1]]
    x <- x[, -1, drop = FALSE]
    rownames(x) <- rn
  } else if (anyDuplicated(x[[1]]) == 0 && !is.numeric(x[[1]])) {
    rn <- x[[1]]
    x <- x[, -1, drop = FALSE]
    rownames(x) <- rn
  }
  x
}

standardize_barcode <- function(x) {
  x <- as.data.frame(x, check.names = FALSE)
  x[] <- lapply(x, function(v) {
    v <- as.character(v)
    v[v %in% c("", "NA", "NaN", "nan", "NULL", "null", "None", "none", "?", "-", "-1")] <- NA
    suppressWarnings(as.integer(v))
  })
  as.matrix(x)
}

load_case_inputs <- function(case_dir) {
  files <- list(
    expression = file.path(case_dir, "expression.csv"),
    barcode_observed = file.path(case_dir, "barcode_observed.csv"),
    barcode_imputed = file.path(case_dir, "barcode_imputed.csv"),
    metadata = file.path(case_dir, "metadata.csv"),
    true_tree = file.path(case_dir, "true_tree.nwk")
  )
  if (!file.exists(files$metadata)) {
    alt <- file.path(case_dir, "cell_metadata.csv")
    if (file.exists(alt)) files$metadata <- alt
  }
  if (!file.exists(files$barcode_imputed)) {
    alt <- file.path(case_dir, "barcode_imputed_fatescape.csv")
    if (file.exists(alt)) files$barcode_imputed <- alt
  }

  expression <- read_table_with_rownames(files$expression)
  barcode_observed <- standardize_barcode(read_table_with_rownames(files$barcode_observed))
  barcode_imputed <- NULL
  if (file.exists(files$barcode_imputed)) barcode_imputed <- standardize_barcode(read_table_with_rownames(files$barcode_imputed))

  metadata <- read.csv(files$metadata, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"cell_id" %in% names(metadata)) {
    if ("cell" %in% names(metadata)) metadata$cell_id <- metadata$cell else stop("metadata must contain cell_id")
  }
  true_tree <- ape::read.tree(files$true_tree)

  list(expression = expression, barcode_observed = barcode_observed, barcode_imputed = barcode_imputed,
       metadata = metadata, true_tree = true_tree, files = files)
}

align_expression_barcode_tree <- function(inputs, barcode) {
  common_cells <- Reduce(intersect, list(rownames(inputs$expression), rownames(barcode), inputs$metadata$cell_id, inputs$true_tree$tip.label))
  if (length(common_cells) < 4) stop("Fewer than four shared cells among expression, barcode, metadata, and true tree. Check cell IDs.")

  true_tree <- ape::keep.tip(inputs$true_tree, common_cells)
  common_cells <- true_tree$tip.label

  expression <- inputs$expression[common_cells, , drop = FALSE]
  barcode <- barcode[common_cells, , drop = FALSE]
  metadata <- inputs$metadata[match(common_cells, inputs$metadata$cell_id), , drop = FALSE]

  list(expression = expression, barcode = barcode, metadata = metadata, true_tree = true_tree,
       ncells = length(common_cells), N_char = ncol(barcode))
}

parse_case_indices <- function(case_indices, case_ids) {
  n <- length(case_ids)
  if (length(case_indices) == 1 && identical(tolower(as.character(case_indices)), "all")) return(seq_len(n))

  if (is.numeric(case_indices)) {
    idx <- as.integer(case_indices)
    return(idx[idx >= 1 & idx <= n])
  }

  case_indices_chr <- as.character(case_indices)

  # Single numeric string, e.g. "5"
  if (length(case_indices_chr) == 1 && grepl("^\\d+$", case_indices_chr)) {
    idx <- as.integer(case_indices_chr)
    return(idx[idx >= 1 & idx <= n])
  }

  # Range, e.g. "5:5" or "1:10"
  if (length(case_indices_chr) == 1 && grepl("^\\d+\\s*:\\s*\\d+$", case_indices_chr)) {
    parts <- strsplit(case_indices_chr, ":", fixed = TRUE)[[1]]
    idx <- seq.int(as.integer(trimws(parts[1])), as.integer(trimws(parts[2])))
    return(idx[idx >= 1 & idx <= n])
  }

  # Comma separated numeric string, e.g. "1,3,5"
  if (length(case_indices_chr) == 1 && grepl(",", case_indices_chr)) {
    parts <- trimws(strsplit(case_indices_chr, ",", fixed = TRUE)[[1]])
    if (all(grepl("^\\d+$", parts))) {
      idx <- as.integer(parts)
      return(idx[idx >= 1 & idx <= n])
    }
  }

  # Otherwise treat as case IDs.
  idx <- match(case_indices_chr, case_ids)
  idx[!is.na(idx)]
}


find_case_dirs <- function(input_root) {
  if (!dir.exists(input_root)) stop("METHOD_INPUT_ROOT does not exist: ", input_root)
  candidate_files <- list.files(input_root, pattern = "^expression\\.csv$", recursive = TRUE, full.names = TRUE)
  dirs <- unique(dirname(candidate_files))
  dirs <- dirs[file.exists(file.path(dirs, "barcode_observed.csv"))]
  dirs <- sort(normalizePath(dirs, winslash = "/", mustWork = TRUE))
  if (length(dirs) == 0) stop("No FateScape case directories found under: ", input_root)
  dirs
}

sanitize_phylo_for_eval <- function(tree, tree_name = "tree") {
  if (is.null(tree)) stop(tree_name, " is NULL.")
  if (!inherits(tree, "phylo")) stop(tree_name, " is not a phylo object.")
  if (is.null(tree$edge)) stop(tree_name, "$edge is NULL.")

  edge <- tree$edge
  if (is.data.frame(edge)) edge <- as.matrix(edge)
  if (!is.matrix(edge)) edge <- matrix(edge, ncol = 2)
  if (ncol(edge) != 2) stop(tree_name, "$edge must have exactly two columns.")

  suppressWarnings(storage.mode(edge) <- "integer")
  edge <- edge[stats::complete.cases(edge), , drop = FALSE]
  tree$edge <- edge

  tree$tip.label <- trimws(as.character(tree$tip.label))
  empty <- is.na(tree$tip.label) | !nzchar(tree$tip.label)
  if (any(empty)) tree <- ape::drop.tip(tree, tree$tip.label[empty])

  if (anyDuplicated(tree$tip.label)) {
    keep <- !duplicated(tree$tip.label)
    drop_labels <- tree$tip.label[!keep]
    tree <- tryCatch(ape::drop.tip(tree, drop_labels), error = function(e) tree)
  }

  tree$edge.length <- NULL
  tree$node.label <- NULL

  ntip <- length(tree$tip.label)
  internal_nodes <- unique(tree$edge[, 1])
  internal_nodes <- internal_nodes[!is.na(internal_nodes) & internal_nodes > ntip]
  tree$Nnode <- length(internal_nodes)

  tree <- tryCatch(ape::collapse.singles(tree), error = function(e) tree)
  tree <- tryCatch(ape::unroot(tree), error = function(e) tree)
  tree <- tryCatch(ape::reorder.phylo(tree, order = "cladewise"), error = function(e) tree)
  tree
}

safe_write_reconstructed_tree <- function(tree, out_dir, filename = "FateScape_reconstructed_tree.nwk") {
  tree_file <- file.path(out_dir, filename)
  tree_to_save <- sanitize_phylo_for_eval(tree, "Tree_Merge")
  ape::write.tree(tree_to_save, file = tree_file)
  message("Reconstructed tree saved to: ", tree_file)
  tree_file
}

compute_tree_metrics <- function(true_tree, pred_tree) {
  true_tree <- sanitize_phylo_for_eval(true_tree, "true_tree")
  pred_tree <- sanitize_phylo_for_eval(pred_tree, "pred_tree")

  common <- intersect(true_tree$tip.label, pred_tree$tip.label)
  if (length(common) < 4) stop("Too few common tips for tree evaluation: ", length(common))

  true_tree2 <- ape::keep.tip(true_tree, common)
  pred_tree2 <- ape::keep.tip(pred_tree, common)
  true_tree2 <- sanitize_phylo_for_eval(true_tree2, "true_tree2")
  pred_tree2 <- sanitize_phylo_for_eval(pred_tree2, "pred_tree2")

  # Deterministically resolve polytomies for metrics that are sensitive to multifurcation.
  true_tree_bin <- tryCatch(ape::multi2di(true_tree2, random = FALSE), error = function(e) true_tree2)
  pred_tree_bin <- tryCatch(ape::multi2di(pred_tree2, random = FALSE), error = function(e) pred_tree2)
  true_tree_bin <- sanitize_phylo_for_eval(true_tree_bin, "true_tree_bin")
  pred_tree_bin <- sanitize_phylo_for_eval(pred_tree_bin, "pred_tree_bin")

  nye <- tryCatch(as.numeric(TreeDist::NyeSimilarity(true_tree_bin, pred_tree_bin, normalize = TRUE)), error = function(e) NA_real_)
  td <- tryCatch(as.numeric(TreeDist::TreeDistance(true_tree_bin, pred_tree_bin)), error = function(e) NA_real_)
  rf_norm <- tryCatch(as.numeric(phangorn::RF.dist(true_tree_bin, pred_tree_bin, normalize = TRUE)), error = function(e) NA_real_)
  rf_raw <- tryCatch(as.numeric(phangorn::RF.dist(true_tree_bin, pred_tree_bin, normalize = FALSE)), error = function(e) NA_real_)

  data.frame(
    n_tips = length(common),
    Nye = nye,
    TreeDistance = td,
    RF_normalized = rf_norm,
    RF_raw = rf_raw,
    stringsAsFactors = FALSE
  )
}



# -------------------------
# 3. MONOCLE3 INFERENCE
# -------------------------
infer_monocle_state <- function(expression, cell_ids, out_dir = NULL) {
  expr_mat <- as.matrix(expression)
  storage.mode(expr_mat) <- "numeric"
  expr_mat[is.na(expr_mat)] <- 0
  expr_mat[expr_mat < 0] <- 0

  # Monocle3 expects genes x cells.
  counts_gene_by_cell <- t(expr_mat)
  colnames(counts_gene_by_cell) <- cell_ids
  gene_ids <- rownames(counts_gene_by_cell)
  if (is.null(gene_ids) || any(gene_ids == "")) gene_ids <- paste0("gene_", seq_len(nrow(counts_gene_by_cell)))
  rownames(counts_gene_by_cell) <- gene_ids

  cell_md <- data.frame(cell_id = cell_ids, row.names = cell_ids, stringsAsFactors = FALSE)
  gene_md <- data.frame(gene_short_name = gene_ids, row.names = gene_ids, stringsAsFactors = FALSE)

  cds <- monocle3::new_cell_data_set(
    counts_gene_by_cell,
    cell_metadata = cell_md,
    gene_metadata = gene_md
  )

  # A small PCA is also computed outside Monocle for robust fallback ordering.
  pca <- stats::prcomp(log1p(expr_mat), scale. = FALSE)
  pc1 <- pca$x[, 1]

  monocle_error <- NA_character_
  fallback_used <- FALSE

  cds <- tryCatch({
    cds <- monocle3::preprocess_cds(
      cds,
      num_dim = min(MONOCLE_NUM_DIM, ncol(expr_mat), nrow(expr_mat) - 1L),
      method = MONOCLE_PREPROCESS_METHOD
    )
    cds <- monocle3::reduce_dimension(cds, reduction_method = MONOCLE_REDUCTION, preprocess_method = MONOCLE_PREPROCESS_METHOD)
    cds <- monocle3::cluster_cells(cds, reduction_method = MONOCLE_REDUCTION, k = min(MONOCLE_K, max(2L, floor(length(cell_ids) / 5))), resolution = MONOCLE_RESOLUTION)
    cds <- monocle3::learn_graph(cds, use_partition = TRUE)

    # Root-cell selection. Use expression-derived PC1/UMAP1 to avoid using true states.
    root_cell <- NULL
    if (MONOCLE_ROOT_MODE == "umap1_min") {
      red <- SingleCellExperiment::reducedDims(cds)[[MONOCLE_REDUCTION]]
      root_cell <- rownames(red)[which.min(red[, 1])]
    } else if (MONOCLE_ROOT_MODE == "largest_low_pc1_cluster") {
      cl_tmp <- as.character(monocle3::clusters(cds, reduction_method = MONOCLE_REDUCTION))
      cl_tab <- sort(table(cl_tmp), decreasing = TRUE)
      candidate_clusters <- names(cl_tab)[seq_len(min(5, length(cl_tab)))]
      means <- tapply(pc1, cl_tmp, mean)
      root_cl <- candidate_clusters[which.min(means[candidate_clusters])]
      root_cell <- names(which.min(ifelse(cl_tmp == root_cl, pc1, Inf)))
    } else {
      root_cell <- names(which.min(pc1))
    }
    cds <- monocle3::order_cells(cds, root_cells = root_cell)
    cds
  }, error = function(e) {
    monocle_error <<- conditionMessage(e)
    if (!ALLOW_MONOCLE_FALLBACK) stop(e)
    fallback_used <<- TRUE
    cds
  })

  # Extract clusters; if unavailable or too few, use kmeans on PCA fallback.
  cl <- tryCatch(as.character(monocle3::clusters(cds, reduction_method = MONOCLE_REDUCTION)), error = function(e) NULL)
  if (is.null(cl) || length(cl) != length(cell_ids) || length(unique(cl)) < MONOCLE_MIN_CLUSTERS) {
    fallback_used <- TRUE
    k <- min(max(MONOCLE_MIN_CLUSTERS, 3L), max(2L, floor(length(cell_ids) / 10)))
    km <- stats::kmeans(pca$x[, 1:min(2, ncol(pca$x)), drop = FALSE], centers = k, nstart = 10)
    cl <- as.character(km$cluster)
  }
  names(cl) <- cell_ids

  # Extract pseudotime. If unavailable/infinite, use PC1-rank fallback.
  pseudotime_vec <- tryCatch(monocle3::pseudotime(cds), error = function(e) rep(NA_real_, length(cell_ids)))
  pseudotime_vec <- as.numeric(pseudotime_vec[cell_ids])
  if (all(!is.finite(pseudotime_vec))) {
    fallback_used <- TRUE
    pseudotime_vec <- rank(pc1, ties.method = "average")
  } else {
    finite_pt <- is.finite(pseudotime_vec)
    if (any(!finite_pt)) {
      pseudotime_vec[!finite_pt] <- max(pseudotime_vec[finite_pt], na.rm = TRUE) + rank(pc1[!finite_pt], ties.method = "average")
    }
  }

  # Define state lineages from Monocle clusters ordered by mean pseudotime.
  # If Monocle partitions are available, define one lineage per partition.
  partitions <- tryCatch(as.character(monocle3::partitions(cds, reduction_method = MONOCLE_REDUCTION)), error = function(e) rep("1", length(cell_ids)))
  if (length(partitions) != length(cell_ids)) partitions <- rep("1", length(cell_ids))
  names(partitions) <- cell_ids

  state_lineages <- list()
  for (ptn in unique(partitions)) {
    cells_p <- cell_ids[partitions[cell_ids] == ptn]
    if (length(cells_p) == 0) next
    cl_p <- unique(cl[cells_p])
    means <- tapply(pseudotime_vec[cells_p], cl[cells_p], mean, na.rm = TRUE)
    cl_order <- names(sort(means))
    cl_order <- cl_order[!is.na(cl_order)]
    if (length(cl_order) >= 1L) state_lineages[[paste0("L", length(state_lineages) + 1L)]] <- as.character(cl_order)
  }
  if (length(state_lineages) == 0L) {
    means <- tapply(pseudotime_vec, cl, mean, na.rm = TRUE)
    state_lineages <- list(L1 = as.character(names(sort(means))))
  }

  state <- data.frame(cell_id = cell_ids, cluster = as.character(cl[cell_ids]), stringsAsFactors = FALSE)

  pseudotime_df <- data.frame(cell_id = cell_ids, pseudotime = pseudotime_vec, monocle_cluster = cl[cell_ids], monocle_partition = partitions[cell_ids], stringsAsFactors = FALSE)

  red <- tryCatch(as.data.frame(SingleCellExperiment::reducedDims(cds)[[MONOCLE_REDUCTION]]), error = function(e) NULL)
  if (is.null(red)) {
    red <- as.data.frame(pca$x[, 1:min(2, ncol(pca$x)), drop = FALSE])
    names(red) <- paste0("PC", seq_len(ncol(red)))
  }
  red$cell_id <- rownames(red)
  if (is.null(red$cell_id) || any(is.na(red$cell_id))) red$cell_id <- cell_ids
  red$monocle_cluster <- cl[red$cell_id]
  red <- red[, c("cell_id", setdiff(names(red), "cell_id")), drop = FALSE]

  cluster_table <- data.frame(cluster = names(table(cl)), n_cells = as.integer(table(cl)), stringsAsFactors = FALSE)
  cluster_pt <- aggregate(pseudotime_vec, by = list(cluster = cl), FUN = mean)
  names(cluster_pt)[2] <- "mean_pseudotime"
  cluster_table <- merge(cluster_table, cluster_pt, by = "cluster", all.x = TRUE)

  if (!is.null(out_dir)) {
    safe_write_csv(state, file.path(out_dir, "monocle_state_assignments.csv"), row.names = FALSE)
    saveRDS(state_lineages, file.path(out_dir, "monocle_state_lineages.rds"))
    safe_write_csv(pseudotime_df, file.path(out_dir, "monocle_pseudotime.csv"), row.names = FALSE)
    safe_write_csv(red, file.path(out_dir, "monocle_reduced_dims.csv"), row.names = FALSE)
    safe_write_csv(cluster_table, file.path(out_dir, "monocle_cluster_table.csv"), row.names = FALSE)
    saveRDS(cds, file.path(out_dir, "monocle_cds.rds"))
  }

  list(
    state = state,
    state_lineages = state_lineages,
    cds = cds,
    pseudotime = pseudotime_df,
    reduced_dims = red,
    cluster_table = cluster_table,
    monocle_error = monocle_error,
    fallback_used = fallback_used,
    n_clusters = length(unique(cl)),
    n_lineages = length(state_lineages)
  )
}

# -------------------------
# 4. RUN ONE CASE
# -------------------------
run_fatescape_monocle_one_case <- function(case_dir, out_dir) {
  case_id <- basename(case_dir)
  safe_dir_create(out_dir)

  message("\n============================================================")
  message("[FateScape-Monocle] Running case: ", case_id)
  message("Input : ", case_dir)
  message("Output: ", out_dir)

  inputs <- load_case_inputs(case_dir)
  if (USE_PREIMPUTED_BARCODE && !is.null(inputs$barcode_imputed)) {
    barcode_for_reconstruction <- inputs$barcode_imputed
    barcode_source <- "preimputed"
  } else {
    barcode_obs <- inputs$barcode_observed
    N_char <- ncol(barcode_obs)
    ncells <- nrow(barcode_obs)
    r_n <- ceiling(N_char * LABEL_SITE_FRACTION)
    barcode_for_reconstruction <- dropout_imputation_alter(barcode_obs, N_char, ncells, r_n)
    barcode_source <- "imputed_in_pipeline"
  }

  aligned <- align_expression_barcode_tree(inputs, barcode_for_reconstruction)
  message("Cells: ", aligned$ncells, "; barcode sites: ", aligned$N_char, "; barcode source: ", barcode_source)

  # Infer state labels and lineages from expression using Monocle3.
  mono <- infer_monocle_state(aligned$expression, aligned$true_tree$tip.label, out_dir = out_dir)
  message("Monocle clusters: ", mono$n_clusters, "; lineages: ", mono$n_lineages, "; fallback: ", mono$fallback_used)

  start_time <- proc.time()

  state_lineages <- mono$state_lineages
  state <- mono$state
  barcodes <- aligned$barcode
  ncells <- aligned$ncells
  N_char <- aligned$N_char

  # Group cells by Monocle3-inferred state lineage.
  sl_info <- state_lineage_info(state_lineages, ncells, state, barcodes)
  state_lineage <- sl_info[[1]]
  cell_lineages <- sl_info[[2]]
  barcodes_lineages <- sl_info[[3]]
  state_labels_lineages <- sl_info[[4]]

  Trees_initial <- initial_tree_construction(state_lineages, barcodes_lineages)

  LiterateRefine <- subtree_refinement(
    Trees_initial,
    state_lineages,
    barcodes_lineages,
    N_char,
    state_labels_lineages,
    lambda1 = LAMBDA_STATE,
    lambda2 = LAMBDA_BARCODE,
    maxIter = MAX_ITER,
    repeat_time = REPEAT_TIME
  )
  bestsubtree <- LiterateRefine[[1]]

  bestsubtree <- drop_duplicated_tips(
    bestsubtree,
    barcodes_lineages,
    cell_lineages,
    state_lineage,
    alpha = DROP_DUP_ALPHA,
    beta = DROP_DUP_BETA
  )

  subtrees_rootbar <- get_subtree_root_barcodes(bestsubtree, state, barcodes, length(state_lineages))
  CM <- common_mutation_matrix(subtrees_rootbar)
  rank_res <- subtrees_rank(CM)
  Nodes_rank <- rank_res[[1]]
  Nodes_weight <- rank_res[[2]]
  subsubtrees <- decompose_subtrees(bestsubtree, state, barcodes, length(state_lineages))

  if (exists("merge_subcell_trees_ordered", mode = "function")) {
    Tree_Merge <- merge_subcell_trees_ordered(CM, Nodes_rank, Nodes_weight, subsubtrees, bind_tree_list = list(), nsubtree = 1)
  } else if (exists("merge_subcell_trees_ward", mode = "function")) {
    Tree_Merge <- merge_subcell_trees_ward(subsubtrees)
  } else {
    stop("No supported subtree merging function found.")
  }

  runtime <- proc.time() - start_time

  # ------------------------------------------------------------
  # Save reconstructed tree immediately after inference.
  # This does not change the FateScape/Monocle reconstruction algorithm;
  # it only prevents a downstream evaluation error from losing the inferred tree.
  # ------------------------------------------------------------
  tree_file <- tryCatch(
    safe_write_reconstructed_tree(Tree_Merge, out_dir, "FateScape_Monocle_reconstructed_tree.nwk"),
    error = function(e) {
      warning("Failed to save reconstructed tree: ", conditionMessage(e))
      NA_character_
    }
  )

  # Evaluate only after tree storage. If evaluation fails, still keep the tree file
  # and write metrics with eval_error.
  metrics <- tryCatch(
    compute_tree_metrics(aligned$true_tree, Tree_Merge),
    error = function(e) {
      warning("Evaluation failed but reconstructed tree was already saved: ", conditionMessage(e))
      data.frame(
        n_tips = length(Tree_Merge$tip.label),
        Nye = NA_real_,
        TreeDistance = NA_real_,
        RF_normalized = NA_real_,
        RF_raw = NA_real_,
        eval_error = conditionMessage(e),
        stringsAsFactors = FALSE
      )
    }
  )
  metrics$case_id <- case_id
  metrics$method <- "FateScape-Monocle"
  metrics$barcode_source <- barcode_source
  metrics$runtime_sec <- as.numeric(runtime[["elapsed"]])
  metrics$lambda_barcode <- LAMBDA_BARCODE
  metrics$lambda_state <- LAMBDA_STATE
  metrics$max_iter <- MAX_ITER
  metrics$repeat_time <- REPEAT_TIME
  metrics$monocle_n_clusters <- mono$n_clusters
  metrics$monocle_n_lineages <- mono$n_lineages
  metrics$monocle_fallback_used <- mono$fallback_used
  metrics$reconstructed_tree_file <- tree_file

  safe_write_csv(metrics, file.path(out_dir, "FateScape_Monocle_metrics.csv"), row.names = FALSE)

  run_params <- data.frame(
    case_id = case_id,
    case_dir = case_dir,
    output_dir = out_dir,
    barcode_source = barcode_source,
    ncells = ncells,
    N_char = N_char,
    lambda_barcode = LAMBDA_BARCODE,
    lambda_state = LAMBDA_STATE,
    max_iter = MAX_ITER,
    repeat_time = REPEAT_TIME,
    drop_dup_alpha = DROP_DUP_ALPHA,
    drop_dup_beta = DROP_DUP_BETA,
    monocle_num_dim = MONOCLE_NUM_DIM,
    monocle_reduction = MONOCLE_REDUCTION,
    monocle_preprocess_method = MONOCLE_PREPROCESS_METHOD,
    monocle_k = MONOCLE_K,
    monocle_resolution = MONOCLE_RESOLUTION,
    monocle_root_mode = MONOCLE_ROOT_MODE,
    monocle_n_clusters = mono$n_clusters,
    monocle_n_lineages = mono$n_lineages,
    monocle_fallback_used = mono$fallback_used,
    runtime_sec = as.numeric(runtime[["elapsed"]]),
    reconstructed_tree_file = tree_file,
    monocle_note = mono$monocle_error,
    stringsAsFactors = FALSE
  )
  safe_write_csv(run_params, file.path(out_dir, "FateScape_Monocle_run_params.csv"), row.names = FALSE)

  if (SAVE_INTERMEDIATE_OBJECTS) {
    saveRDS(
      list(
        aligned = aligned,
        monocle = mono,
        state_lineage = state_lineage,
        cell_lineages = cell_lineages,
        barcodes_lineages = barcodes_lineages,
        state_labels_lineages = state_labels_lineages,
        Trees_initial = Trees_initial,
        LiterateRefine = LiterateRefine,
        bestsubtree = bestsubtree,
        subtrees_rootbar = subtrees_rootbar,
        CM = CM,
        subsubtrees = subsubtrees,
        Tree_Merge = Tree_Merge,
        metrics = metrics
      ),
      file = file.path(out_dir, "FateScape_Monocle_intermediate_objects.rds")
    )
  }

  list(status = "success", metrics = metrics, runtime_sec = as.numeric(runtime[["elapsed"]]))
}

# -------------------------
# 5. BATCH DRIVER
# -------------------------
METHOD_INPUT_ROOT <- normalizePath(METHOD_INPUT_ROOT, winslash = "/", mustWork = TRUE)
METHOD_OUTPUT_ROOT <- safe_dir_create(METHOD_OUTPUT_ROOT)
SUMMARY_DIR <- safe_dir_create(file.path(METHOD_OUTPUT_ROOT, "_batch_summary"))

case_dirs <- find_case_dirs(METHOD_INPUT_ROOT)
case_ids <- basename(case_dirs)
selected_idx <- parse_case_indices(CASE_INDICES, case_ids)
if (length(selected_idx) == 0) stop("CASE_INDICES did not select any valid cases.")
selected_case_dirs <- case_dirs[selected_idx]
selected_case_ids <- case_ids[selected_idx]

batch_plan <- data.frame(
  case_index = selected_idx,
  case_id = selected_case_ids,
  case_dir = selected_case_dirs,
  output_dir = file.path(METHOD_OUTPUT_ROOT, selected_case_ids),
  stringsAsFactors = FALSE
)

safe_write_csv(batch_plan, file.path(SUMMARY_DIR, "FateScape_Monocle_batch_plan.csv"), row.names = FALSE)
message("Batch plan written to: ", file.path(SUMMARY_DIR, "FateScape_Monocle_batch_plan.csv"))
message("Number of selected cases: ", nrow(batch_plan))

batch_records <- list()
metrics_records <- list()
error_records <- list()

for (i in seq_len(nrow(batch_plan))) {
  case_id <- batch_plan$case_id[i]
  case_dir <- batch_plan$case_dir[i]
  out_dir <- batch_plan$output_dir[i]
  metrics_file <- file.path(out_dir, "FateScape_Monocle_metrics.csv")

  if (SKIP_FINISHED && !FORCE_RERUN && file.exists(metrics_file)) {
    message("\n[Skip] Existing metrics found for case: ", case_id)
    m <- read.csv(metrics_file, stringsAsFactors = FALSE)
    metrics_records[[length(metrics_records) + 1L]] <- m
    batch_records[[length(batch_records) + 1L]] <- data.frame(case_id = case_id, case_dir = case_dir, output_dir = out_dir, status = "skipped_existing", error = NA_character_, stringsAsFactors = FALSE)
    next
  }

  result <- tryCatch(
    run_fatescape_monocle_one_case(case_dir, out_dir),
    error = function(e) {
      msg <- conditionMessage(e)
      message("[ERROR] ", case_id, ": ", msg)
      if (SAVE_FAILED_INTERMEDIATE_OBJECTS) {
        safe_dir_create(out_dir)
        saveRDS(list(error = msg, case_dir = case_dir), file = file.path(out_dir, "FateScape_Monocle_failed_case.rds"))
      }
      if (STOP_ON_ERROR) stop(e)
      list(status = "failed", error = msg)
    }
  )

  if (identical(result$status, "success")) {
    metrics_records[[length(metrics_records) + 1L]] <- result$metrics
    batch_records[[length(batch_records) + 1L]] <- data.frame(case_id = case_id, case_dir = case_dir, output_dir = out_dir, status = "success", error = NA_character_, stringsAsFactors = FALSE)
  } else {
    err <- result$error
    error_records[[length(error_records) + 1L]] <- data.frame(case_id = case_id, case_dir = case_dir, output_dir = out_dir, error = err, stringsAsFactors = FALSE)
    batch_records[[length(batch_records) + 1L]] <- data.frame(case_id = case_id, case_dir = case_dir, output_dir = out_dir, status = "failed", error = err, stringsAsFactors = FALSE)
  }

  safe_write_csv(do.call(rbind, batch_records), file.path(SUMMARY_DIR, "FateScape_Monocle_batch_manifest.csv"), row.names = FALSE)
  if (length(metrics_records) > 0) safe_write_csv(do.call(rbind, metrics_records), file.path(SUMMARY_DIR, "FateScape_Monocle_batch_metrics_all.csv"), row.names = FALSE)
  if (length(error_records) > 0) safe_write_csv(do.call(rbind, error_records), file.path(SUMMARY_DIR, "FateScape_Monocle_batch_errors.csv"), row.names = FALSE)
}

message("\nDone.")
message("Batch summary directory: ", SUMMARY_DIR)
if (length(metrics_records) > 0) message("Metrics: ", file.path(SUMMARY_DIR, "FateScape_Monocle_batch_metrics_all.csv"))
if (length(error_records) > 0) message("Errors : ", file.path(SUMMARY_DIR, "FateScape_Monocle_batch_errors.csv"))
