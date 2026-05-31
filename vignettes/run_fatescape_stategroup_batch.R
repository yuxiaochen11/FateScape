#!/usr/bin/env Rscript


# -------------------------
# 0. USER SETTINGS
# -------------------------
METHOD_INPUT_ROOT <- "FateScape/data/Sensitivity_analysis"
METHOD_OUTPUT_ROOT <- "FateScape/results/Sensitivity_analysis_stategroup"

# Optional: if FateScape is not installed as an R package, set package root here.
# Example: FATESCAPE_LOCAL_PATH <- "D:/projects/lineage tracing/FateScape"

CASE_INDICES <- "all"

# If TRUE, do not rerun cases with an existing FateScape_StateGroup_metrics.csv unless FORCE_RERUN=TRUE.
SKIP_FINISHED <- TRUE
FORCE_RERUN <- FALSE
STOP_ON_ERROR <- FALSE

# Reconstruction parameters.
USE_PREIMPUTED_BARCODE <- TRUE
LABEL_SITE_FRACTION <- 0.70
LAMBDA_BARCODE <- 0.90
LAMBDA_STATE <- 0.10
MAX_ITER <- 100
REPEAT_TIME <- 10
DROP_DUP_ALPHA <- 1.5
DROP_DUP_BETA <- 1.5

# StateGroup-specific handling:
# Cell-state groups with fewer than this many cells are not sent into
# FateScape subtree refinement/integration, because the original refinement
# assumes enough internal nodes for subtree swapping. Instead, their cells are
# bound directly to the merged tree root as individual tips.
MIN_STATEGROUP_TREE_CELLS <- 3

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
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  devtools::load_all(FATESCAPE_LOCAL_PATH)
} else {
  library(FateScape)
}

suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
  library(TreeTools)
  library(TreeDist)
})

set.seed(GLOBAL_SEED)

# Explicitly source patched core files after load_all/library to avoid stale installed namespace.
source_patched_fatescape_core <- function(package_root) {
  if (is.null(package_root) || !dir.exists(package_root)) return(invisible(FALSE))
  rdir <- file.path(package_root, "R")
  core_files <- c(
    "utils.R",
    "barcode_imputation.R",
    "likelihood_calculation.R",
    "subtree_refinement.R",
    "subtree_integration.R"
  )
  for (ff in core_files) {
    f <- file.path(rdir, ff)
    if (file.exists(f)) {
      source(f, local = .GlobalEnv)
      message("Sourced patched core file: ", normalizePath(f, winslash = "/", mustWork = FALSE))
    } else {
      warning("Patched core file not found: ", f)
    }
  }
  if (exists("improved_hamming_distance", mode = "function")) {
    test1 <- improved_hamming_distance(c(NA, 1, 0), c(1, 1, 0), alpha = 0, beta = 1)
    test2 <- improved_hamming_distance(c(NA, 1, 0), c(NA, 1, 0), alpha = 0, beta = 1)
    message("Missing-safe hamming test passed: ", test1, ", ", test2)
  }
  invisible(TRUE)
}

source_patched_fatescape_core(FATESCAPE_LOCAL_PATH)

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
    v[is.na(v) | v %in% c("", "NA", "NaN", "nan", "NULL", "null", "None", "none", "?", "-1")] <- "-"
    v
  })
  as.matrix(x)
}

load_case_inputs <- function(case_dir) {
  files <- list(
    expression = file.path(case_dir, "expression.csv"),
    barcode_observed = file.path(case_dir, "barcode_observed.csv"),
    barcode_imputed = file.path(case_dir, "barcode_imputed.csv"),
    metadata = file.path(case_dir, "metadata.csv"),
    true_tree = file.path(case_dir, "true_tree.nwk"),
    state_tree = file.path(case_dir, "state_tree.nwk"),
    state_lineages = file.path(case_dir, "state_lineages.rds")
  )

  # Compatibility with canonical folder names or older converter names.
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
  if (file.exists(files$barcode_imputed)) {
    barcode_imputed <- standardize_barcode(read_table_with_rownames(files$barcode_imputed))
  }

  metadata <- read.csv(files$metadata, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"cell_id" %in% names(metadata)) {
    if ("cell" %in% names(metadata)) {
      metadata$cell_id <- metadata$cell
    } else {
      stop("metadata must contain cell_id")
    }
  }
  if (!"cluster" %in% names(metadata)) {
    if ("state_label" %in% names(metadata)) {
      metadata$cluster <- as.character(metadata$state_label)
    } else if ("state" %in% names(metadata)) {
      metadata$cluster <- as.character(metadata$state)
    } else {
      stop("metadata must contain cluster, state_label, or state")
    }
  }
  state <- data.frame(
    cell_id = metadata$cell_id,
    cluster = as.character(metadata$cluster),
    stringsAsFactors = FALSE
  )

  true_tree <- ape::read.tree(files$true_tree)
  state_tree <- if (file.exists(files$state_tree)) ape::read.tree(files$state_tree) else NULL
  state_lineages <- if (file.exists(files$state_lineages)) readRDS(files$state_lineages) else NULL

  list(
    expression = expression,
    barcode_observed = barcode_observed,
    barcode_imputed = barcode_imputed,
    metadata = metadata,
    state = state,
    true_tree = true_tree,
    state_tree = state_tree,
    state_lineages = state_lineages,
    files = files
  )
}

infer_default_state_lineages <- function(state_tree) {
  if (is.null(state_tree)) {
    stop("state_lineages.rds is missing and state_tree.nwk is unavailable; cannot infer state lineages.")
  }
  n_tips <- length(state_tree$tip.label)
  terminal_nodes <- seq_len(n_tips)
  root <- length(state_tree$tip.label) + 1L
  lineages <- list()
  for (tip in terminal_nodes) {
    path <- ape::nodepath(state_tree, from = root, to = tip)[[1]]
    lineages[[paste0("L", tip)]] <- path
  }
  lineages
}

make_state_group_lineages <- function(state) {
  # Reviewer-requested ablation: group cells by discrete cell-state clusters.
  #
  # IMPORTANT: state_lineage_info() internally indexes state_lineages by
  # "L1", "L2", ... rather than by arbitrary names. Therefore the state-group
  # ablation must name groups L1, L2, ... . Otherwise every group becomes empty
  # and initial_tree_construction()/as.phylo.hclust() fails with
  # "data object must contain taxa names".
  if (!"cluster" %in% names(state)) stop("state must contain a cluster column")
  states <- unique(as.character(state$cluster))
  states <- states[!is.na(states) & nzchar(states)]
  if (length(states) == 0) stop("No valid cell-state labels found.")

  suppressWarnings(num_states <- as.numeric(states))
  if (all(!is.na(num_states))) {
    states <- states[order(num_states)]
  } else {
    states <- sort(states)
  }

  out <- lapply(states, function(s) as.character(s))
  names(out) <- paste0("L", seq_along(states))
  attr(out, "state_group_labels") <- states
  out
}

quote_tip_for_newick <- function(x) {
  x <- as.character(x)
  if (grepl("[\\s,:;()\\[\\]]", x)) {
    x <- paste0("'", gsub("'", "''", x), "'")
  }
  x
}

make_one_tip_phylo <- function(tip_label) {
  # Construct a valid ape::phylo object for a singleton group. This is only used
  # when a cell-state group contains one cell; it does not change the original
  # FateScape algorithm for multi-cell groups.
  tip_label <- as.character(tip_label)
  tree <- list(
    edge = matrix(c(2L, 1L), ncol = 2, byrow = TRUE),
    tip.label = tip_label,
    Nnode = 1L
  )
  class(tree) <- "phylo"
  tree
}

make_two_tip_phylo <- function(tip_labels) {
  tip_labels <- as.character(tip_labels)
  txt <- paste0("(", paste(vapply(tip_labels, quote_tip_for_newick, character(1)), collapse = ","), ");")
  ape::read.tree(text = txt)
}

make_star_phylo <- function(tip_labels) {
  tip_labels <- unique(as.character(tip_labels))
  tip_labels <- tip_labels[!is.na(tip_labels) & nzchar(tip_labels)]
  if (length(tip_labels) == 0) stop("Cannot construct a star tree with zero tips.")
  if (length(tip_labels) == 1) return(make_one_tip_phylo(tip_labels[1]))
  if (length(tip_labels) == 2) return(make_two_tip_phylo(tip_labels))
  txt <- paste0("(", paste(vapply(tip_labels, quote_tip_for_newick, character(1)), collapse = ","), ");")
  ape::read.tree(text = txt)
}

bind_direct_tips_to_root <- function(tree, direct_tips) {
  direct_tips <- unique(as.character(direct_tips))
  direct_tips <- direct_tips[!is.na(direct_tips) & nzchar(direct_tips)]
  if (length(direct_tips) == 0) return(tree)

  if (is.null(tree) || !inherits(tree, "phylo") || is.null(tree$tip.label) || length(tree$tip.label) == 0) {
    return(make_star_phylo(direct_tips))
  }

  direct_tips <- setdiff(direct_tips, as.character(tree$tip.label))
  if (length(direct_tips) == 0) return(tree)

  # Attach each small-group cell independently to the current root. This follows
  # the intended StateGroup ablation: small groups are not refined as subtrees,
  # but are still retained as individual leaves in the final tree.
  tree$tip.label <- as.character(tree$tip.label)
  tree$edge <- as.matrix(tree$edge)
  storage.mode(tree$edge) <- "integer"
  if (is.null(tree$edge.length)) tree$edge.length <- rep(1, nrow(tree$edge))

  for (tip in direct_tips) {
    root_node <- length(tree$tip.label) + 1L
    tree <- tryCatch(
      ape::bind.tip(tree, tip.label = tip, where = root_node, position = 0),
      error = function(e) {
        warning("bind.tip failed for direct small-group tip ", tip,
                "; falling back to root star over existing tips plus remaining direct tips. Error: ",
                conditionMessage(e))
        make_star_phylo(c(tree$tip.label, tip, setdiff(direct_tips, tip)))
      }
    )
    if (is.null(tree$edge.length)) tree$edge.length <- rep(1, nrow(tree$edge))
  }

  tree$edge.length <- NULL
  tree$node.label <- NULL
  tree
}

initial_tree_construction_stategroup_safe <- function(state_lineages, barcodes_lineages) {
  # Same behavior as initial_tree_construction() for groups with >= 3 cells,
  # but robust to empty/singleton/two-cell state groups and missing rownames.
  trees_initial <- list()

  for (j in seq_along(state_lineages)) {
    lineage_label <- paste0("L", j)
    barcode_lineage <- barcodes_lineages[[lineage_label]]

    if (is.null(barcode_lineage)) {
      warning("Missing barcode matrix for ", lineage_label, "; skipping this group.")
      next
    }

    barcode_lineage <- as.matrix(barcode_lineage)
    n <- nrow(barcode_lineage)
    if (is.null(n) || n == 0) {
      warning("No cells found for ", lineage_label, "; skipping this group.")
      next
    }

    tip_labels <- rownames(barcode_lineage)
    if (is.null(tip_labels) || any(is.na(tip_labels)) || any(!nzchar(tip_labels))) {
      tip_labels <- paste0("cell_", seq_len(n))
      rownames(barcode_lineage) <- tip_labels
    }

    if (n == 1) {
      trees_initial[[lineage_label]] <- make_one_tip_phylo(tip_labels[1])
      next
    }

    if (n == 2) {
      trees_initial[[lineage_label]] <- make_two_tip_phylo(tip_labels)
      next
    }

    dist_matrix <- hamming_distance(barcode_lineage)
    if (is.null(attr(dist_matrix, "Labels"))) {
      attr(dist_matrix, "Labels") <- tip_labels
    }

    tree_WD <- nodes_clustering(dist_matrix)
    tree_WD$tip.label <- as.character(tree_WD$tip.label)
    trees_initial[[lineage_label]] <- tree_WD
  }

  if (length(trees_initial) == 0) {
    stop("No valid state-group initial trees were constructed. Check state labels and cell IDs.")
  }

  trees_initial
}

subtree_refinement_stategroup_safe <- function(Trees_initial, state_lineages, barcodes_lineages, N_char,
                                               state_labels_lineages, lambda1, lambda2,
                                               maxIter = 100, repeat_time = 10) {
  # StateGroup ablation can contain singleton or two-cell state groups. The original
  # subtree_refinement() assumes enough internal nodes for subtree swapping and can
  # fail with "subscript out of bounds" on these small groups. This wrapper keeps the
  # original refinement for groups with >=3 cells, but safely returns the initial tree
  # for small or failed groups.
  ptm <- proc.time()
  bestsubtree <- list()
  bestsubtreescore <- list()

  for (lineage_label in names(Trees_initial)) {
    tree <- Trees_initial[[lineage_label]]
    n_tip <- length(tree$tip.label)

    if (is.null(tree) || !inherits(tree, "phylo") || is.null(tree$edge) || n_tip < 3) {
      bestsubtree[[lineage_label]] <- tree
      bestsubtreescore[[lineage_label]] <- NA_real_
      next
    }

    # The original subtree_refinement() indexes lineages internally as L1, L2, ... .
    # Therefore, for a one-group refinement call, rename the current state group to L1,
    # run the original function, then map the result back to the original label.
    one_Trees <- list(L1 = tree)
    one_state_lineages <- list(L1 = state_lineages[[lineage_label]])
    one_barcodes <- list(L1 = barcodes_lineages[[lineage_label]])
    one_state_labels <- list(L1 = state_labels_lineages[[lineage_label]])

    refined <- tryCatch({
      subtree_refinement(
        one_Trees,
        one_state_lineages,
        one_barcodes,
        N_char,
        one_state_labels,
        lambda1 = lambda1,
        lambda2 = lambda2,
        maxIter = maxIter,
        repeat_time = repeat_time
      )
    }, error = function(e) {
      warning("StateGroup refinement failed for ", lineage_label,
              "; using initial local tree. Error: ", conditionMessage(e))
      NULL
    })

    if (is.null(refined) || is.null(refined[[1]]) || is.null(refined[[1]][["L1"]])) {
      bestsubtree[[lineage_label]] <- tree
      bestsubtreescore[[lineage_label]] <- NA_real_
    } else {
      bestsubtree[[lineage_label]] <- refined[[1]][["L1"]]
      score_val <- tryCatch(refined[[2]][["L1"]], error = function(e) NA_real_)
      bestsubtreescore[[lineage_label]] <- score_val
    }
  }

  list(bestsubtree = bestsubtree,
       bestsubtreescore = bestsubtreescore,
       total_time = proc.time() - ptm)
}


align_case <- function(inputs, barcode) {
  common_cells <- Reduce(intersect, list(rownames(barcode), inputs$state$cell_id, inputs$true_tree$tip.label))
  if (length(common_cells) < 4) {
    stop("Fewer than four shared cells among barcode, metadata, and true tree. Check cell IDs.")
  }

  barcode <- barcode[common_cells, , drop = FALSE]
  inputs$state <- inputs$state[match(common_cells, inputs$state$cell_id), , drop = FALSE]
  inputs$true_tree <- ape::keep.tip(inputs$true_tree, common_cells)

  # Reorder according to true tree tip labels.
  common_cells <- inputs$true_tree$tip.label
  barcode <- barcode[common_cells, , drop = FALSE]
  inputs$state <- inputs$state[match(common_cells, inputs$state$cell_id), , drop = FALSE]

  inputs$barcode <- barcode
  inputs$ncells <- nrow(barcode)
  inputs$N_char <- ncol(barcode)
  inputs
}

sanitize_phylo_for_eval <- function(tree, tree_name = "tree") {
  if (is.null(tree)) stop(tree_name, " is NULL.")
  if (!inherits(tree, "phylo")) stop(tree_name, " is not a phylo object.")
  if (is.null(tree$edge)) stop(tree_name, "$edge is NULL.")

  edge <- tree$edge
  if (is.data.frame(edge)) edge <- as.matrix(edge)
  if (!is.matrix(edge)) edge <- matrix(edge, ncol = 2)
  if (ncol(edge) != 2) stop(tree_name, "$edge must have exactly two columns.")
  edge <- edge[stats::complete.cases(edge), , drop = FALSE]
  storage.mode(edge) <- "integer"
  tree$edge <- edge

  tree$tip.label <- trimws(as.character(tree$tip.label))
  tree$tip.label[is.na(tree$tip.label)] <- ""

  tree$edge.length <- NULL
  tree$node.label <- NULL

  ntip <- length(tree$tip.label)
  internal_nodes <- unique(edge[, 1])
  internal_nodes <- internal_nodes[!is.na(internal_nodes) & internal_nodes > ntip]
  tree$Nnode <- max(1L, length(internal_nodes))

  tree <- tryCatch(ape::collapse.singles(tree), error = function(e) tree)
  tree <- tryCatch(ape::unroot(tree), error = function(e) tree)
  tree <- tryCatch(ape::reorder.phylo(tree, order = "cladewise"), error = function(e) tree)
  tree
}

save_reconstructed_tree_first <- function(tree, out_dir, file_name = "FateScape_StateGroup_reconstructed_tree.nwk") {
  safe_dir_create(out_dir)
  tree_file <- file.path(out_dir, file_name)
  tree_to_save <- tryCatch(
    sanitize_phylo_for_eval(tree, "Tree_Merge"),
    error = function(e) {
      warning("Tree sanitization before saving failed; trying raw tree write: ", conditionMessage(e))
      tree
    }
  )
  tree_to_save$edge.length <- NULL
  tree_to_save$node.label <- NULL
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

  true_tree_bin <- tryCatch(ape::multi2di(true_tree2, random = FALSE), error = function(e) true_tree2)
  pred_tree_bin <- tryCatch(ape::multi2di(pred_tree2, random = FALSE), error = function(e) pred_tree2)

  true_tree_bin <- sanitize_phylo_for_eval(true_tree_bin, "true_tree_bin")
  pred_tree_bin <- sanitize_phylo_for_eval(pred_tree_bin, "pred_tree_bin")

  rf_raw <- tryCatch(as.numeric(phangorn::RF.dist(true_tree_bin, pred_tree_bin, normalize = FALSE)), error = function(e) NA_real_)
  rf_norm <- tryCatch(as.numeric(phangorn::RF.dist(true_tree_bin, pred_tree_bin, normalize = TRUE)), error = function(e) NA_real_)
  nye <- tryCatch(as.numeric(TreeDist::NyeSimilarity(true_tree_bin, pred_tree_bin, normalize = TRUE)), error = function(e) NA_real_)
  td <- tryCatch(as.numeric(TreeDist::TreeDistance(true_tree_bin, pred_tree_bin)), error = function(e) NA_real_)

  data.frame(
    n_tips = length(common),
    Nye = nye,
    TreeDistance = td,
    RF_normalized = rf_norm,
    RF_raw = rf_raw,
    eval_error = NA_character_,
    stringsAsFactors = FALSE
  )
}
parse_case_indices <- function(case_indices, case_ids) {
  n <- length(case_ids)
  if (length(case_indices) == 1 && identical(tolower(as.character(case_indices)), "all")) return(seq_len(n))

  if (is.numeric(case_indices)) {
    idx <- as.integer(case_indices)
    return(idx[idx >= 1 & idx <= n])
  }

  case_indices_chr <- as.character(case_indices)

  if (length(case_indices_chr) == 1 && grepl("^\\d+$", case_indices_chr)) {
    idx <- as.integer(case_indices_chr)
    return(idx[idx >= 1 & idx <= n])
  }

  if (length(case_indices_chr) == 1 && grepl("^\\d+\\s*:\\s*\\d+$", case_indices_chr)) {
    parts <- strsplit(case_indices_chr, ":", fixed = TRUE)[[1]]
    idx <- seq.int(as.integer(trimws(parts[1])), as.integer(trimws(parts[2])))
    return(idx[idx >= 1 & idx <= n])
  }

  if (length(case_indices_chr) == 1 && grepl(",", case_indices_chr)) {
    parts <- trimws(strsplit(case_indices_chr, ",", fixed = TRUE)[[1]])
    if (all(grepl("^\\d+$", parts))) {
      idx <- as.integer(parts)
      return(idx[idx >= 1 & idx <= n])
    }
  }

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

# -------------------------
# 3. RUN ONE CASE
# -------------------------
run_fatescape_one_case <- function(case_dir, out_dir) {
  case_id <- basename(case_dir)
  safe_dir_create(out_dir)

  message("\n============================================================")
  message("[FateScape-StateGroup] Running case: ", case_id)
  message("Input : ", case_dir)
  message("Output: ", out_dir)

  inputs <- load_case_inputs(case_dir)
  # Do not use the true root-to-terminal state lineages in this ablation.
  # State-group strategy is constructed after cell alignment from metadata-derived state labels.

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

  inputs <- align_case(inputs, barcode_for_reconstruction)
  message("Cells: ", inputs$ncells, "; barcode sites: ", inputs$N_char, "; barcode source: ", barcode_source)

  start_time <- proc.time()

  state <- inputs$state
  barcodes <- inputs$barcode
  ncells <- inputs$ncells
  N_char <- inputs$N_char

  # Reviewer-requested alternative grouping: construct one local tree per discrete state.
  # This forms a non-overlapping partition and does not use the true state-lineage paths.
  state_lineages <- make_state_group_lineages(state)
  saveRDS(state_lineages, file = file.path(out_dir, "state_group_lineages.rds"))
  safe_write_csv(
    data.frame(cell_id = state$cell_id, cluster = state$cluster, stringsAsFactors = FALSE),
    file.path(out_dir, "state_group_assignments.csv"),
    row.names = FALSE
  )
  message("State-group lineages: ", length(state_lineages), " singleton groups.")

  # Group cells by discrete state groups.
  sl_info <- state_lineage_info(state_lineages, ncells, state, barcodes)
  state_lineage <- sl_info[[1]]
  cell_lineages <- sl_info[[2]]
  barcodes_lineages <- sl_info[[3]]
  state_labels_lineages <- sl_info[[4]]

  group_sizes <- vapply(seq_along(state_lineages), function(j) {
    lineage_label <- paste0("L", j)
    cell_key <- paste0("cell_", lineage_label)
    length(cell_lineages[[cell_key]])
  }, integer(1))
  safe_write_csv(
    data.frame(
      lineage = paste0("L", seq_along(state_lineages)),
      state_group = unname(unlist(state_lineages)),
      n_cells = group_sizes,
      stringsAsFactors = FALSE
    ),
    file.path(out_dir, "state_group_sizes.csv"),
    row.names = FALSE
  )
  message("State-group cell counts: min=", min(group_sizes), "; median=", stats::median(group_sizes), "; max=", max(group_sizes))

  # ------------------------------------------------------------
  # Directly place very small state groups on the final tree.
  # ------------------------------------------------------------
  # The original FateScape refinement/integration assumes non-trivial local
  # subtrees. For StateGroup ablation, some states can contain only 1--2 cells.
  # These small groups are not informative enough for local subtree refinement;
  # here they are retained as individual tips and attached directly to the root
  # of the merged tree after all sufficiently large state groups are merged.
  lineage_names_all <- paste0("L", seq_along(state_lineages))
  small_lineages <- lineage_names_all[group_sizes < MIN_STATEGROUP_TREE_CELLS]
  active_lineages <- lineage_names_all[group_sizes >= MIN_STATEGROUP_TREE_CELLS]

  direct_small_cells <- unique(unlist(cell_lineages[paste0("cell_", small_lineages)], use.names = FALSE))
  direct_small_cells <- as.character(direct_small_cells)
  direct_small_cells <- direct_small_cells[!is.na(direct_small_cells) & nzchar(direct_small_cells)]

  # Write small-group table safely. When there are no small groups,
  # all columns must have length 0; otherwise data.frame() errors with
  # "arguments imply differing number of rows: 0, 1" because a scalar
  # handling label is recycled incompatibly with zero-length columns.
  if (length(small_lineages) == 0) {
    small_group_df <- data.frame(
      lineage = character(0),
      state_group = character(0),
      n_cells = integer(0),
      handling = character(0),
      stringsAsFactors = FALSE
    )
  } else {
    small_group_df <- data.frame(
      lineage = small_lineages,
      state_group = as.character(unname(unlist(state_lineages[small_lineages]))),
      n_cells = as.integer(group_sizes[match(small_lineages, lineage_names_all)]),
      handling = rep("direct_root_tip", length(small_lineages)),
      stringsAsFactors = FALSE
    )
  }
  safe_write_csv(
    small_group_df,
    file.path(out_dir, "state_group_direct_small_groups.csv"),
    row.names = FALSE
  )
  message("State groups used for local tree construction: ", length(active_lineages),
          "; directly attached small-group cells: ", length(direct_small_cells))

  if (length(active_lineages) == 0) {
    warning("No state group has >= ", MIN_STATEGROUP_TREE_CELLS,
            " cells. Constructing a root star tree over all cells.")
    Tree_Merge <- make_star_phylo(state$cell_id)
    LiterateRefine <- list(bestsubtree = list(), bestsubtreescore = list(), total_time = NA)
    bestsubtree <- list()
    subtrees_rootbar <- NULL
    CM <- NULL
    subsubtrees <- NULL
  } else {
    # Compact active state groups to L1, L2, ... because the original FateScape
    # integration functions iterate over 1:length(state_lineages) and expect
    # these contiguous names.
    compact_names <- paste0("L", seq_along(active_lineages))
    state_lineages_active <- state_lineages[active_lineages]
    names(state_lineages_active) <- compact_names

    barcodes_lineages_active <- barcodes_lineages[active_lineages]
    names(barcodes_lineages_active) <- compact_names

    state_labels_lineages_active <- state_labels_lineages[active_lineages]
    names(state_labels_lineages_active) <- compact_names

    cell_lineages_active <- cell_lineages[paste0("cell_", active_lineages)]
    names(cell_lineages_active) <- paste0("cell_", compact_names)

    Trees_initial <- initial_tree_construction_stategroup_safe(state_lineages_active, barcodes_lineages_active)

    valid_lineages <- names(Trees_initial)
    if (length(valid_lineages) < length(state_lineages_active)) {
      state_lineages_active <- state_lineages_active[valid_lineages]
      barcodes_lineages_active <- barcodes_lineages_active[valid_lineages]
      state_labels_lineages_active <- state_labels_lineages_active[valid_lineages]
      cell_lineages_active <- cell_lineages_active[paste0("cell_", valid_lineages)]
    }

    LiterateRefine <- subtree_refinement_stategroup_safe(
      Trees_initial,
      state_lineages_active,
      barcodes_lineages_active,
      N_char,
      state_labels_lineages_active,
      lambda1 = LAMBDA_STATE,
      lambda2 = LAMBDA_BARCODE,
      maxIter = MAX_ITER,
      repeat_time = REPEAT_TIME
    )
    bestsubtree <- LiterateRefine[[1]]

    # State groups are a disjoint partition of cells by construction. Therefore
    # duplicate-tip resolution is not required here; skipping it also avoids
    # sending artificial small groups into duplicate-tip pruning.
    all_local_tips <- unlist(lapply(bestsubtree, function(tr) as.character(tr$tip.label)), use.names = FALSE)
    if (anyDuplicated(all_local_tips)) {
      warning("Duplicated tips detected unexpectedly in StateGroup active subtrees; applying drop_duplicated_tips().")
      bestsubtree <- drop_duplicated_tips(
        bestsubtree,
        barcodes_lineages_active,
        cell_lineages_active,
        state_lineage,
        alpha = DROP_DUP_ALPHA,
        beta = DROP_DUP_BETA
      )
    }

    if (length(bestsubtree) == 1) {
      Tree_Merge <- bestsubtree[[1]]
      subtrees_rootbar <- NULL
      CM <- NULL
      subsubtrees <- bestsubtree
    } else {
      subtrees_rootbar <- get_subtree_root_barcodes(bestsubtree, state, barcodes, length(bestsubtree))
      CM <- common_mutation_matrix(subtrees_rootbar)
      rank_res <- subtrees_rank(CM)
      Nodes_rank <- rank_res[[1]]
      Nodes_weight <- rank_res[[2]]
      subsubtrees <- decompose_subtrees(bestsubtree, state, barcodes, length(bestsubtree))

      if (exists("merge_subcell_trees_ordered", mode = "function")) {
        Tree_Merge <- merge_subcell_trees_ordered(
          CM,
          Nodes_rank,
          Nodes_weight,
          subsubtrees,
          bind_tree_list = list(),
          nsubtree = 1
        )
      } else if (exists("merge_subcell_trees_ward", mode = "function")) {
        Tree_Merge <- merge_subcell_trees_ward(subsubtrees)
      } else {
        stop("No supported subtree merging function found: expected merge_subcell_trees_ordered or merge_subcell_trees_ward.")
      }
    }

    # Attach small state-group cells directly as individual tips to the final root.
    Tree_Merge <- bind_direct_tips_to_root(Tree_Merge, direct_small_cells)
  }

  runtime <- proc.time() - start_time

  # Save reconstructed tree immediately after inference, before evaluation.
  reconstructed_tree_file <- save_reconstructed_tree_first(
    Tree_Merge,
    out_dir,
    file_name = "FateScape_StateGroup_reconstructed_tree.nwk"
  )

  metrics <- tryCatch(
    compute_tree_metrics(inputs$true_tree, Tree_Merge),
    error = function(e) {
      warning("Evaluation failed but reconstructed tree was saved: ", conditionMessage(e))
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
  metrics$method <- "FateScape-StateGroup"
  metrics$grouping_strategy <- "state_cluster_partition"
  metrics$barcode_source <- barcode_source
  metrics$runtime_sec <- as.numeric(runtime[["elapsed"]])
  metrics$lambda_barcode <- LAMBDA_BARCODE
  metrics$lambda_state <- LAMBDA_STATE
  metrics$max_iter <- MAX_ITER
  metrics$repeat_time <- REPEAT_TIME
  metrics$reconstructed_tree_file <- reconstructed_tree_file

  safe_write_csv(metrics, file.path(out_dir, "FateScape_StateGroup_metrics.csv"), row.names = FALSE)

  run_params <- data.frame(
    case_id = case_id,
    case_dir = case_dir,
    output_dir = out_dir,
    barcode_source = barcode_source,
    grouping_strategy = "state_cluster_partition",
    n_state_groups = length(state_lineages),
    min_stategroup_tree_cells = MIN_STATEGROUP_TREE_CELLS,
    n_direct_small_group_cells = length(direct_small_cells),
    ncells = ncells,
    N_char = N_char,
    lambda_barcode = LAMBDA_BARCODE,
    lambda_state = LAMBDA_STATE,
    max_iter = MAX_ITER,
    repeat_time = REPEAT_TIME,
    drop_dup_alpha = DROP_DUP_ALPHA,
    drop_dup_beta = DROP_DUP_BETA,
    runtime_sec = as.numeric(runtime[["elapsed"]]),
    reconstructed_tree_file = reconstructed_tree_file,
    stringsAsFactors = FALSE
  )
  safe_write_csv(run_params, file.path(out_dir, "FateScape_StateGroup_run_params.csv"), row.names = FALSE)

  if (SAVE_INTERMEDIATE_OBJECTS) {
    saveRDS(
      list(
        inputs = inputs,
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
        direct_small_cells = direct_small_cells,
        metrics = metrics
      ),
      file = file.path(out_dir, "FateScape_StateGroup_intermediate_objects.rds")
    )
  }

  list(status = "success", metrics = metrics, runtime_sec = as.numeric(runtime[["elapsed"]]))
}

# -------------------------
# 4. BATCH DRIVER
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

safe_write_csv(batch_plan, file.path(SUMMARY_DIR, "FateScape_StateGroup_batch_plan.csv"), row.names = FALSE)

message("Batch plan written to: ", file.path(SUMMARY_DIR, "FateScape_StateGroup_batch_plan.csv"))
message("Number of selected cases: ", nrow(batch_plan))

batch_records <- list()
metrics_records <- list()
error_records <- list()

for (i in seq_len(nrow(batch_plan))) {
  case_id <- batch_plan$case_id[i]
  case_dir <- batch_plan$case_dir[i]
  out_dir <- batch_plan$output_dir[i]
  metrics_file <- file.path(out_dir, "FateScape_StateGroup_metrics.csv")

  if (SKIP_FINISHED && !FORCE_RERUN && file.exists(metrics_file)) {
    message("\n[Skip] Existing metrics found for case: ", case_id)
    m <- read.csv(metrics_file, stringsAsFactors = FALSE)
    metrics_records[[length(metrics_records) + 1L]] <- m
    batch_records[[length(batch_records) + 1L]] <- data.frame(
      case_id = case_id,
      case_dir = case_dir,
      output_dir = out_dir,
      status = "skipped_existing",
      error = NA_character_,
      stringsAsFactors = FALSE
    )
    next
  }

  result <- tryCatch(
    run_fatescape_one_case(case_dir, out_dir),
    error = function(e) {
      msg <- conditionMessage(e)
      message("[ERROR] ", case_id, ": ", msg)
      if (SAVE_FAILED_INTERMEDIATE_OBJECTS) {
        safe_dir_create(out_dir)
        saveRDS(list(error = msg, case_dir = case_dir), file = file.path(out_dir, "FateScape_StateGroup_failed_case.rds"))
      }
      if (STOP_ON_ERROR) stop(e)
      list(status = "failed", error = msg)
    }
  )

  if (identical(result$status, "success")) {
    metrics_records[[length(metrics_records) + 1L]] <- result$metrics
    batch_records[[length(batch_records) + 1L]] <- data.frame(
      case_id = case_id,
      case_dir = case_dir,
      output_dir = out_dir,
      status = "success",
      error = NA_character_,
      stringsAsFactors = FALSE
    )
  } else {
    err <- result$error
    error_records[[length(error_records) + 1L]] <- data.frame(
      case_id = case_id,
      case_dir = case_dir,
      output_dir = out_dir,
      error = err,
      stringsAsFactors = FALSE
    )
    batch_records[[length(batch_records) + 1L]] <- data.frame(
      case_id = case_id,
      case_dir = case_dir,
      output_dir = out_dir,
      status = "failed",
      error = err,
      stringsAsFactors = FALSE
    )
  }

  # Incremental writes so that partial runs are preserved if interrupted.
  batch_manifest_tmp <- do.call(rbind, batch_records)
  safe_write_csv(batch_manifest_tmp, file.path(SUMMARY_DIR, "FateScape_StateGroup_batch_manifest.csv"), row.names = FALSE)

  if (length(metrics_records) > 0) {
    metrics_all_tmp <- do.call(rbind, metrics_records)
    safe_write_csv(metrics_all_tmp, file.path(SUMMARY_DIR, "FateScape_StateGroup_batch_metrics_all.csv"), row.names = FALSE)
  }

  if (length(error_records) > 0) {
    errors_tmp <- do.call(rbind, error_records)
    safe_write_csv(errors_tmp, file.path(SUMMARY_DIR, "FateScape_StateGroup_batch_errors.csv"), row.names = FALSE)
  }
}

batch_manifest <- do.call(rbind, batch_records)
safe_write_csv(batch_manifest, file.path(SUMMARY_DIR, "FateScape_StateGroup_batch_manifest.csv"), row.names = FALSE)

if (length(metrics_records) > 0) {
  metrics_all <- do.call(rbind, metrics_records)
  safe_write_csv(metrics_all, file.path(SUMMARY_DIR, "FateScape_StateGroup_batch_metrics_all.csv"), row.names = FALSE)
  message("Batch metrics written to: ", file.path(SUMMARY_DIR, "FateScape_StateGroup_batch_metrics_all.csv"))
}

if (length(error_records) > 0) {
  errors_all <- do.call(rbind, error_records)
  safe_write_csv(errors_all, file.path(SUMMARY_DIR, "FateScape_StateGroup_batch_errors.csv"), row.names = FALSE)
  message("Batch errors written to: ", file.path(SUMMARY_DIR, "FateScape_StateGroup_batch_errors.csv"))
}

message("\nBatch completed.")
print(table(batch_manifest$status))
message("Summary directory: ", SUMMARY_DIR)
