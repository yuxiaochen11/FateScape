# Robust missing-safe subtree integration utilities for FateScape
# -----------------------------------------------------------------------------
# This file replaces subtree_integration.R. It is designed to be robust to:
#   1. retained NA / '-' / '?' barcode entries;
#   2. pruned trees with unary or missing children;
#   3. empty candidate split sets;
#   4. all-NA / all-zero common-mutation matrices;
#   5. edge cases in ape::drop.tip(), ape::bind.tree(), and phytools::Subtree().
# -----------------------------------------------------------------------------

# ---- Missing-safe barcode helpers ------------------------------------------------
missing_barcode_tokens <- c("", "NA", "NaN", "-", "?", "missing", "Missing", "MISSING", "NULL", "null")

is_missing_barcode_value <- function(x) {
  x_chr <- as.character(x)
  is.na(x_chr) | x_chr %in% missing_barcode_tokens
}

barcode_to_numeric_zero_for_missing <- function(x) {
  x_chr <- as.character(x)
  x_chr[is_missing_barcode_value(x_chr)] <- "0"
  suppressWarnings(as.integer(x_chr))
}

barcode_to_numeric_na_for_missing <- function(x) {
  x_chr <- as.character(x)
  x_chr[is_missing_barcode_value(x_chr)] <- NA_character_
  suppressWarnings(as.integer(x_chr))
}

all_zero_barcode <- function(x) {
  x_num <- barcode_to_numeric_zero_for_missing(x)
  length(x_num) > 0 && all(!is.na(x_num) & x_num == 0)
}

barcode_diff_count <- function(a, b) {
  a <- barcode_to_numeric_na_for_missing(a)
  b <- barcode_to_numeric_na_for_missing(b)
  n <- min(length(a), length(b))
  if (n == 0) return(0L)
  a <- a[seq_len(n)]
  b <- b[seq_len(n)]
  both_na <- is.na(a) & is.na(b)
  one_na <- xor(is.na(a), is.na(b))
  both_obs <- !is.na(a) & !is.na(b)
  sum(one_na | (both_obs & a != b), na.rm = TRUE)
}

safe_improved_hamming_distance <- function(barcode1, barcode2, alpha = 0, beta = 1) {
  a <- barcode_to_numeric_na_for_missing(barcode1)
  b <- barcode_to_numeric_na_for_missing(barcode2)
  n <- min(length(a), length(b))
  if (n == 0) return(0)
  a <- a[seq_len(n)]
  b <- b[seq_len(n)]

  dist <- 0
  for (i in seq_len(n)) {
    ai <- a[i]
    bi <- b[i]

    # Two missing sites are uninformative; one missing site is penalized.
    if (is.na(ai) && is.na(bi)) {
      next
    } else if (is.na(ai) || is.na(bi)) {
      dist <- dist + beta
      next
    }

    if (ai == bi && ai != 0) {
      dist <- dist - alpha
    } else if ((ai == 0 && bi != 0) || (ai != 0 && bi == 0)) {
      dist <- dist + 1
    } else if (ai != bi && ai != 0 && bi != 0) {
      dist <- dist + beta
    }
  }
  if (is.na(dist) || !is.finite(dist)) dist <- Inf
  dist
}

# Prefer the global improved_hamming_distance if it has already been replaced by
# a missing-safe implementation; otherwise fall back to safe version locally.
.safe_hamming <- function(a, b, alpha = 0, beta = 1) {
  out <- tryCatch(improved_hamming_distance(a, b, alpha, beta), error = function(e) NA_real_)
  if (length(out) != 1 || is.na(out) || !is.finite(out)) {
    out <- safe_improved_hamming_distance(a, b, alpha, beta)
  }
  out
}

# ---- Tree/node helpers -----------------------------------------------------------
get_first_element <- function(x) {
  if (length(x) < 2) return(NA_character_)
  x[[2]]
}

get_tips <- function(x) {
  if (is.null(x) || is.null(x$tip.label)) return(character(0))
  unlist(x$tip.label)
}

get_root_barcode <- function(x) {
  if (is.null(x) || is.null(x$root.barcode)) return(matrix(0, nrow = 1, ncol = if (exists("N_char")) N_char else 0))
  x$root.barcode
}

subtree_root_bar <- function(x) {
  if (length(x) == 0) return(list())
  lapply(x, get_root_barcode)
}

safe_subtree <- function(tree, node) {
  if (is.null(tree) || is.na(node) || length(node) != 1) return(NULL)
  out <- tryCatch(Subtree(tree, node), error = function(e) NULL)
  if (is.null(out) || is.null(out$tip.label) || length(out$tip.label) == 0) return(NULL)
  out
}

safe_drop_tip <- function(tree, tips, ...) {
  if (is.null(tree)) return(NULL)
  tips <- unique(tips[!is.na(tips)])
  tips <- intersect(tips, tree$tip.label)
  if (length(tips) == 0) return(tree)
  if (length(tips) >= length(tree$tip.label)) return(NULL)
  tryCatch(drop.tip(tree, tips, ...), error = function(e) NULL)
}

safe_bind_root <- function(base_tree, add_tree) {
  if (is.null(base_tree)) return(add_tree)
  if (is.null(add_tree)) return(base_tree)
  out <- tryCatch(bind.tree(base_tree, add_tree, where = "root", position = 0.1), error = function(e) NULL)
  if (is.null(out)) {
    out <- tryCatch(bind.tree(base_tree, add_tree, where = "root"), error = function(e) base_tree)
  }
  out$root.edge <- 1
  out
}

node_name_from_index <- function(idx, tip_labels, barcode_names = NULL) {
  if (is.na(idx) || length(idx) != 1) return(NA_character_)
  ntips <- length(tip_labels)
  candidates <- character(0)
  if (idx <= ntips && idx >= 1) {
    candidates <- c(candidates, tip_labels[idx])
  } else {
    candidates <- c(candidates, paste("node", idx, sep = "_"), paste("node", idx, sep = ""), as.character(idx))
  }
  if (!is.null(barcode_names)) {
    hit <- candidates[candidates %in% barcode_names]
    if (length(hit) > 0) return(hit[1])
  }
  candidates[1]
}

barcode_row_by_node <- function(barances, idx, tip_labels) {
  rn <- rownames(barances)
  nm <- node_name_from_index(idx, tip_labels, rn)
  if (!is.na(nm) && nm %in% rn) return(barances[nm, , drop = FALSE])
  if (!is.na(idx) && idx >= 1 && idx <= nrow(barances)) return(barances[idx, , drop = FALSE])
  matrix(0, nrow = 1, ncol = ncol(barances))
}

# ---- Duplicate-tip resolution ----------------------------------------------------
sibling_cells_distance <- function(cell, tree, lineage, alpha, beta) {
  if (is.null(tree) || !(cell %in% tree$tip.label)) return(Inf)
  tips <- tree$tip.label
  edges <- tree$edge
  cell_index <- which(tips == cell)
  if (length(cell_index) != 1) return(Inf)

  parent_node <- edges[edges[, 2] == cell_index, 1]
  if (length(parent_node) < 1 || is.na(parent_node[1])) return(Inf)
  parent_node <- parent_node[1]

  siblings <- edges[edges[, 1] == parent_node, 2]
  siblings <- siblings[!is.na(siblings)]
  if (length(siblings) < 2) return(Inf)

  ances_res <- tryCatch(
    ancestor_inference(tree, N_char, barcodes_lineages[[lineage]],
                       state_labels_lineages[[lineage]], state_lineages),
    error = function(e) NULL
  )
  if (is.null(ances_res) || length(ances_res) < 1) return(Inf)

  bar_all <- ances_res[[1]]
  valid_sibs <- siblings[siblings >= 1 & siblings <= nrow(bar_all)]
  if (length(valid_sibs) < 2) return(Inf)
  sibling_barcodes <- bar_all[valid_sibs[1:2], , drop = FALSE]

  dist <- .safe_hamming(sibling_barcodes[1, ], sibling_barcodes[2, ], alpha, beta)
  if (is.na(dist) || !is.finite(dist)) Inf else dist
}

drop_non_minimum_tips <- function(tree_list, dist_list, cell) {
  dist <- suppressWarnings(as.numeric(dist_list$dist))
  dist[is.na(dist) | !is.finite(dist)] <- Inf
  if (length(dist) == 0) return(tree_list)

  if (all(is.infinite(dist))) {
    min_lineage <- dist_list$lineage[1]
  } else {
    min_lineage <- dist_list$lineage[which.min(dist)]
  }
  if (length(min_lineage) != 1 || is.na(min_lineage)) min_lineage <- dist_list$lineage[1]

  drop_lineages <- dist_list$lineage[dist_list$lineage != min_lineage]
  for (j in drop_lineages) {
    if (!is.null(tree_list[[j]]) && cell %in% tree_list[[j]]$tip.label) {
      new_tree <- safe_drop_tip(tree_list[[j]], cell)
      if (!is.null(new_tree)) tree_list[[j]] <- new_tree
    }
  }
  tree_list
}

drop_duplicated_tips <- function(tree_list, barcodes, cell_lineages, state_lineage, alpha, beta) {
  if (length(state_lineage) == 0) return(tree_list)
  for (i in seq_along(state_lineage)) {
    if (length(state_lineage[[i]]) > 1) {
      cell <- names(state_lineage)[i]
      if (is.null(cell) || is.na(cell) || !nzchar(cell)) next
      lineages <- state_lineage[[i]]
      lineages <- lineages[lineages %in% names(tree_list)]
      if (length(lineages) <= 1) next
      dist_vec <- sapply(lineages, function(lin) {
        sibling_cells_distance(cell, tree_list[[lin]], lineage = lin, alpha, beta)
      })
      dist_list <- data.frame(lineage = lineages, dist = dist_vec)
      tree_list <- drop_non_minimum_tips(tree_list, dist_list, cell)
    }
  }
  tree_list
}

# ---- Subtree decomposition -------------------------------------------------------
prune_cell_tree <- function(sub_tree, state, barcodes, subid, itr_time = 1) {
  if (is.null(sub_tree) || length(sub_tree$tip.label) == 0) return(list())
  sub_tree <- tryCatch(Preorder(sub_tree), error = function(e) sub_tree)
  ntips <- length(sub_tree$tip.label)

  st_barleaves <- barcodes[sub_tree$tip.label, , drop = FALSE]
  st_stateleaves <- state$cluster[match(sub_tree$tip.label, state$cell_id)]

  # If the subtree has only one tip, return it directly.
  if (ntips <= 1 || is.null(sub_tree$edge) || nrow(sub_tree$edge) == 0) {
    sub_tree$root.barcode <- st_barleaves[1, , drop = FALSE]
    sub_tree$root.edge <- 1
    nm <- paste(itr_time, subid, "single", sub_tree$tip.label[1], sep = "_")
    out <- list(sub_tree)
    names(out) <- nm
    return(out)
  }

  ances <- tryCatch(
    ancestor_inference(sub_tree, N_char, st_barleaves, st_stateleaves, state_lineages)[[1]],
    error = function(e) NULL
  )
  if (is.null(ances) || nrow(ances) < ntips) {
    sub_tree$root.barcode <- st_barleaves[1, , drop = FALSE]
    sub_tree$root.edge <- 1
    nm <- paste(itr_time, subid, "whole", sep = "_")
    out <- list(sub_tree)
    names(out) <- nm
    return(out)
  }

  internal_idx <- (ntips + 1):min(nrow(ances), 2 * ntips - 1)
  st_barances <- ances[internal_idx, , drop = FALSE]

  bar_ances <- matrix(barcode_to_numeric_zero_for_missing(st_barances), nrow = nrow(st_barances))
  rownames(bar_ances) <- rownames(st_barances)
  bar_leaves <- matrix(barcode_to_numeric_zero_for_missing(st_barleaves), nrow = nrow(st_barleaves))
  rownames(bar_leaves) <- rownames(st_barleaves)
  barances <- rbind(bar_leaves, bar_ances)

  candidate_split_nodes_name <- rownames(barances)[which(rowSums(barances, na.rm = TRUE) == 0)]
  # Internal nodes only; tips with all-zero barcodes should not be used as split nodes.
  candidate_split_nodes_name <- setdiff(candidate_split_nodes_name, sub_tree$tip.label)
  id <- suppressWarnings(as.integer(vapply(strsplit(candidate_split_nodes_name, "_"), get_first_element, character(1))))
  id <- unique(id[!is.na(id) & id > ntips])

  subsubtrees_list <- list()
  subsubtree_left <- NULL

  if (length(id) == 0) {
    root_idx <- setdiff(sub_tree$edge[, 1], sub_tree$edge[, 2])
    root_idx <- if (length(root_idx) > 0) root_idx[1] else ntips + 1
    sub_tree$root.barcode <- barcode_row_by_node(barances, root_idx, sub_tree$tip.label)
    sub_tree$root.edge <- 1
    key <- paste(itr_time, subid, "whole", sep = "_")
    subsubtrees_list[[key]] <- sub_tree
    return(subsubtrees_list)
  }

  for (split_node in id) {
    children_index <- sub_tree$edge[sub_tree$edge[, 1] == split_node, 2]
    children_index <- children_index[!is.na(children_index)]
    if (length(children_index) < 2) next
    children_index <- children_index[seq_len(2)]

    # If children are themselves split nodes, skip to avoid nested duplicates.
    if (any(children_index %in% id)) next

    for (child_idx in children_index) {
      subtree_child <- safe_subtree(sub_tree, child_idx)
      if (is.null(subtree_child)) next
      subtree_child$root.barcode <- barcode_row_by_node(barances, child_idx, sub_tree$tip.label)
      subtree_child$root.edge <- 1
      key <- paste(itr_time, subid, "node", child_idx, sep = "_")
      subsubtrees_list[[key]] <- subtree_child
    }
  }

  if (length(subsubtrees_list) == 0) {
    root_idx <- setdiff(sub_tree$edge[, 1], sub_tree$edge[, 2])
    root_idx <- if (length(root_idx) > 0) root_idx[1] else ntips + 1
    sub_tree$root.barcode <- barcode_row_by_node(barances, root_idx, sub_tree$tip.label)
    sub_tree$root.edge <- 1
    key <- paste(itr_time, subid, "whole", sep = "_")
    subsubtrees_list[[key]] <- sub_tree
    return(subsubtrees_list)
  }

  if (length(id) != 1) {
    drop_tips_vec <- unique(unlist(lapply(subsubtrees_list, get_tips)))
    subsubtree_left <- safe_drop_tip(sub_tree, drop_tips_vec, trim.internal = TRUE, subtree = FALSE, root.edge = 1)
    if (!is.null(subsubtree_left) && length(subsubtree_left$tip.label) > 0) {
      barleft <- barcodes[subsubtree_left$tip.label, , drop = FALSE]
      stateleft <- state$cluster[match(subsubtree_left$tip.label, state$cell_id)]
      left_ances <- tryCatch(
        ancestor_inference(subsubtree_left, N_char, barleft, stateleft, state_lineages)[[1]],
        error = function(e) NULL
      )
      if (!is.null(left_ances) && nrow(left_ances) > length(subsubtree_left$tip.label)) {
        left_internal <- (length(subsubtree_left$tip.label) + 1):nrow(left_ances)
        left_barances <- left_ances[left_internal, , drop = FALSE]
        leftbarances <- matrix(barcode_to_numeric_zero_for_missing(left_barances), nrow = nrow(left_barances))
        rownames(leftbarances) <- rownames(left_barances)
        subsubtree_left$root.barcode <- leftbarances[1, , drop = FALSE]
      } else {
        subsubtree_left$root.barcode <- barleft[1, , drop = FALSE]
      }
      subsubtree_left$root.edge <- 1
      key_left <- paste(itr_time, subid, "left", sep = "_")
      subsubtrees_list[[key_left]] <- subsubtree_left
    }
  }

  # Optional recursive refinement only if the last subtree has a zero root barcode.
  roots <- subtree_root_bar(subsubtrees_list)
  if (length(roots) == 0) return(subsubtrees_list)
  last_root <- roots[[length(roots)]]
  if (!all_zero_barcode(last_root) || is.null(subsubtree_left) || itr_time >= 5) {
    return(subsubtrees_list)
  }
  c(subsubtrees_list[-length(subsubtrees_list)], prune_cell_tree(subsubtree_left, state, barcodes, subid, itr_time + 1))
}

# ---- Common mutation and ranking -------------------------------------------------
count_common_mutations <- function(barcode1, barcode2, N_char = NULL) {
  a <- barcode_to_numeric_na_for_missing(barcode1)
  b <- barcode_to_numeric_na_for_missing(barcode2)
  if (is.null(N_char)) N_char <- min(length(a), length(b))
  if (N_char == 0) return(0L)
  a <- a[seq_len(N_char)]
  b <- b[seq_len(N_char)]
  sum(!is.na(a) & !is.na(b) & a == b & a != 0, na.rm = TRUE)
}

count_uncommon_mutations <- function(barcode1, barcode2, N_char = NULL) {
  a <- barcode_to_numeric_na_for_missing(barcode1)
  b <- barcode_to_numeric_na_for_missing(barcode2)
  if (is.null(N_char)) N_char <- min(length(a), length(b))
  if (N_char == 0) return(0L)
  a <- a[seq_len(N_char)]
  b <- b[seq_len(N_char)]
  sum(!is.na(a) & !is.na(b) & a != b & a != 0 & b != 0, na.rm = TRUE)
}

get_subtree_root_barcodes <- function(bestsubtree, state, barcodes, l_sl) {
  root_bar <- list()
  if (length(bestsubtree) == 0) return(root_bar)
  n_iter <- min(l_sl, length(bestsubtree))
  for (i in seq_len(n_iter)) {
    if (is.null(bestsubtree[[i]])) next
    pruned <- prune_cell_tree(bestsubtree[[i]], state, barcodes, i, 1)
    root_bar <- c(root_bar, subtree_root_bar(pruned))
  }
  root_bar
}

common_mutation_matrix <- function(subtrees_rootbar) {
  if (length(subtrees_rootbar) == 0) return(matrix(numeric(0), nrow = 0, ncol = 0))
  CM <- matrix(0, nrow = length(subtrees_rootbar), ncol = length(subtrees_rootbar))
  rownames(CM) <- colnames(CM) <- names(subtrees_rootbar)
  for (i in seq_along(subtrees_rootbar)) {
    for (j in seq_along(subtrees_rootbar)) {
      CM[i, j] <- count_common_mutations(subtrees_rootbar[[i]], subtrees_rootbar[[j]])
      if (i >= j) CM[i, j] <- NA
    }
  }
  CM
}

subtrees_rank <- function(CM) {
  Nodes <- rownames(CM)
  N <- length(Nodes)
  Nodes_rank <- matrix(NA_character_, nrow = N, ncol = max(N, 1))
  Nodes_weight <- matrix(NA_real_, nrow = N, ncol = max(N, 1))
  rownames(Nodes_rank) <- rownames(Nodes_weight) <- Nodes
  if (N == 0) return(list(Nodes_rank, Nodes_weight))

  for (i in seq_len(N)) {
    vals <- CM[i, ]
    vals[is.na(vals)] <- -Inf
    num_positive <- sum(vals > 0, na.rm = TRUE)
    if (num_positive <= 0) next
    ordered_indices <- order(vals, decreasing = TRUE)
    keep <- ordered_indices[seq_len(num_positive)]
    Nodes_rank[i, seq_len(num_positive)] <- Nodes[keep]
    Nodes_weight[i, seq_len(num_positive)] <- vals[keep]
  }
  list(Nodes_rank, Nodes_weight)
}

decompose_subtrees <- function(bestsubtree, state, barcodes, l_sl) {
  subtrees <- list()
  if (length(bestsubtree) == 0) return(subtrees)
  n_iter <- min(l_sl, length(bestsubtree))
  for (i in seq_len(n_iter)) {
    if (is.null(bestsubtree[[i]])) next
    subtrees <- c(subtrees, prune_cell_tree(bestsubtree[[i]], state, barcodes, i, 1))
  }
  subtrees
}

# ---- Merge functions -------------------------------------------------------------
merge_list_at_root <- function(tree_list) {
  tree_list <- tree_list[!vapply(tree_list, is.null, logical(1))]
  if (length(tree_list) == 0) return(NULL)
  final_tree <- tree_list[[1]]
  final_tree$root.edge <- 1
  if (length(tree_list) > 1) {
    for (i in 2:length(tree_list)) {
      final_tree <- safe_bind_root(final_tree, tree_list[[i]])
    }
  }
  final_tree
}

select_best_group <- function(Nodes_rank, Nodes_weight, Nodes) {
  W <- Nodes_weight
  W[is.na(W)] <- -Inf
  if (length(W) == 0 || all(!is.finite(W)) || max(W, na.rm = TRUE) <= 0) {
    return(Nodes[1])
  }
  idx <- which(W == max(W, na.rm = TRUE), arr.ind = TRUE)[1, , drop = FALSE]
  node1 <- rownames(W)[idx[1, 1]]
  node2 <- Nodes_rank[idx[1, 1], idx[1, 2]]
  group <- unique(na.omit(c(node1, node2, Nodes_rank[node1, ], Nodes_rank[node2, ])))
  group[group %in% Nodes]
}

merge_subcell_trees_ordered <- function(CM, Nodes_rank, Nodes_weight, subsubtrees, bind_tree_list = list(), nsubtree = 1) {
  if (length(subsubtrees) == 0) return(NULL)
  for (i in seq_along(subsubtrees)) {
    if (!is.null(subsubtrees[[i]]) && !is.null(subsubtrees[[i]]$edge)) {
      subsubtrees[[i]]$edge.length <- rep(1, nrow(subsubtrees[[i]]$edge))
    }
  }
  Nodes <- rownames(Nodes_rank)
  if (length(Nodes) == 0) return(merge_list_at_root(subsubtrees))

  Nodes_remain <- Nodes
  while (length(Nodes_remain) > 0) {
    sub_CM <- CM[Nodes_remain, Nodes_remain, drop = FALSE]
    rank_weight <- subtrees_rank(sub_CM)
    group <- select_best_group(rank_weight[[1]], rank_weight[[2]], Nodes_remain)
    if (length(group) == 0) group <- Nodes_remain[1]
    bind_tree_list[[length(bind_tree_list) + 1]] <- merge_list_at_root(subsubtrees[group])
    Nodes_remain <- setdiff(Nodes_remain, group)
  }
  merge_list_at_root(bind_tree_list)
}

group_subcell_trees <- function(CM, Nodes_rank, Nodes_weight, subsubtrees, bind_tree_list = list(), nsubtree = 1) {
  if (length(subsubtrees) == 0) return(list())
  Nodes <- rownames(Nodes_rank)
  if (length(Nodes) == 0) return(list(merge_list_at_root(subsubtrees), character(0)))
  group <- select_best_group(Nodes_rank, Nodes_weight, Nodes)
  if (length(group) == 0) group <- Nodes[1]
  bind_tree_list[[nsubtree]] <- merge_list_at_root(subsubtrees[group])
  Nodes_remain <- setdiff(Nodes, group)
  list(bind_tree_list, Nodes_remain)
}

safe_barcode_distance_matrix <- function(barcodes_mat) {
  n <- nrow(barcodes_mat)
  D <- matrix(0, n, n)
  rownames(D) <- colnames(D) <- rownames(barcodes_mat)
  if (n <= 1) return(D)
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      D[i, j] <- D[j, i] <- barcode_diff_count(barcodes_mat[i, ], barcodes_mat[j, ])
    }
  }
  D
}

merge_subcell_trees_ward <- function(subtrees_rootbar, subsubtrees) {
  if (length(subsubtrees) == 0) return(NULL)
  if (length(subsubtrees) == 1) return(subsubtrees[[1]])

  for (i in seq_along(subsubtrees)) {
    if (!is.null(subsubtrees[[i]]) && !is.null(subsubtrees[[i]]$edge)) {
      subsubtrees[[i]]$edge.length <- rep(1, nrow(subsubtrees[[i]]$edge))
    }
  }

  b <- do.call(rbind, subtrees_rootbar)
  rownames(b) <- names(subtrees_rootbar)
  D <- safe_barcode_distance_matrix(b)
  # If all distances are zero, binding all subtrees at root is deterministic and safe.
  if (all(D == 0, na.rm = TRUE)) return(merge_list_at_root(subsubtrees))

  tree_WD <- tryCatch(nodes_clustering(D), error = function(e) NULL)
  if (is.null(tree_WD)) return(merge_list_at_root(subsubtrees))

  for (i in names(subsubtrees)) {
    position <- which(i == tree_WD$tip.label)
    if (length(position) != 1) next
    tree_WD <- tryCatch(bind.tree(tree_WD, subsubtrees[[i]], position), error = function(e) tree_WD)
  }
  tree_WD
}

# ---- Downstream helper functions retained with NA guards -------------------------
get_subtree_celltypes <- function(tree_groups) {
  subtree_celltypes_list <- list()
  for (j in seq_along(tree_groups)) {
    if (j == 1) {
      for (i in seq_along(tree_groups[[j]])) {
        subtree_celltype <- unique(na.omit(name_index$cell_type[match(tree_groups[[1]][[i]]$tip.label, name_index$id)]))
        subtree_celltypes_list[[paste(j, i, sep = "_")]] <- subtree_celltype
      }
    } else {
      for (i in tree_groups[[j]]) {
        subtree_celltype <- unique(na.omit(name_index$cell_type[match(subsubtrees[[i]]$tip.label, name_index$id)]))
        subtree_celltypes_list[[paste(j, i, sep = "_")]] <- subtree_celltype
      }
    }
  }
  names(subtree_celltypes_list) <- paste("subtree", seq_along(subtree_celltypes_list), sep = "_")
  subtree_celltypes_list
}

get_unique_celltypes <- function(subtree_celltypes) {
  all_elements <- unique(unlist(subtree_celltypes))
  element_counts <- table(unlist(subtree_celltypes))
  unique_elements <- list()
  for (element in names(element_counts)) {
    if (element_counts[element] == 1) {
      vec_index <- which(sapply(subtree_celltypes, function(x) element %in% x))
      if (length(vec_index) == 1) unique_elements[[element]] <- names(subtree_celltypes)[vec_index]
    }
  }
  unique_elements
}

phylogenetic_tree <- function(tree, N_char, barcodes, ncells = NULL, Nnodes = NULL, edges = NULL) {
  if (is.null(ncells)) ncells <- length(tree$tip.label)
  if (is.null(Nnodes)) Nnodes <- 2 * ncells - 2
  edges <- data.frame(par = tree$edge[, 1], child = tree$edge[, 2], length = rep(0, nrow(tree$edge)))
  barcodes_ances <- ancestor_inference_time_scale(tree, N_char, barcodes)
  for (par in unique(edges$par)) {
    children <- edges$child[edges$par == par]
    if (length(children) >= 1) {
      for (child in children) {
        if (par <= nrow(barcodes_ances) && child <= nrow(barcodes_ances)) {
          H_D <- .safe_hamming(barcodes_ances[par, ], barcodes_ances[child, ], 0, 1)
          edges$length[edges$child == child] <- ifelse(H_D == 0 || !is.finite(H_D), 0.5, H_D)
        }
      }
    }
  }
  tree$edge.length <- edges$length
  tree
}
