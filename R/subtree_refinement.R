#' Missing-safe barcode difference count
#'
#' Counts barcode differences without propagating NA values. Two missing values are
#' treated as uninformative and not different; one missing and one observed value
#' is counted as a difference; two observed unequal values are counted as a
#' difference.
#'
#' @param x,y Barcode vectors of the same length.
#' @return Integer number of differing sites.
barcode_diff_count <- function(x, y) {
  x <- as.vector(x)
  y <- as.vector(y)
  if (length(x) != length(y)) {
    stop("barcode_diff_count(): x and y must have the same length.")
  }

  x_na <- is.na(x)
  y_na <- is.na(y)
  one_na <- xor(x_na, y_na)
  both_obs <- !x_na & !y_na
  sum(one_na | (both_obs & x != y), na.rm = TRUE)
}

#' Calculate the probability of each node being selected as the root for subtree exchange
#'
#' @param tree Sub-cell division tree.
#' @param edges Edge matrix of the sub-cell division tree.
#' @param cell_state_labels Vector of cell state labels.
#' @param barcodes Matrix of lineage barcodes.
#' @param state_lineages List of state lineage paths.
#' @param state_threshold Positive state shifts smaller than this threshold are rewarded.
#' @param barcode_threshold Barcode differences smaller than this threshold are rewarded.
#' @param state_weight Weight for state-shift based proposal adjustment.
#' @param barcode_weight Weight for barcode-difference based proposal adjustment.
#' @param undefined_state_penalty Penalty when parent state is undefined.
#'
#' @return A probability vector for each node.
nodes_sampling <- function(tree, edges, cell_state_labels, barcodes, state_lineages,
                           state_threshold = 4, barcode_threshold = 15,
                           state_weight = 5, barcode_weight = 5,
                           undefined_state_penalty = 10) {
  total_nodes <- 2 * length(tree$tip.label) - 1
  renew_prob <- rep(0, total_nodes)
  parent_nodes <- unique(edges[, 1])

  for (par in parent_nodes) {
    children_nodes <- edges[edges[, 1] == par, 2]
    if (length(children_nodes) < 2) next

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
          # Fixed bug: use state_dist_2 for the second child.
          if (!is.na(state_dist_2) && state_dist_2 > 0 && state_dist_2 < state_threshold) {
            renew_prob[children_nodes[2]] <- renew_prob[children_nodes[2]] - (state_weight / state_dist_2)
          } else {
            renew_prob[children_nodes[2]] <- renew_prob[children_nodes[2]] + state_weight
          }
        }
      }
    }

    diff_sites <- barcode_diff_count(child_barcodes_1, child_barcodes_2)
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
    # Fallback: sample only non-root nodes if all proposal weights vanish.
    renew_prob <- rep(0, total_nodes)
    non_root <- setdiff(seq_len(total_nodes), unique(edges[, 1])[1])
    renew_prob[non_root] <- 1 / length(non_root)
  } else {
    renew_prob <- renew_prob / total
  }

  return(renew_prob)
}

#' Swap subtrees in a sub-cell division tree
#'
#' @param tree Sub-cell division tree.
#' @param edges Edge matrix of the tree.
#' @param cell_state_labels Vector of cell state labels.
#' @param barcodes Matrix of lineage barcodes.
#' @param state_lineages List of state lineage paths.
#' @param tree_renew_prob Optional probability vector for selecting nodes for subtree exchange.
#'
#' @return The tree after swapping selected subtrees.
subtrees_swapping <- function(tree, edges, cell_state_labels, barcodes, state_lineages, tree_renew_prob = NULL) {
  tree <- Preorder(tree)
  tree$node.label <- paste0("node", 1:tree$Nnode)

  if (is.null(tree_renew_prob)) {
    tree_renew_prob <- nodes_sampling(tree, edges, cell_state_labels, barcodes, state_lineages)
  }

  total_nodes <- tree$Nnode + length(tree$tip.label)
  subtree_nodes <- sample(seq_len(total_nodes), 2, prob = tree_renew_prob)
  subtree1 <- Subtree(tree, subtree_nodes[1])
  subtree2 <- Subtree(tree, subtree_nodes[2])

  guard <- 0L
  while (length(intersect(subtree1$tip.label, subtree2$tip.label)) > 0) {
    guard <- guard + 1L
    if (guard > 100L) {
      stop("subtrees_swapping(): failed to sample two non-overlapping subtrees after 100 attempts.")
    }
    subtree_nodes <- sample(seq_len(total_nodes), 2, prob = tree_renew_prob)
    subtree1 <- Subtree(tree, subtree_nodes[1])
    subtree2 <- Subtree(tree, subtree_nodes[2])
  }

  subtree1$name <- subtree_nodes[1]
  subtree1$edge.length <- rep(1, nrow(subtree1$edge))
  subtree2$name <- subtree_nodes[2]
  subtree2$edge.length <- rep(1, nrow(subtree2$edge))

  root1 <- tree$edge[tree$edge[, 2] == subtree1$name, 1]
  root2 <- tree$edge[tree$edge[, 2] == subtree2$name, 1]

  if (root1 > root2) {
    temp <- subtree2
    subtree2 <- subtree1
    subtree1 <- temp
  }

  root1 <- tree$edge[tree$edge[, 2] == subtree1$name, 1]
  root2 <- tree$edge[tree$edge[, 2] == subtree2$name, 1]
  root1_label <- paste0("node", root1 - Ntip(tree))
  root2_label <- paste0("node", root2 - Ntip(tree))

  subtree1$root.edge <- 1
  subtree2$root.edge <- 1

  if (length(subtree1$tip.label) > 1) {
    tree_prune1 <- drop.tip(tree, subtree1$tip.label, trim.internal = FALSE, subtree = FALSE, root.edge = 1)
  } else {
    tree_prune1 <- tree
    tree_prune1$tip.label[tree_prune1$tip.label == subtree1$tip.label] <- "NA"
  }
  if (length(subtree2$tip.label) > 1) {
    tree_prune2 <- drop.tip(tree_prune1, subtree2$tip.label, trim.internal = FALSE, subtree = FALSE, root.edge = 1)
  } else {
    tree_prune2 <- tree_prune1
    tree_prune2$tip.label[tree_prune2$tip.label == subtree2$tip.label] <- "NA"
  }

  where1 <- Ntip(tree_prune2) + match(root1_label, tree_prune2$node.label)
  where2_label <- root2_label
  if (is.na(where1)) stop("subtrees_swapping(): failed to locate first graft position.")
  tree_new1 <- bind.tree(tree_prune2, subtree2, where = where1)

  where2 <- Ntip(tree_new1) + match(where2_label, tree_new1$node.label)
  if (is.na(where2)) stop("subtrees_swapping(): failed to locate second graft position.")
  tree_new2 <- bind.tree(tree_new1, subtree1, where = where2)

  tree_final <- drop.tip(tree_new2, tree_new2$tip.label[grepl("node", tree_new2$tip.label, fixed = TRUE)], trim.internal = TRUE)
  if ("NA" %in% tree_final$tip.label) {
    tree_final <- drop.tip(tree_final, "NA")
  }
  tree_final$node.label <- NULL
  tree_final$edge.length <- rep(1, length(tree_final$edge.length))

  return(tree_final)
}

#' Refine sub-cell division trees based on lineage barcodes and cell states
#'
#' @param Trees_initial List of initial sub-cell division trees (by state lineage).
#' @param state_lineages List of state lineage paths.
#' @param barcodes_lineages List of barcode matrices per state lineage.
#' @param N_char Total number of sites.
#' @param state_labels_lineages List of state labels per state lineage.
#' @param lambda1 Weight for cell state score.
#' @param lambda2 Weight for barcode score.
#' @param maxIter Maximum number of iterations.
#' @param repeat_time Number of iterations required for convergence.
#'
#' @return A list containing bestsubtree, bestsubtreescore and total_time.
#' @export
subtree_refinement <- function(Trees_initial, state_lineages, barcodes_lineages, N_char,
                               state_labels_lineages, lambda1, lambda2,
                               maxIter = 100, repeat_time = 10) {
  ptm_FateScape <- proc.time()
  bestsubtreescore <- list()
  bestsubtree <- list()

  for (j in seq_along(state_lineages)) {
    lineage_label <- paste0("L", j)
    tree <- Trees_initial[[lineage_label]]
    tree$edge.length <- rep(1, nrow(tree$edge))

    score <- composite_score(
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
    barcodes <- ances_res[[1]]
    cell_state_labels <- ances_res[[2]]

    score_all <- c()
    iter_repeat <- repeat_time
    best_tree <- tree
    best_tree_list_FateScape <- list()
    best_score_list_FateScape <- c()
    tree_list <- list()

    for (i in seq_len(maxIter)) {
      renew_prob <- nodes_sampling(tree, tree$edge, cell_state_labels, barcodes, state_lineages)
      tree_new <- subtrees_swapping(tree, tree$edge, cell_state_labels, barcodes, state_lineages,
                                    tree_renew_prob = renew_prob)

      score <- composite_score(
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
      barcodes <- ances_res[[1]]
      cell_state_labels <- ances_res[[2]]

      score_all <- c(score_all, maxscore)

      if (score > maxscore) {
        maxscore <- score
        best_tree <- tree_new
        tree <- tree_new
        tree_list[[length(tree_list) + 1]] <- tree
        ances_res <- ancestor_inference(
          tree,
          N_char,
          barcodes_lineages[[lineage_label]],
          state_labels_lineages[[lineage_label]],
          state_lineages
        )
        barcodes <- ances_res[[1]]
        cell_state_labels <- ances_res[[2]]
      }

      if (i %% iter_repeat == 0) {
        recent_scores <- score_all[(i - iter_repeat + 1):i]
        if (length(unique(recent_scores)) == 1) {
          best_tree_list_FateScape[[length(best_tree_list_FateScape) + 1]] <- best_tree
          best_score_list_FateScape <- c(best_score_list_FateScape, maxscore)

          tree <- Trees_initial[[lineage_label]]
          tree$edge.length <- rep(1, nrow(tree$edge))
          score <- composite_score(
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
          barcodes <- ances_res[[1]]
          cell_state_labels <- ances_res[[2]]
        }
      }
    }

    # Fallback if convergence checkpoint was not triggered during maxIter.
    if (length(best_score_list_FateScape) == 0) {
      best_score_list_FateScape <- maxscore
      best_tree_list_FateScape <- list(best_tree)
    }

    bestsubtreescore[[lineage_label]] <- max(best_score_list_FateScape)
    best_index <- match(max(best_score_list_FateScape), best_score_list_FateScape)
    bestsubtree[[lineage_label]] <- best_tree_list_FateScape[[best_index]]
  }

  total_time <- proc.time() - ptm_FateScape
  return(list(bestsubtree = bestsubtree, bestsubtreescore = bestsubtreescore, total_time = total_time))
}
