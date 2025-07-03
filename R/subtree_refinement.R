#' Calculate the probability of each node being selected as the root for subtree exchange
#'
#' @param tree Sub-cell division tree.
#' @param edges Edge matrix of the sub-cell division tree.
#' @param cell_state_labels Vector of cell state labels.
#' @param barcodes Matrix of lineage barcodes.
#' @param state_lineages List of state lineage paths.
#'
#' @return A probability vector for each node.
nodes_sampling<- function(tree, edges, cell_state_labels, barcodes, state_lineages) {
  # Initialize probability vector for all nodes in the tree.
  total_nodes <- 2 * length(tree$tip.label) - 1
  renew_prob <- rep(0, total_nodes)

  # Get unique parent nodes from the edge matrix.
  parent_nodes <- unique(edges[, 1])

  for (par in parent_nodes) {
    # Get children nodes for current parent.
    children_nodes <- edges[edges[, 1] == par, 2]

    # Determine first child index.
    if (children_nodes[1] > length(tree$tip.label)) {
      child_id_1 <- children_nodes[1]
    } else {
      cell_id_1 <- tree$tip.label[children_nodes[1]]
      child_id_1 <- match(cell_id_1, rownames(barcodes))
    }

    # Determine second child index.
    if (children_nodes[2] > length(tree$tip.label)) {
      child_id_2 <- children_nodes[2]
    } else {
      cell_id_2 <- tree$tip.label[children_nodes[2]]
      child_id_2 <- match(cell_id_2, rownames(barcodes))
    }

    # Retrieve barcodes for both children.
    child_barcodes_1 <- barcodes[child_id_1, ]
    child_barcodes_2 <- barcodes[child_id_2, ]

    # Get state labels.
    parent_state  <- cell_state_labels[par]
    child_state_1 <- cell_state_labels[child_id_1]
    child_state_2 <- cell_state_labels[child_id_2]

    if (parent_state == 0) {
      # If parent state is undefined (0), add a high penalty.
      renew_prob[children_nodes[1]] <- renew_prob[children_nodes[1]] + 10
      renew_prob[children_nodes[2]] <- renew_prob[children_nodes[2]] + 10
    } else {
      # For each state lineage, adjust probabilities based on state shifts.
      for (lineage in state_lineages) {
        if ((parent_state %in% lineage) && (child_state_1 %in% lineage)) {
          state_dist_1 <- match(child_state_1, lineage) - match(parent_state, lineage)
          if (state_dist_1 > 0 && state_dist_1 < 4) {
            renew_prob[children_nodes[1]] <- renew_prob[children_nodes[1]] - (5 / state_dist_1)
          } else {
            renew_prob[children_nodes[1]] <- renew_prob[children_nodes[1]] + 5
          }
        }
      }
      for (lineage in state_lineages) {
        if ((parent_state %in% lineage) && (child_state_2 %in% lineage)) {
          state_dist_2 <- match(child_state_2, lineage) - match(parent_state, lineage)
          # NOTE: Using state_dist_1 (from first child) in condition as in the original code.
          if (state_dist_1 > 0 && state_dist_1 < 4) {
            renew_prob[children_nodes[2]] <- renew_prob[children_nodes[2]] - (5 / state_dist_2)
          } else {
            renew_prob[children_nodes[2]] <- renew_prob[children_nodes[2]] + 5
          }
        }
      }
    }

    # Adjust probabilities based on barcode differences between the two children.
    diff_sites <- sum(child_barcodes_1 != child_barcodes_2)
    if (diff_sites < 15) {
      renew_prob[children_nodes[1]] <- renew_prob[children_nodes[1]] - (5 / (diff_sites + 1))
      renew_prob[children_nodes[2]] <- renew_prob[children_nodes[2]] - (5 / (diff_sites + 1))
    } else {
      renew_prob[children_nodes[1]] <- renew_prob[children_nodes[1]] + 5
      renew_prob[children_nodes[2]] <- renew_prob[children_nodes[2]] + 5
    }
  }

  # Remove negative probabilities and normalize.
  renew_prob[renew_prob < 0] <- 0
  renew_prob <- renew_prob / sum(renew_prob)

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
  # Reorder tree in preorder and label internal nodes.
  tree <- Preorder(tree)
  tree$node.label <- paste0("node", 1:tree$Nnode)

  # Compute renewal probabilities if not provided.
  if (is.null(tree_renew_prob)) {
    tree_renew_prob <-  nodes_sampling(tree, edges, cell_state_labels, barcodes, state_lineages)
  }

  # Sample two nodes based on renewal probabilities.
  total_nodes <- tree$Nnode + length(tree$tip.label)
  subtree_nodes <- sample(total_nodes, 2, prob = tree_renew_prob)
  subtree1 <- Subtree(tree, subtree_nodes[1])
  subtree2 <- Subtree(tree, subtree_nodes[2])

  # Resample if the selected subtrees overlap.
  while (length(intersect(subtree1$tip.label, subtree2$tip.label)) > 0) {
    subtree_nodes <- sample(total_nodes, 2, prob = tree_renew_prob)
    subtree1 <- Subtree(tree, subtree_nodes[1])
    subtree2 <- Subtree(tree, subtree_nodes[2])
  }

  # Set subtree names and standardize edge lengths.
  subtree1$name <- subtree_nodes[1]
  subtree1$edge.length <- rep(1, nrow(subtree1$edge))
  subtree2$name <- subtree_nodes[2]
  subtree2$edge.length <- rep(1, nrow(subtree2$edge))

  # Identify the parent nodes (roots) of the selected subtrees.
  root1 <- tree$edge[tree$edge[, 2] == subtree1$name, 1]
  root2 <- tree$edge[tree$edge[, 2] == subtree2$name, 1]

  # Ensure consistent ordering of subtrees.
  if (root1 > root2) {
    temp <- subtree2
    subtree2 <- subtree1
    subtree1 <- temp
  }

  # Recompute parent nodes after potential swap.
  root1 <- tree$edge[tree$edge[, 2] == subtree1$name, 1]
  root2 <- tree$edge[tree$edge[, 2] == subtree2$name, 1]
  root1_label <- paste0("node", root1 - Ntip(tree))
  root2_label <- paste0("node", root2 - Ntip(tree))

  subtree1$root.edge <- 1
  subtree2$root.edge <- 1

  # Prune the tree to remove the selected subtrees.
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

  # Reconstruct the tree by binding the subtrees at the determined positions.
  tree_new1 <- bind.tree(tree_prune2, subtree2, where = Ntip(tree_prune2) + match(root1_label, tree_prune2$node.label))
  tree_new2 <- bind.tree(tree_new1, subtree1, where = Ntip(tree_new1) + match(root2_label, tree_new1$node.label))

  # Remove any temporary tip labels.
  tree_final <- drop.tip(tree_new2, tree_new2$tip.label[grepl("node", tree_new2$tip.label, fixed = TRUE)], trim.internal = TRUE)
  tree_final <- drop.tip(tree_final, "NA")
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
#' @return A list containing:
#' \describe{
#'   \item{bestsubtree}{List of refined sub-cell division trees (one per state lineage).}
#'   \item{bestsubtreescore}{List of corresponding best scores for each refined tree.}
#'   \item{total_time}{Total running time.}
#' }
#' @export
subtree_refinement <- function(Trees_initial, state_lineages, barcodes_lineages, N_char,
                           state_labels_lineages, lambda1, lambda2, maxIter = 100, repeat_time = 10) {
  ptm_LineageCast <- proc.time()
  bestsubtreescore <- list()
  bestsubtree <- list()

  # Process each state lineage individually.
  for (j in seq_along(state_lineages)) {
    lineage_label <- paste0("L", j)
    tree <- Trees_initial[[lineage_label]]
    tree$edge.length <- rep(1, nrow(tree$edge))
    edges <- tree$edge

    # Calculate the initial tree score.
    score <- combined_likelihood(tree = tree, barcodes = barcodes_lineages[[lineage_label]], N_char = N_char,
                       cell_state_labels = state_labels_lineages[[lineage_label]], state_lineages = state_lineages,
                       state_score = NULL, barcode_score = NULL, lambda_1 = lambda1, lambda_2 = lambda2)
    maxscore <- score

    # Perform ancestral inference to update barcodes and cell state labels.
    ances_res <- ancestor_inference(tree = tree, N_char = N_char, barcodes = barcodes_lineages[[lineage_label]],
                                cell_state_labels = state_labels_lineages[[lineage_label]], state_lineages = state_lineages)
    barcodes <- ances_res[[1]]
    cell_state_labels <- ances_res[[2]]

    score_all <- c()
    iter_repeat <- repeat_time
    best_tree <- tree
    best_tree_list_LineageCast <- list()
    best_score_list_LineageCast <- c()
    tree_list <- list()

    # Iterative refinement.
    for (i in 1:maxIter) {
      renew_prob <-  nodes_sampling(tree, tree$edge, cell_state_labels, barcodes, state_lineages)
      tree_new <- subtrees_swapping(tree, tree$edge, cell_state_labels, barcodes, state_lineages, tree_renew_prob = NULL)

      score <- combined_likelihood(tree_new, barcodes_lineages[[lineage_label]], N_char = N_char,
                         cell_state_labels = state_labels_lineages[[lineage_label]], state_lineages = state_lineages,
                         state_score = NULL, barcode_score = NULL, lambda_1 = lambda1, lambda_2 = lambda2)

      ances_res <- ancestor_inference(tree_new, N_char = N_char, barcodes = barcodes_lineages[[lineage_label]],
                                  cell_state_labels = state_labels_lineages[[lineage_label]], state_lineages = state_lineages)
      barcodes <- ances_res[[1]]
      cell_state_labels <- ances_res[[2]]

      score_all <- c(score_all, maxscore)

      if (score > maxscore) {
        maxscore <- score
        best_tree <- tree_new
        tree <- tree_new
        tree_list[[length(tree_list) + 1]] <- tree
        edges <- tree$edge
        ances_res <- ancestor_inference(tree, N_char, barcodes_lineages[[lineage_label]],
                                    state_labels_lineages[[lineage_label]], state_lineages)
        barcodes <- ances_res[[1]]
        cell_state_labels <- ances_res[[2]]
      }

      if (i %% iter_repeat == 0) {
        # Check for convergence: if the score did not change during the last repeat_time iterations.
        if (length(unique(score_all[(i - iter_repeat + 1):i])) == 1) {
          best_tree_list_LineageCast[[length(best_tree_list_LineageCast) + 1]] <- best_tree
          best_score_list_LineageCast <- c(best_score_list_LineageCast, maxscore)
          # Reinitialize the tree from the initial one for this lineage.
          tree <- Trees_initial[[lineage_label]]
          tree$edge.length <- rep(1, nrow(tree$edge))
          edges <- tree$edge
          score <- combined_likelihood(tree, barcodes_lineages[[lineage_label]], N_char,
                             state_labels_lineages[[lineage_label]], state_lineages,
                             state_score = NULL, barcode_score = NULL, lambda_1 = lambda1, lambda_2 = lambda2)
          maxscore <- score
          ances_res <- ancestor_inference(tree, N_char, barcodes_lineages[[lineage_label]],
                                      state_labels_lineages[[lineage_label]], state_lineages)
          barcodes <- ances_res[[1]]
          cell_state_labels <- ances_res[[2]]
        }
      }
    }

    bestsubtreescore[[lineage_label]] <- max(best_score_list_LineageCast)
    best_index <- match(max(best_score_list_LineageCast), best_score_list_LineageCast)
    bestsubtree[[lineage_label]] <- best_tree_list_LineageCast[[best_index]]
  }

  total_time <- proc.time() - ptm_LineageCast
  return(list(bestsubtree = bestsubtree, bestsubtreescore = bestsubtreescore, total_time = total_time))
}

