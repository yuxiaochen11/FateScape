#' Find the distance between sibling cells
#'
#' @param cell A tip label (cell ID) for which to find the sibling distance.
#' @param tree A sub-cell division tree.
#' @param lineage The state lineage identifier (e.g., "L1") that the cell belongs to.
#' @param alpha Alpha parameter for the improved Hamming distance.
#' @param beta Beta parameter for the improved Hamming distance.
#'
#' @return The improved Hamming distance between the sibling cells.
#' @export
sibling_cells_distance <- function(cell, tree, lineage, alpha, beta) {
  # Get tip labels and edge matrix from the tree.
  tips <- tree$tip.label
  edges <- tree$edge

  # Identify the parent node for the specified cell.
  cell_index <- which(tips == cell)
  parent_node <- edges[edges[, 2] == cell_index, 1]

  # Find sibling nodes: all children of the parent node.
  siblings <- edges[edges[, 1] == parent_node, 2]

  # Perform ancestral inference for the given lineage.
  # Note: N_char, barcodes_lineages, state_labels_lineages, and state_lineages are assumed to be globally defined.
  ances_res <- ancestor_inference(tree, N_char, barcodes_lineages[[lineage]],
                              state_labels_lineages[[lineage]], state_lineages)

  # Extract barcodes for the sibling cells.
  sibling_barcodes <- ances_res[[1]][siblings, ]

  # Compute the improved Hamming distance between the two sibling barcodes.
  dist <- improved_hamming_distance(sibling_barcodes[1, ], sibling_barcodes[2, ], alpha, beta)

  return(dist)
}


#' Drop non-minimum tips from tree list
#'
#' @param tree_list List of trees.
#' @param dist_list Data frame with columns "lineage" and "dist" containing distances.
#' @param cell Cell ID (tip) to be dropped.
#'
#' @return Modified tree_list after dropping the specified tip from trees not in the best (minimum-distance) lineage.
#' @export
drop_non_minimum_tips <- function(tree_list, dist_list, cell) {
  dist <- dist_list$dist
  min_lineage <- dist_list[which.min(dist), 1]
  if (length(min_lineage) != 1) {
    min_lineage <- sample(min_lineage, 1)
  }
  drop_lineages <- dist_list[dist_list$lineage != min_lineage, 1]
  for (j in drop_lineages) {
    tree_list[[j]] <- drop.tip(tree_list[[j]], cell)
  }
  return(tree_list)
}


#' Drop duplicated leaf nodes from tree list
#'
#' @param tree_list List of raw sub-cell division trees.
#' @param barcodes Matrix of lineage barcodes.
#' @param cell_lineages Cell lineages (unused in the current function but kept for compatibility).
#' @param state_lineage A list where each element is a vector of state lineage labels for a cell.
#' @param alpha Alpha parameter for computing sibling distance.
#' @param beta Beta parameter for computing sibling distance.
#'
#' @return Modified tree_list with duplicated leaf nodes dropped.
#' @export
drop_duplicated_tips <- function(tree_list, barcodes, cell_lineages, state_lineage, alpha, beta) {
  for (i in seq_along(state_lineage)) {
    if (length(state_lineage[[i]]) > 1) {
      cell <- names(state_lineage)[i]
      lineages <- state_lineage[[i]]
      # Compute sibling distances for each lineage using sapply.
      dist_vec <- sapply(lineages, function(lin) {
        sibling_cells_distance(cell, tree_list[[lin]], lineage = lin, alpha, beta)
      })
      dist_list <- data.frame(lineage = lineages, dist = dist_vec)
      tree_list <- drop_non_minimum_tips(tree_list, dist_list, cell)
    }
  }
  return(tree_list)
}

#' Get first element (actually returns the second element)
#'
#' @param x A list.
#'
#' @return The second element of x.
#' @export
get_first_element <- function(x) {
  return(x[[2]])
}


#' Get tip labels from a tree
#'
#' @param x A tree object.
#'
#' @return A vector of tip labels.
#' @export
get_tips <- function(x) {
  return(unlist(x$tip.label))
}


#' Get the root barcode from a tree
#'
#' @param x A tree object with a root.barcode element.
#'
#' @return The root barcode.
#' @export
get_root_barcode <- function(x) {
  return(x$root.barcode)
}


#' Get root barcodes from a list of trees
#'
#' @param x A list of tree objects.
#'
#' @return A list of root barcodes.
#' @export
subtree_root_bar <- function(x) {
  return(lapply(x, get_root_barcode))
}


#' Subdivide a sub-cell division tree into sub-subtrees
#'
#' @param sub_tree A sub-cell division tree.
#' @param state A data frame with cell state information.
#' @param barcodes A matrix of lineage barcodes.
#' @param subid Identifier for the sub-cell division tree.
#' @param itr_time Iteration time (for recursive calls).
#'
#' @return A list of sub-subtrees.
#' @export
prune_cell_tree <- function(sub_tree, state, barcodes, subid, itr_time) {
  sub_tree <- Preorder(sub_tree)
  ntips <- length(sub_tree$tip.label)

  st_barleaves <- barcodes[sub_tree$tip.label, ]
  st_stateleaves <- state$cluster[match(sub_tree$tip.label, state$cell_id)]

  # Get ancestral barcodes for internal nodes of sub_tree
  st_barances <- ancestor_inference(sub_tree, N_char, st_barleaves, st_stateleaves, state_lineages)[[1]][(ntips + 1):(2 * ntips - 1), , drop = FALSE]

  bar_ances <- matrix(as.integer(replace(st_barances, st_barances == "-", 0)), nrow = nrow(st_barances))
  bar_leaves <- matrix(as.integer(replace(st_barleaves, st_barleaves == "-", 0)), nrow = nrow(st_barleaves))
  barances <- rbind(bar_leaves, bar_ances)
  rownames(barances) <- c(rownames(st_barleaves), rownames(st_barances))

  candidate_split_nodes_name <- rownames(barances)[which(rowSums(barances) == 0)]
  # Extract the numeric id from the candidate names using get_first_element
  id <- as.integer(lapply(strsplit(candidate_split_nodes_name, "_"), get_first_element))

  subsubtrees_list <- list()
  if (length(id) == 0) {
    sub_tree$root.barcode <- barances[1, , drop = FALSE]
    sub_tree$root.edge <- 1
    subsubtrees_list[[rownames(barances)[1]]] <- sub_tree
  } else {
    for (i in 1:length(id)) {
      children_index <- sub_tree$edge[sub_tree$edge[, 1] == id[i], 2]
      if ((children_index[1] %in% id) || (children_index[2] %in% id)) {
        next
      } else {
        children_names <- character(2)
        for (k in 1:2) {
          if (children_index[k] < ntips) {
            children_names[k] <- rownames(barances)[children_index[k]]
          } else {
            children_names[k] <- paste("node", children_index[k], sep = "_")
          }
        }
        subtree_1 <- Subtree(sub_tree, children_index[1])
        subtree_1$root.barcode <- barances[children_names[1], , drop = FALSE]
        subtree_1$root.edge <- 1
        key1 <- paste(itr_time, subid, "node", children_index[1], sep = "_")
        subsubtrees_list[[key1]] <- subtree_1

        subtree_2 <- Subtree(sub_tree, children_index[2])
        subtree_2$root.barcode <- barances[children_names[2], , drop = FALSE]
        subtree_2$root.edge <- 1
        key2 <- paste(itr_time, subid, "node", children_index[2], sep = "_")
        subsubtrees_list[[key2]] <- subtree_2
      }
    }

    if (length(id) != 1) {
      drop_tips_vec <- unlist(lapply(subsubtrees_list, get_tips))
      subsubtree_left <- drop.tip(sub_tree, drop_tips_vec, trim.internal = TRUE, subtree = FALSE, root.edge = 1)
      if (!is.null(subsubtree_left)) {
        barleft <- barcodes[subsubtree_left$tip.label, ]
        stateleft <- state$cluster[match(subsubtree_left$tip.label, state$cell_id)]
        left_barances <- ancestor_inference(subsubtree_left, N_char, barleft, stateleft, state_lineages)[[1]][(length(subsubtree_left$tip.label) + 1):(2 * length(subsubtree_left$tip.label) - 1), , drop = FALSE]
        leftbarances <- matrix(as.integer(replace(left_barances, left_barances == "-", 0)), nrow = nrow(left_barances))
        rownames(leftbarances) <- rownames(left_barances)
        subsubtree_left$root.barcode <- leftbarances[1, , drop = FALSE]
        subsubtree_left$root.edge <- 1
        key_left <- paste(subid, rownames(leftbarances)[1], sep = "_")
        subsubtrees_list[[key_left]] <- subsubtree_left
      }
    }
  }

  m0 <- t(as.matrix(rep(0, N_char), nrow = 1, ncol = N_char))
  # If the root barcode of the last subtree in the list is not all zeros, return the list;
  # otherwise, if subsubtree_left is not null, recursively subdivide it.
  if (sum(subtree_root_bar(subsubtrees_list)[[length(subsubtrees_list)]] == m0) != N_char) {
    return(subsubtrees_list)
  } else if (is.null(subsubtree_left)) {
    return(subsubtrees_list)
  } else {
    itr_time <- itr_time + 1
    return(c(subsubtrees_list[-length(subsubtrees_list)], prune_cell_tree(subsubtree_left, state, barcodes, subid, itr_time)))
  }
}


#' Count common mutations between two barcodes
#'
#' @param barcode1 A vector representing the first barcode.
#' @param barcode2 A vector representing the second barcode.
#' @param N_char Optional. Number of sites to consider.
#'
#' @return Count of sites with the same mutation (excluding zeros).
#' @export
count_common_mutations <- function(barcode1, barcode2, N_char = NULL) {
  if (is.null(N_char)) {
    N_char <- length(barcode1)
  }
  CM <- 0
  for (i in 1:N_char) {
    if (barcode1[i] == barcode2[i] && barcode1[i] != 0) {
      CM <- CM + 1
    }
  }
  return(CM)
}


#' Count uncommon mutations between two barcodes
#'
#' @param barcode1 A vector representing the first barcode.
#' @param barcode2 A vector representing the second barcode.
#' @param N_char Optional. Number of sites to consider.
#'
#' @return Count of sites with different mutations (both nonzero).
#' @export
count_uncommon_mutations <- function(barcode1, barcode2, N_char = NULL) {
  if (is.null(N_char)) {
    N_char <- length(barcode1)
  }
  UM <- 0
  for (i in 1:N_char) {
    if (barcode1[i] != barcode2[i] && (barcode1[i] * barcode2[i] != 0)) {
      UM <- UM + 1
    }
  }
  return(UM)
}


#' Get subtree root barcodes for each best subtree in a lineage
#'
#' @param bestsubtree A list of best subtrees.
#' @param state Data frame with cell state information.
#' @param barcodes Matrix of lineage barcodes.
#' @param l_sl Length (number) of state lineages.
#'
#' @return A vector of root barcodes for the sub-subtrees.
#' @export
get_subtree_root_barcodes <- function(bestsubtree, state, barcodes, l_sl) {
  root_bar <- c()
  for (i in 1:l_sl) {
    root_bar <- c(root_bar, subtree_root_bar(prune_cell_tree(bestsubtree[[i]], state, barcodes, i, 1)))
  }
  return(root_bar)
}


#' Build a common mutation matrix between root barcodes of sub-subtrees
#'
#' @param subtrees_rootbar A list of root barcodes for sub-subtrees.
#'
#' @return A matrix of common mutation counts between each pair of sub-subtree root barcodes.
#' @export
common_mutation_matrix <- function(subtrees_rootbar) {
  CM <- matrix(0, nrow = length(subtrees_rootbar), ncol = length(subtrees_rootbar))
  rownames(CM) <- colnames(CM) <- names(subtrees_rootbar)
  for (i in 1:length(subtrees_rootbar)) {
    for (j in 1:length(subtrees_rootbar)) {
      CM[i, j] <- count_common_mutations(subtrees_rootbar[[i]], subtrees_rootbar[[j]])
      if (i >= j) {
        CM[i, j] <- NA
      }
    }
  }
  return(CM)
}


#' Rank subtrees based on common mutation counts
#'
#' @param CM A common mutation matrix.
#'
#' @return A list containing the nodes rank matrix and nodes weight matrix.
#' @export
subtrees_rank <- function(CM) {
  Nodes <- rownames(CM)
  N <- length(Nodes)
  Nodes_rank <- matrix(NA, nrow = N, ncol = N)
  Nodes_weight <- matrix(NA, nrow = N, ncol = N)
  rownames(Nodes_rank) <- rownames(Nodes_weight) <- Nodes
  for (i in 1:(N - 1)) {
    num_positive <- sum(CM[i, ] > 0, na.rm = TRUE)
    if (num_positive == 0) {
      next
    } else {
      ordered_indices <- order(CM[i, ], decreasing = TRUE)
      Nodes_rank[i, 1:num_positive] <- Nodes[ordered_indices][1:num_positive]
      Nodes_weight[i, 1:num_positive] <- CM[i, ][ordered_indices][1:num_positive]
    }
  }
  return(list(Nodes_rank, Nodes_weight))
}


#' Generate sub-subtrees from best subtrees
#'
#' @param bestsubtree A list of best subtrees.
#' @param state Data frame with cell state information.
#' @param barcodes Matrix of lineage barcodes.
#' @param l_sl Length (number) of state lineages.
#'
#' @return A vector (list) of sub-subtrees.
#' @export
decompose_subtrees <- function(bestsubtree, state, barcodes, l_sl) {
  subtrees <- c()
  for (i in 1:l_sl) {
    subtrees <- c(subtrees, prune_cell_tree(bestsubtree[[i]], state, barcodes, i, 1))
  }
  return(subtrees)
}


#' Merge sub-cell division trees in order based on nodes rank and weight
#'
#' @param CM Common mutation matrix.
#' @param Nodes_rank Matrix of nodes ranking.
#' @param Nodes_weight Matrix of nodes weight.
#' @param subsubtrees List of sub-subtrees.
#' @param bind_tree_list List of trees to bind.
#' @param nsubtree Current number/index for subtrees.
#'
#' @return A merged cell division tree.
#' @export
merge_subcell_trees_ordered <- function(CM, Nodes_rank, Nodes_weight, subsubtrees, bind_tree_list, nsubtree) {
  for (i in 1:length(subsubtrees)) {
    subsubtrees[[i]]$edge.length <- rep(1, nrow(subsubtrees[[i]]$edge))
  }
  Nodes <- rownames(Nodes_rank)
  Nodes_remain <- Nodes
  N <- length(Nodes)
  node1 <- Nodes_rank[which.max(Nodes_weight)]
  node2 <- Nodes[which.max(Nodes_weight) %% N]
  group <- na.omit(unique(c(node1, node2, Nodes_rank[node1, ], Nodes_rank[node2, ])))
  ngroup <- length(group)
  bind_tree_list[[nsubtree]] <- subsubtrees[[group[1]]]
  bind_tree_list[[nsubtree]]$root.edge <- 1
  if (ngroup > 1) {
    for (i in 2:ngroup) {
      bind_tree_list[[nsubtree]] <- bind.tree(bind_tree_list[[nsubtree]], subsubtrees[[group[i]]], where = "root", .1)
      bind_tree_list[[nsubtree]]$root.edge <- 1
    }
  }
  Nodes_remain <- Nodes_remain[!Nodes_remain %in% group]
  if (length(Nodes_remain) == 0) {
    final_tree <- bind_tree_list[[1]]
    final_tree$root.edge <- 1
    if (length(bind_tree_list) > 1) {
      for (i in 2:length(bind_tree_list)) {
        final_tree <- bind.tree(final_tree, bind_tree_list[[i]], where = "root", .1)
        final_tree$root.edge <- 1
      }
    }
    return(final_tree)
  } else if (length(Nodes_remain) != 0 && sum(CM[Nodes_remain, Nodes_remain] > 0, na.rm = TRUE) > 0) {
    new_rank_weight <- subtrees_rank(CM[Nodes_remain, Nodes_remain])
    return(merge_subcell_trees_ordered(CM, Nodes_rank = new_rank_weight[[1]], Nodes_weight = new_rank_weight[[2]],
                                       subsubtrees, bind_tree_list = bind_tree_list, nsubtree = length(bind_tree_list) + 1))
  } else {
    final_tree <- bind_tree_list[[1]]
    final_tree$root.edge <- 1
    if (length(bind_tree_list) > 1) {
      for (i in 2:length(bind_tree_list)) {
        final_tree <- bind.tree(final_tree, bind_tree_list[[i]], where = "root", .1)
        final_tree$root.edge <- 1
      }
    }
    for (j in 1:length(Nodes_remain)) {
      final_tree <- bind.tree(final_tree, subsubtrees[[Nodes_remain[j]]], where = "root", .1)
      final_tree$root.edge <- 1
    }
    return(final_tree)
  }
}


#' Group sub-cell division trees for merging based on common mutations
#'
#' @param CM Common mutation matrix.
#' @param Nodes_rank Matrix of nodes ranking.
#' @param Nodes_weight Matrix of nodes weight.
#' @param subsubtrees List of sub-subtrees.
#' @param bind_tree_list List of trees to bind.
#' @param nsubtree Current number/index for subtrees.
#'
#' @return A list containing the merged tree group and remaining nodes.
#' @export
group_subcell_trees <- function(CM, Nodes_rank, Nodes_weight, subsubtrees, bind_tree_list, nsubtree) {
  Nodes <- rownames(Nodes_rank)
  Nodes_remain <- Nodes
  N <- length(Nodes)
  node1 <- Nodes_rank[which.max(Nodes_weight)]
  node2 <- Nodes[which.max(Nodes_weight) %% N]
  group <- na.omit(unique(c(node1, node2, Nodes_rank[node1, ], Nodes_rank[node2, ])))
  ngroup <- length(group)
  bind_tree_list[[nsubtree]] <- subsubtrees[[group[1]]]
  bind_tree_list[[nsubtree]]$root.edge <- 1
  if (ngroup > 1) {
    for (i in 2:ngroup) {
      bind_tree_list[[nsubtree]] <- bind.tree(bind_tree_list[[nsubtree]], subsubtrees[[group[i]]], where = "root", .1)
      bind_tree_list[[nsubtree]]$root.edge <- 1
    }
  }
  Nodes_remain <- Nodes_remain[!Nodes_remain %in% group]
  if (length(Nodes_remain) == 0) {
    return(bind_tree_list)
  } else if (length(Nodes_remain) != 0 && sum(CM[Nodes_remain, Nodes_remain] > 0, na.rm = TRUE) > 0) {
    new_rank_weight <- subtrees_rank(CM[Nodes_remain, Nodes_remain])
    return(group_subcell_trees(CM, Nodes_rank = new_rank_weight[[1]], Nodes_weight = new_rank_weight[[2]],
                               subsubtrees, bind_tree_list = bind_tree_list, nsubtree = length(bind_tree_list) + 1))
  } else if (length(Nodes_remain) != 0 && sum(CM[Nodes_remain, Nodes_remain] > 0, na.rm = TRUE) == 0) {
    return(list(bind_tree_list, Nodes_remain))
  }
}


#' Merge sub-cell division trees using Ward's method
#'
#' @param subtrees_rootbar A list of root barcodes for sub-subtrees.
#' @param subsubtrees A list of sub-subtrees.
#'
#' @return A merged cell division tree.
#' @export
merge_subcell_trees_ward <- function(subtrees_rootbar, subsubtrees) {
  for (i in 1:length(subsubtrees)) {
    subsubtrees[[i]]$edge.length <- rep(1, nrow(subsubtrees[[i]]$edge))
  }
  b <- do.call(rbind, subtrees_rootbar)
  rownames(b) <- names(subtrees_rootbar)
  barcodes_lineage <- b
  D <- hamming_distance(barcodes_lineage)
  tree_WD <- nodes_clustering(D)
  for (i in names(subsubtrees)) {
    position <- which(i == tree_WD$tip.label)
    tree_WD <- bind.tree(tree_WD, subsubtrees[[i]], position)
  }
  return(tree_WD)
}


#' Get subtree cell types
#'
#' @param tree_groups A list containing two groups of trees. For the first group, each element is a tree object;
#'   for the second group, each element is an index referring to a tree in the global "subsubtrees" list.
#'
#' @return A list of cell type vectors for each subtree, with names "subtree_1", "subtree_2", etc.
#' @export
get_subtree_celltypes <- function(tree_groups) {
  subtree_celltypes_list <- list()

  # Process group 1 and group 2 separately.
  for (j in 1:2) {
    if (j == 1) {
      for (i in seq_along(tree_groups[[j]])) {
        # Retrieve the cell types for the tips of the tree.
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

  # Rename list elements to "subtree_1", "subtree_2", etc.
  names(subtree_celltypes_list) <- paste("subtree", seq_along(subtree_celltypes_list), sep = "_")

  return(subtree_celltypes_list)
}


#' Get unique cell types from subtrees
#'
#' @param subtree_celltypes A list of cell type vectors (e.g., obtained from get_subtree_celltypes).
#'
#' @return A list where each element corresponds to a cell type that appears in exactly one subtree.
#'   The element names indicate the cell type and the value is the name of the subtree that contains it.
#' @export
get_unique_celltypes <- function(subtree_celltypes) {
  all_elements <- unique(unlist(subtree_celltypes))
  element_counts <- table(unlist(subtree_celltypes))
  unique_elements <- list()

  for (element in names(element_counts)) {
    if (element_counts[element] == 1) {
      vec_index <- which(sapply(subtree_celltypes, function(x) element %in% x))
      if (length(vec_index) == 1) {
        unique_elements[[element]] <- names(subtree_celltypes)[vec_index]
      }
    }
  }

  return(unique_elements)
}


#' Compute a time-scaled phylogenetic tree based on lineage barcodes
#'
#' @param tree A sub-cell division tree.
#' @param N_char Total number of sites.
#' @param barcodes A matrix of lineage barcodes.
#' @param ncells Optional. Number of cells (defaults to the number of tip labels in the tree).
#' @param Nnodes Optional. Number of nodes (defaults to 2 * ncells - 2).
#' @param edges Optional. Edge information (ignored if NULL).
#'
#' @return A time-scaled phylogenetic tree.
#' @export
phylogenetic_tree <- function(tree, N_char, barcodes, ncells = NULL, Nnodes = NULL, edges = NULL) {
  # Set number of cells if not provided.
  if (is.null(ncells)) {
    ncells <- length(tree$tip.label)
  }

  # Set number of nodes if not provided.
  if (is.null(Nnodes)) {
    Nnodes <- 2 * ncells - 2
  }

  # Create an edge data frame with parent, child, and initial edge lengths set to 0.
  edges <- data.frame(
    par = tree$edge[, 1],
    child = tree$edge[, 2],
    length = rep(0, nrow(tree$edge))
  )

  # Perform ancestral inference to obtain ancestral barcodes.
  barcodes_ances <- ancestor_inference_time_scale(tree, N_char, barcodes)

  # Iterate over each parent node in the edge data frame.
  for (par in edges$par) {
    children <- edges$child[edges$par == par]
    if (length(children) == 2) {
      for (j in 1:2) {
        child <- children[j]
        # Calculate the improved Hamming distance between parent's and child's barcodes.

        H_D <-  improved_hamming_distance(barcodes_ances[par, ], barcodes_ances[child, ],0,1)
        # Set the edge length: if no difference then 0.5, else proportional to H_D.
        if (H_D == 0) {
          edges$length[edges$child == child] <- 0.5
        } else {
          edges$length[edges$child == child] <- H_D * 1
        }
      }
    }
  }

  # Update the tree edge lengths with computed distances.
  tree$edge.length <- edges$length

  # Compute and return a time-scaled tree using the chronos function.
  time_scale_tree <- chronos(tree)
  return(time_scale_tree)
}


