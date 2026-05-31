#' Get all descendants of a node in a phylogenetic tree (recursive)
#'
#' @param tree A phylogenetic tree of class `"phylo"`, with edge matrix.
#' @param node An integer giving the node number from which to collect descendants.
#'
#' @return An integer vector of all descendant node/leaf numbers (excluding the input node).
#' @export
getDescendants <- function(tree, node) {
  children <- tree$edge[tree$edge[, 1] == node, 2]
  if (length(children) == 0) return(NULL)
  return(c(
    children,
    unlist(lapply(children, function(child) getDescendants(tree, child)))
  ))
}


#' Compute depths of all nodes in a phylogenetic tree
#'
#' @param tree A phylogenetic tree of class `"phylo"`.
#'
#' @return A data.frame with columns:
#'   - `node`: node or tip indices (1 to Ntip + Nnode),
#'   - `depth`: distance (in number of edges) from the root.
#' @export
compute_node_depths <- function(tree) {
  root <- Ntip(tree) + 1
  all_nodes <- 1:(Ntip(tree) + tree$Nnode)
  depths <- rep(NA, length(all_nodes))
  depths[root] <- 0

  edge_df <- as.data.frame(tree$edge)
  colnames(edge_df) <- c("parent", "child")

  while (any(is.na(depths))) {
    for (i in seq_len(nrow(edge_df))) {
      p <- edge_df$parent[i]
      c <- edge_df$child[i]
      if (!is.na(depths[p]) && is.na(depths[c])) {
        depths[c] <- depths[p] + 1
      }
    }
  }

  node_df <- data.frame(node = all_nodes, depth = depths)
  return(node_df)
}


#' Compute normalized entropy of a given cell state at a specific depth
#'
#' @param tree A phylogenetic tree of class `"phylo"`.
#' @param leaf_states A named vector mapping tip labels to state labels.
#' @param state_label A single value indicating which state to compute entropy for.
#' @param node_depth_df A data.frame as returned by `compute_node_depths()`.
#' @param depth_level An integer depth at which to compute entropy.
#'
#' @return A numeric between 0 and 1 giving the normalized Shannon entropy
#'   of counts of the given state among descendants at that depth, or `NA`
#'   if there are no nodes or no counts.
#' @export
compute_entropy_at_depth <- function(tree, leaf_states, state_label, node_depth_df, depth_level) {
  nodes_d <- node_depth_df$node[node_depth_df$depth == depth_level]
  if (length(nodes_d) == 0) return(NA)

  n_s_list <- sapply(nodes_d, function(v) {
    descendants <- getDescendants(tree, v)
    leaves <- descendants[descendants <= Ntip(tree)]
    sum(leaf_states[as.character(tree$tip.label[leaves])] == state_label, na.rm = TRUE)
  })

  N_s <- sum(n_s_list)
  if (N_s == 0) return(NA)

  probs <- n_s_list / N_s
  probs <- probs[probs > 0]
  H <- -sum(probs * log(probs))
  H_norm <- if (length(probs) == 1) {
    0
  } else {
    H / log(length(probs))
  }

  return(H_norm)
}


#' Compute the entropy trajectory (path) over all depths
#'
#' @param tree A phylogenetic tree of class `"phylo"`.
#' @param leaf_states A named vector mapping tip labels to state labels.
#' @param state_label A single value indicating which state to compute entropy for.
#'
#' @return A data.frame with columns:
#'   - `depth`: integer depths from 1 to max-1,
#'   - `entropy`: normalized entropy at each depth.
#' @export
compute_entropy_path <- function(tree, leaf_states, state_label) {
  node_depth_df <- compute_node_depths(tree)
  max_depth <- max(node_depth_df$depth)
  entropy_path <- sapply(seq_len(max_depth), function(d) {
    compute_entropy_at_depth(tree, leaf_states, state_label, node_depth_df, d)
  })
  data.frame(
    depth   = seq_len(max_depth - 1),
    entropy = entropy_path[-length(entropy_path)]
  )
}


#' Get all descendants of a node in a phylogenetic tree (iterative)
#'
#' @param tree A phylogenetic tree of class `"phylo"`.
#' @param node An integer giving the node number from which to collect descendants.
#'
#' @return An integer vector of all descendant node/leaf numbers (excluding the input node).
#' @export
getDescendants_iter <- function(tree, node) {
  descendants <- integer(0)
  queue <- node
  while (length(queue) > 0) {
    current <- queue[1]
    queue <- queue[-1]
    children <- tree$edge[tree$edge[, 1] == current, 2]
    descendants <- c(descendants, children)
    queue <- c(queue, children)
  }
  return(descendants)
}


#' Compute normalized entropy at a given depth using iterative descendant search
#'
#' @inheritParams compute_entropy_at_depth
#' @param node_depth_df A data.frame as returned by `compute_node_depths()`.
#' @param depth_level An integer depth at which to compute entropy.
#'
#' @return A numeric between 0 and 1 giving the normalized Shannon entropy
#'   or `NA` if not applicable.
#' @export
compute_entropy_at_depth_alter <- function(tree, leaf_states, state_label, node_depth_df, depth_level) {
  nodes_d <- node_depth_df$node[node_depth_df$depth == depth_level]
  if (length(nodes_d) == 0) return(NA)

  n_s_list <- sapply(nodes_d, function(v) {
    descendants <- getDescendants_iter(tree, v)
    leaves <- descendants[descendants <= Ntip(tree)]
    sum(leaf_states[as.character(tree$tip.label[leaves])] == state_label, na.rm = TRUE)
  })

  N_s <- sum(n_s_list)
  if (N_s == 0) return(NA)

  probs <- n_s_list / N_s
  probs <- probs[probs > 0]
  H <- -sum(probs * log(probs))
  H_norm <- if (length(probs) == 1) {
    0
  } else {
    H / log(length(probs))
  }
  return(H_norm)
}


#' Compute the entropy trajectory using iterative descendant search
#'
#' @inheritParams compute_entropy_path
#'
#' @return A data.frame with columns `depth` and `entropy`.
#' @export
compute_entropy_path_alter <- function(tree, leaf_states, state_label) {
  node_depth_df <- compute_node_depths(tree)
  max_depth <- max(node_depth_df$depth)
  entropy_path <- sapply(seq_len(max_depth), function(d) {
    message("Computing depth: ", d)
    compute_entropy_at_depth_alter(tree, leaf_states, state_label, node_depth_df, d)
  })
  data.frame(
    depth   = seq_len(max_depth - 1),
    entropy = entropy_path[-length(entropy_path)]
  )
}


#' Compute discrete entropy slope over a sliding interval
#'
#' @param tree A phylogenetic tree of class `"phylo"`.
#' @param leaf_states A named vector mapping tip labels to state labels.
#' @param state_label A single value indicating which state to compute entropy for.
#' @param interval Integer specifying the step over which to compute slope (default 3).
#'
#' @return A numeric vector of entropy differences (`H[t+interval] - H[t]`) or `NA`
#'   if the number of depths ≤ interval.
#' @export
compute_entropy_slope <- function(tree, leaf_states, state_label, interval = 3) {
  entropy_df <- compute_entropy_path(tree, leaf_states, state_label)
  valid <- !is.na(entropy_df$entropy)
  H <- entropy_df$entropy[valid]
  n <- length(H)

  if (n <= interval) return(NA)

  idx <- seq_len(n - interval)
  slope <- H[idx + interval] - H[idx]

  return(slope)
}
