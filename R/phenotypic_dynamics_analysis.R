getDescendants <- function(tree, node) {
  children <- tree$edge[tree$edge[, 1] == node, 2]
  if (length(children) == 0) return(NULL)
  return(c(children, unlist(lapply(children, function(child) getDescendants(tree, child)))))
}

compute_node_depths <- function(tree) {
  root <- Ntip(tree) + 1
  all_nodes <- 1:(Ntip(tree) + tree$Nnode)
  depths <- rep(NA, length(all_nodes))
  depths[root] <- 0

  edge_df <- as.data.frame(tree$edge)
  colnames(edge_df) <- c("parent", "child")

  while (any(is.na(depths))) {
    for (i in 1:nrow(edge_df)) {
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
  if (length(probs) == 1) {H_norm <- 0}else{
    H_norm <- H / log(length(probs))
  }

  return(H_norm)
}

compute_entropy_path <- function(tree, leaf_states, state_label) {
  node_depth_df <- compute_node_depths(tree)
  max_depth <- max(node_depth_df$depth)
  entropy_path <- sapply(1:max_depth, function(d) {
    compute_entropy_at_depth(tree, leaf_states, state_label, node_depth_df, d)
  })
  return(data.frame(depth = 1:(max_depth-1), entropy = entropy_path[-length(entropy_path)]))
}

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
  H_norm <- if (length(probs) == 1) 0 else H / log(length(probs))
  return(H_norm)
}

compute_entropy_path_alter <- function(tree, leaf_states, state_label) {
  node_depth_df <- compute_node_depths(tree)
  max_depth <- max(node_depth_df$depth)
  entropy_path <- sapply(1:max_depth, function(d) {
    cat("Computing depth:", d, "\n")
    compute_entropy_at_depth_alter(tree, leaf_states, state_label, node_depth_df, d)
  })
  return(data.frame(depth = 1:(max_depth - 1), entropy = entropy_path[-length(entropy_path)]))
}

compute_entropy_slope <- function(tree, leaf_states, state_label, interval = 3) {
  entropy_df <- compute_entropy_path(tree, leaf_states, state_label)
  valid <- !is.na(entropy_df$entropy)
  H <- entropy_df$entropy[valid]
  n <- length(H)

  if (n <= interval) return(NA)

  idx <- 1:(n - interval)
  slope <- H[idx + interval] - H[idx]

  return(slope)
}
