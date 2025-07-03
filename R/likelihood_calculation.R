#' Calculate state shift probability for cell states
#'
#' @param tree Sub-cell division tree.
#' @param edges Edges of the sub-cell division tree.
#' @param barcodes Lineage barcodes.
#' @param cell_state_labels Vector of cell state labels.
#' @param state_lineages List of state lineages.
#'
#' @return Score of cell state shifts.
state_transition_likelihood <- function(tree, edges, barcodes, cell_state_labels, state_lineages) {
  max_step <- max(sapply(state_lineages, length))
  # Create a data frame to count occurrences for each state shift step.
  shift_state_prob <- data.frame(step = 0:max_step, prob = rep(0, max_step + 1))

  # Loop through each edge of the tree.
  for (i in 1:nrow(edges)) {
    edge <- edges[i, ]

    # Determine the child state index based on whether the node is a tip or internal node.
    if (edge[2] > length(tree$tip.label)) {
      child_state_id <- edge[2]
    } else {
      cell_id <- tree$tip.label[edge[2]]
      child_state_id <- match(cell_id, rownames(barcodes))
    }

    parent_state <- cell_state_labels[edge[1]]
    child_state <- cell_state_labels[child_state_id]

    # Process only if parent state is non-zero.
    if (parent_state != 0) {
      # For each state lineage, check if both states belong to it.
      for (lineage in state_lineages) {
        if ((parent_state %in% lineage) && (child_state %in% lineage)) {
          state_distance <- match(child_state, lineage) - match(parent_state, lineage)
          shift_state_prob[shift_state_prob$step == state_distance, "prob"] <-
            shift_state_prob[shift_state_prob$step == state_distance, "prob"] + 1
        }
      }
    }
  }

  # Remove zero step (no state shift) and normalize probabilities.
  shift_state_prob <- shift_state_prob[shift_state_prob$step != 0, ]
  shift_state_prob$prob <- shift_state_prob$prob / sum(shift_state_prob$prob)

  # Calculate the score as a weighted average of the inverse of the shift steps.
  score <- sum((1 / shift_state_prob$step) * shift_state_prob$prob)
  if (is.nan(score)) { score <- 1 }
  return(score)
}


#' Calculate barcode shift probability
#'
#' @param tree Sub-cell division tree.
#' @param N_char Total number of sites.
#' @param edges Edges of the sub-cell division tree.
#' @param barcodes Lineage barcodes.
#'
#' @return Score of barcode shifts.
barcode_similarity_likelihood <- function(tree, N_char, edges, barcodes) {
  # Create a data frame to count occurrences for each number of differing sites.
  barcode_dist <- data.frame(diff_sites = 0:N_char, prob = rep(0, N_char + 1))

  # Get the unique parent nodes from the tree edges.
  parent_nodes <- unique(edges[, 1])

  for (par in parent_nodes) {
    children_nodes <- edges[edges[, 1] == par, 2]

    # Determine the first child index.
    if (children_nodes[1] > length(tree$tip.label)) {
      child_id_1 <- children_nodes[1]
    } else {
      cell_id_1 <- tree$tip.label[children_nodes[1]]
      child_id_1 <- match(cell_id_1, rownames(barcodes))
    }

    # Determine the second child index.
    if (children_nodes[2] > length(tree$tip.label)) {
      child_id_2 <- children_nodes[2]
    } else {
      cell_id_2 <- tree$tip.label[children_nodes[2]]
      child_id_2 <- match(cell_id_2, rownames(barcodes))
    }

    child_barcodes_1 <- barcodes[child_id_1, ]
    child_barcodes_2 <- barcodes[child_id_2, ]

    diff_sites_num <- sum(child_barcodes_1 != child_barcodes_2)
    barcode_dist[barcode_dist$diff_sites == diff_sites_num, "prob"] <-
      barcode_dist[barcode_dist$diff_sites == diff_sites_num, "prob"] + 1
  }

  # Normalize the counts to get probabilities.
  barcode_dist$prob <- barcode_dist$prob / sum(barcode_dist$prob)
  # Add 0.5 to each difference count to avoid division by zero.
  barcode_dist$diff_sites <- barcode_dist$diff_sites + 0.5
  score <- sum((1 / barcode_dist$diff_sites) * barcode_dist$prob)
  return(score)
}


#' Calculate combined tree score based on cell state and barcode shift scores
#'
#' @param tree Sub-cell division tree.
#' @param barcodes Lineage barcodes.
#' @param N_char Total number of sites.
#' @param cell_state_labels Vector of cell state labels.
#' @param state_lineages List of state lineages.
#' @param state_score Precomputed cell state score (optional).
#' @param barcode_score Precomputed barcode score (optional).
#' @param lambda_1 Weight for the cell state score.
#' @param lambda_2 Weight for the barcode score.
#'
#' @return Combined score of cell state and barcode shifts.
combined_likelihood <- function(tree, barcodes, N_char, cell_state_labels, state_lineages,
                      state_score = NULL, barcode_score = NULL, lambda_1 = lambda1, lambda_2 = lambda2) {
  # Perform ancestral inference to update barcodes and cell state labels.
  ances_res <- ancestor_inference(tree, N_char, barcodes, cell_state_labels, state_lineages)
  barcodes1 <- ances_res[[1]]
  cell_state_labels <- ances_res[[2]]
  edges <- tree$edge

  # Compute state shift probability if not provided.
  if (is.null(state_score)) {
    state_score <- state_transition_likelihood(tree, edges, barcodes, cell_state_labels, state_lineages)
  }
  # Compute barcode shift probability if not provided.
  if (is.null(barcode_score)) {
    barcode_score <- barcode_similarity_likelihood(tree, N_char, edges, barcodes1)
  }

  # Calculate the combined score using the specified weights.
  Score <- lambda_1 * state_score + lambda_2 * barcode_score
  return(Score)
}
