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

#' Missing-safe barcode equality count
#'
#' Counts shared observed barcode states. Missing values are ignored.
#'
#' @param x,y Barcode vectors of the same length.
#' @return Integer number of shared observed states.
barcode_shared_count <- function(x, y) {
  x <- as.vector(x)
  y <- as.vector(y)
  if (length(x) != length(y)) {
    stop("barcode_shared_count(): x and y must have the same length.")
  }

  both_obs <- !is.na(x) & !is.na(y)
  sum(both_obs & x == y, na.rm = TRUE)
}

#' Calculate state shift score for cell states
#'
#' @param tree Sub-cell division tree.
#' @param edges Edges of the sub-cell division tree.
#' @param barcodes Lineage barcodes.
#' @param cell_state_labels Vector of cell state labels.
#' @param state_lineages List of state lineages.
#'
#' @return Score of cell state shifts.
state_score <- function(tree, edges, barcodes, cell_state_labels, state_lineages) {
  max_step <- max(sapply(state_lineages, length))
  shift_state_prob <- data.frame(step = 0:max_step, prob = rep(0, max_step + 1))

  for (i in seq_len(nrow(edges))) {
    edge <- edges[i, ]

    if (edge[2] > length(tree$tip.label)) {
      child_state_id <- edge[2]
    } else {
      cell_id <- tree$tip.label[edge[2]]
      child_state_id <- match(cell_id, rownames(barcodes))
    }

    parent_state <- cell_state_labels[edge[1]]
    child_state <- cell_state_labels[child_state_id]

    if (is.na(parent_state) || is.na(child_state) || parent_state == 0) {
      next
    }

    for (lineage in state_lineages) {
      if ((parent_state %in% lineage) && (child_state %in% lineage)) {
        state_distance <- match(child_state, lineage) - match(parent_state, lineage)
        if (!is.na(state_distance) && state_distance >= 0 && state_distance <= max_step) {
          shift_state_prob[shift_state_prob$step == state_distance, "prob"] <-
            shift_state_prob[shift_state_prob$step == state_distance, "prob"] + 1
        }
      }
    }
  }

  shift_state_prob <- shift_state_prob[shift_state_prob$step != 0, ]
  total <- sum(shift_state_prob$prob)
  if (total <= 0 || is.na(total)) {
    return(1)
  }

  shift_state_prob$prob <- shift_state_prob$prob / total
  score <- sum((1 / shift_state_prob$step) * shift_state_prob$prob)
  if (is.nan(score) || is.na(score)) score <- 1
  return(score)
}

# Backward-compatible alias with clearer naming.
state_continuity_score <- state_score

#' Calculate barcode shift score
#'
#' @param tree Sub-cell division tree.
#' @param N_char Total number of sites.
#' @param edges Edges of the sub-cell division tree.
#' @param barcodes Lineage barcodes.
#'
#' @return Score of barcode shifts.
barcode_score <- function(tree, N_char, edges, barcodes) {
  barcode_dist <- data.frame(diff_sites = 0:N_char, prob = rep(0, N_char + 1))
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

    diff_sites_num <- barcode_diff_count(child_barcodes_1, child_barcodes_2)
    diff_sites_num <- max(0, min(N_char, diff_sites_num))
    barcode_dist[barcode_dist$diff_sites == diff_sites_num, "prob"] <-
      barcode_dist[barcode_dist$diff_sites == diff_sites_num, "prob"] + 1
  }

  total <- sum(barcode_dist$prob)
  if (total <= 0 || is.na(total)) {
    return(1 / 0.5)
  }

  barcode_dist$prob <- barcode_dist$prob / total
  barcode_dist$diff_sites <- barcode_dist$diff_sites + 0.5
  score <- sum((1 / barcode_dist$diff_sites) * barcode_dist$prob)
  if (is.nan(score) || is.na(score)) score <- 0
  return(score)
}

# Backward-compatible alias with clearer naming.
barcode_consistency_score <- barcode_score

#' Calculate combined tree score based on cell state and barcode scores
#'
#' This function is retained under its original name for backward compatibility.
#' It computes a composite tree score.
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
#' @return Composite score of cell state and barcode shifts.
composition_score <- function(tree, barcodes, N_char, cell_state_labels, state_lineages,
                                state_score = NULL, barcode_score = NULL,
                                lambda_1 = lambda1, lambda_2 = lambda2) {
  ances_res <- ancestor_inference(tree, N_char, barcodes, cell_state_labels, state_lineages)
  barcodes1 <- ances_res[[1]]
  cell_state_labels <- ances_res[[2]]
  edges <- tree$edge

  if (is.null(state_score)) {
    state_score <- state_score(tree, edges, barcodes, cell_state_labels, state_lineages)
  }
  if (is.null(barcode_score)) {
    barcode_score <- barcode_score(tree, N_char, edges, barcodes1)
  }

  Score <- lambda_1 * state_score + lambda_2 * barcode_score
  return(Score)
}

# Clearer alias for revised manuscript/code.
combined_tree_score <- composition_score
