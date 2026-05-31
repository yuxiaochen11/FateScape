#' Compute Hamming Distance Matrix for Barcodes
#'
#' @param barcodes A matrix where each row represents a barcode.
#'
#' @return A Hamming distance matrix computed from the input barcodes.
#' @export
hamming_distance<- function(barcodes) {
  # Convert each row of the barcode matrix to a character vector and store in a list.
  barcode_list <- lapply(seq_len(nrow(barcodes)), function(i) as.character(barcodes[i, ]))
  names(barcode_list) <- rownames(barcodes)

  # Get the unique symbols present in the barcodes.
  unique_levels <- unique(unlist(barcode_list))

  # Create a phylogenetic data object using the custom levels.
  phy_data <- phyDat(barcode_list, type = "USER", levels = unique_levels)

  # Compute the Hamming distance matrix.
  distance_matrix <- dist.hamming(phy_data)

  return(distance_matrix)
}



#' Check whether a barcode entry is missing/dropout
#'
#' @param x A scalar barcode entry.
#' @return TRUE if x is NA, empty, or a dropout symbol.
#' @keywords internal
is_missing_barcode_entry <- function(x) {
  if (length(x) == 0L) return(TRUE)
  if (is.na(x)) return(TRUE)
  x_chr <- trimws(as.character(x))
  x_chr %in% c("", "NA", "NaN", "-", "?", "missing", "Missing", "MISSING")
}

#' Convert a barcode entry to numeric allele safely
#'
#' @param x A scalar barcode entry.
#' @return Numeric allele value, or NA_real_ for missing/non-convertible values.
#' @keywords internal
barcode_entry_numeric <- function(x) {
  if (is_missing_barcode_entry(x)) return(NA_real_)
  suppressWarnings(as.numeric(as.character(x)))
}

#' Calculate improved Hamming distance between two lineage barcodes
#'
#' Missing-safe version. Missing entries (NA, "-", "?", empty strings)
#' are treated as unknown/dropout states. Two missing entries do not add
#' distance, whereas one missing and one observed entry adds beta.
#'
#' @param barcode1 A vector representing the first lineage barcode.
#' @param barcode2 A vector representing the second lineage barcode.
#' @param alpha A positive number; if two sites have the same nonzero mutation,
#'   the total distance is reduced by alpha.
#' @param beta A number greater than 1; if two sites have different mutations or
#'   one site is missing and the other observed, the total distance is increased.
#' @param N_char Optional. Number of sites to consider.
#'
#' @return A finite improved Hamming distance.
#' @export
improved_hamming_distance <- function(barcode1, barcode2, alpha, beta, N_char = NULL) {
  if (is.null(N_char)) {
    N_char <- min(length(barcode1), length(barcode2))
  }
  if (length(barcode1) < N_char || length(barcode2) < N_char) {
    stop("barcode1 and barcode2 must have at least N_char entries.")
  }

  dist <- 0

  for (i in seq_len(N_char)) {
    b1_missing <- is_missing_barcode_entry(barcode1[i])
    b2_missing <- is_missing_barcode_entry(barcode2[i])

    # Both missing: unknown in both cells; do not create artificial distance.
    if (b1_missing && b2_missing) {
      next
    }

    # One missing and one observed: penalize as an uncertain difference.
    if (xor(b1_missing, b2_missing)) {
      dist <- dist + beta
      next
    }

    b1 <- barcode_entry_numeric(barcode1[i])
    b2 <- barcode_entry_numeric(barcode2[i])

    # If non-numeric but non-missing values occur, compare as strings.
    if (is.na(b1) || is.na(b2)) {
      s1 <- as.character(barcode1[i])
      s2 <- as.character(barcode2[i])
      if (identical(s1, s2) && !s1 %in% c("0", "0.0")) {
        dist <- dist - alpha
      } else if (!identical(s1, s2)) {
        dist <- dist + beta
      }
      next
    }

    # Same nonzero mutation supports shared ancestry.
    if (b1 == b2 && b1 != 0) {
      dist <- dist - alpha
    } else if (b1 != b2 && (b1 == 0 || b2 == 0)) {
      # Wild-type versus mutation.
      dist <- dist + 1
    } else if (b1 != b2 && b1 != 0 && b2 != 0) {
      # Conflicting nonzero mutations.
      dist <- dist + beta
    }
  }

  if (is.na(dist) || !is.finite(dist)) {
    stop("improved_hamming_distance produced a non-finite value; check barcode encoding.")
  }
  return(dist)
}


#' Define Hamming distance between lineage barcodes with dropout (For imputation)
#'
#' @param barcode_1 A character vector representing the first barcode.
#' @param barcode_2 A character vector representing the second barcode.
#' @param N_char Optional. Number of positions to consider in the barcode.
#'
#' @return The Hamming distance between the two barcodes with dropout.
#' @export
distance_with_dropout <- function(barcode_1, barcode_2, N_char = NULL) {
  # Check if both barcodes have the same length
  if (length(barcode_1) != length(barcode_2)) {
    stop("The two barcodes are not of the same length")
  }

  # Use the full length if N_char is not provided
  if (is.null(N_char)) {
    N_char <- length(barcode_1)
  }

  # Consider only the first N_char positions
  barcode_1 <- barcode_1[1:N_char]
  barcode_2 <- barcode_2[1:N_char]

  # Compute Hamming distance:
  # For each position, if the two values differ and neither is a dropout ("-"), count as 1.
  distance <- sum((barcode_1 != barcode_2) & (barcode_1 != "-") & (barcode_2 != "-"))

  return(distance)
}



#' Convert tree information into Newick (.nwk) format
#'
#' @param root_node The current root node.
#' @param edges A two-column matrix or data frame of edges (parent in column 1, child in column 2).
#'
#' @return A character string in Newick format.
#' @export
format_tree_newick <- function(root_node, edges) {
  parent <- root_node
  children <- edges[edges[, 1] == parent, 2]

  if (!(children[1] %in% edges[, 1]) & !(children[2] %in% edges[, 1])) {
    # Both children are terminal (leaf nodes)
    return(paste0("(", children[1], ",", children[2], ")", parent))
  } else if ((children[1] %in% edges[, 1]) & !(children[2] %in% edges[, 1])) {
    # First child is internal, second is terminal
    return(paste0("(", format_tree_newick(children[1], edges), ",", children[2], ")", parent))
  } else if (!(children[1] %in% edges[, 1]) & (children[2] %in% edges[, 1])) {
    # First child is terminal, second is internal
    return(paste0("(", children[1], ",", format_tree_newick(children[2], edges), ")", parent))
  } else {
    # Both children are internal nodes
    return(paste0("(", format_tree_newick(children[1], edges), ",", format_tree_newick(children[2], edges), ")", parent))
  }
}


#' Convert tree information into Newick (.nwk) format including time information
#'
#' @param root_node The current root node.
#' @param edges A two-column matrix or data frame of edges (parent in column 1, child in column 2).
#' @param Time A vector of time values corresponding to vertices; vertex names are taken from V(graph)$name.
#'
#' @return A character string in Newick format with branch lengths (time information).
#' @export
format_tree_newick_time <- function(root_node, edges, Time) {
  parent <- root_node
  children <- edges[edges[, 1] == parent, 2]

  # Map node names to indices in the global graph's vertex list.
  id1 <- match(children[1], V(graph)$name)
  id2 <- match(children[2], V(graph)$name)
  id_par <- match(parent, V(graph)$name)

  time1 <- Time[id1]
  time2 <- Time[id2]
  time_par <- Time[id_par]

  #print(c(time1, time2, time_par))

  if (!(children[1] %in% edges[, 1]) & !(children[2] %in% edges[, 1])) {
    return(paste0("(", children[1], ":", time1, ",", children[2], ":", time2, ")", parent, ":", time_par))
  } else if ((children[1] %in% edges[, 1]) & !(children[2] %in% edges[, 1])) {
    return(paste0("(", format_tree_newick_time(children[1], edges, Time), ",", children[2], ":", time2, ")", parent, ":", time_par))
  } else if (!(children[1] %in% edges[, 1]) & (children[2] %in% edges[, 1])) {
    return(paste0("(", children[1], ":", time1, ",", format_tree_newick_time(children[2], edges, Time), ")", parent, ":", time_par))
  } else {
    return(paste0("(", format_tree_newick_time(children[1], edges, Time), ",", format_tree_newick_time(children[2], edges, Time), ")", parent, ":", time_par))
  }
}

#' Find the depth of a given node in a tree
#'
#' @param tree A tree object containing "edge" and "edge.length" elements.
#' @param node The node (integer) whose depth is to be computed.
#' @param ncells Number of leaf nodes in the tree; the root is assumed to be ncells + 1.
#'
#' @return The depth (sum of edge lengths) from the given node to the root.
#' @export
find_node_depth <- function(tree, node, ncells) {
  root_node <- ncells + 1
  tree_edges <- tree$edge
  # Combine edges with their corresponding lengths into a matrix:
  # Column 1: parent node, Column 2: child node, Column 3: edge length.
  edge_frame <- cbind(tree_edges, tree$edge.length)

  depth <- 0
  current_node <- node

  # Traverse upward from the current node to the root.
  while (current_node != root_node) {
    # Find the row in edge_frame where the child equals the current node
    row_index <- which(edge_frame[, 2] == current_node)
    # Add the corresponding edge length to the depth
    depth <- depth + edge_frame[row_index, 3]
    # Update current_node to its parent (first column)
    current_node <- edge_frame[row_index, 1]
  }

  return(depth)
}


# FateScape missing-safe hotfix
# Source this file AFTER all FateScape R/*.R files have been sourced.
# It overrides barcode distance functions that fail when barcode entries contain NA.

missing_barcode_tokens <- c(NA_character_, "", "NA", "NaN", "nan", "?", "-", "-1")

is_missing_barcode_entry <- function(x) {
  x_chr <- as.character(x)
  is.na(x_chr) | x_chr %in% c("", "NA", "NaN", "nan", "?", "-", "-1")
}

barcode_entry_numeric <- function(x) {
  if (is_missing_barcode_entry(x)) return(NA_real_)
  suppressWarnings(as.numeric(as.character(x)))
}

improved_hamming_distance <- function(barcode1, barcode2, alpha, beta, N_char = NULL) {
  if (is.null(N_char)) {
    N_char <- min(length(barcode1), length(barcode2))
  }
  if (length(barcode1) < N_char || length(barcode2) < N_char) {
    stop("barcode1 and barcode2 must have at least N_char entries.")
  }

  dist <- 0
  for (i in seq_len(N_char)) {
    b1_missing <- is_missing_barcode_entry(barcode1[i])
    b2_missing <- is_missing_barcode_entry(barcode2[i])

    if (b1_missing && b2_missing) {
      next
    }
    if (xor(b1_missing, b2_missing)) {
      dist <- dist + beta
      next
    }

    b1 <- barcode_entry_numeric(barcode1[i])
    b2 <- barcode_entry_numeric(barcode2[i])

    if (is.na(b1) || is.na(b2)) {
      s1 <- as.character(barcode1[i])
      s2 <- as.character(barcode2[i])
      if (identical(s1, s2) && !s1 %in% c("0", "0.0")) {
        dist <- dist - alpha
      } else if (!identical(s1, s2)) {
        dist <- dist + beta
      }
      next
    }

    if (b1 == b2 && b1 != 0) {
      dist <- dist - alpha
    } else if (b1 != b2 && (b1 == 0 || b2 == 0)) {
      dist <- dist + 1
    } else if (b1 != b2 && b1 != 0 && b2 != 0) {
      dist <- dist + beta
    }
  }

  if (is.na(dist) || !is.finite(dist)) {
    stop("improved_hamming_distance produced a non-finite value; check barcode encoding.")
  }
  dist
}

barcode_to_numeric_na_for_missing <- function(x) {
  x_chr <- as.character(x)
  x_chr[is_missing_barcode_entry(x_chr)] <- NA_character_
  suppressWarnings(as.integer(as.numeric(x_chr)))
}

barcode_to_numeric_zero_for_missing <- function(x) {
  x_chr <- as.character(x)
  x_chr[is_missing_barcode_entry(x_chr)] <- "0"
  suppressWarnings(as.integer(as.numeric(x_chr)))
}

safe_improved_hamming_distance <- function(barcode1, barcode2, alpha = 0, beta = 1) {
  improved_hamming_distance(barcode1, barcode2, alpha, beta)
}

.message <- "Loaded FateScape missing-safe hotfix: improved_hamming_distance now handles NA/?/-/-1 missing barcodes."
message(.message)


