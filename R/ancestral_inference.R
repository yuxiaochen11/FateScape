#' Infer ancestral lineage barcodes and cell states
#'
#' @param tree A sub-cell division tree.
#' @param N_char Number of characters in the lineage barcode.
#' @param barcodes Lineage barcodes of the leaf nodes.
#' @param cell_state_labels Cell states of the leaf nodes.
#' @param state_lineages A list of state lineages.
#'
#' @return A list containing ancestral lineage barcodes and cell states.
#' @export
ancestor_inference <- function(tree, N_char, barcodes, cell_state_labels, state_lineages) {
  # Retrieve all internal (parent) nodes and sort them in decreasing order
  nodes_internal <- unique(tree$edge[, 1])
  nodes_internal <- sort(nodes_internal, decreasing = TRUE)

  # Create a matrix for internal nodes with default barcode values (0)
  barcodes2 <- matrix(0, nrow = length(nodes_internal), ncol = N_char)
  rownames(barcodes2) <- paste("node", sort(nodes_internal, decreasing = FALSE), sep = "_")

  # Append the internal nodes' barcode matrix to the existing barcodes
  barcodes <- rbind(barcodes, barcodes2)

  # Append default cell state (0) for the internal nodes
  cell_state_labels <- c(cell_state_labels, rep(0, length(nodes_internal)))

  # Iterate over each internal node to infer barcode and cell state
  for (node in nodes_internal) {
    # Get the indices of child nodes for the current parent node
    children_index <- tree$edge[tree$edge[, 1] == node, 2]

    children <- c()
    # Process each child node index
    for (cell_id in children_index) {
      if (cell_id > length(tree$tip.label)) {
        # For internal child nodes, use the cell_id as is
        cell_row <- cell_id
      } else {
        # For tip nodes, obtain the corresponding label and find its row in the barcodes matrix
        cell_id <- tree$tip.label[cell_id]
        cell_row <- match(cell_id, rownames(barcodes))
      }
      children <- c(children, cell_row)
    }

    # Extract barcode and cell state information for the child nodes
    barcode_children <- barcodes[children, ]
    cell_state_children <- cell_state_labels[children]

    # Infer the parent's cell state based on the state lineages of its children
    for (lineage in state_lineages) {
      if ((cell_state_children[1] %in% lineage) && (cell_state_children[2] %in% lineage)) {
        cell_state_labels[node] <- lineage[min(match(cell_state_children[1], lineage),
                                               match(cell_state_children[2], lineage))]
      }
    }

    # Infer the parent's barcode for each character position
    for (i in 1:N_char) {
      characters <- barcode_children[, i]
      if (length(unique(characters)) == 1) {
        barcodes[node, i] <- unique(characters)
      } else {
        barcodes[node, i] <- "0"
      }
    }
  }

  # Return the inferred barcodes and cell states
  return(list(barcodes, cell_state_labels))
}


#' Infer ancestral barcodes with time scaling
#'
#' @param tree A sub-cell division tree.
#' @param N_char Total number of sites.
#' @param barcodes A matrix of lineage barcodes for the leaf nodes.
#'
#' @return A matrix of barcodes for both leaf and internal nodes, where each internal node's barcode is inferred from its children.
#' @export
ancestor_inference_time_scale <- function(tree, N_char, barcodes) {
  # Identify internal nodes (parent nodes) and sort them in decreasing order.
  internal_nodes <- sort(unique(tree$edge[, 1]), decreasing = TRUE)

  # Create a character matrix for internal nodes with placeholder "0" values.
  internal_barcodes <- matrix("0", nrow = length(internal_nodes), ncol = N_char)
  rownames(internal_barcodes) <- paste("node", sort(internal_nodes, decreasing = FALSE), sep = "_")

  # Combine leaf barcodes with internal node placeholders.
  barcodes_matrix <- rbind(barcodes, internal_barcodes)

  # Iterate over each internal node to infer its barcode from its children.
  for (node in internal_nodes) {
    # Find indices of children for the current node.
    children_indices <- tree$edge[tree$edge[, 1] == node, 2]
    # For each site (column), assign the parent's barcode based on its children.
    for (i in 1:N_char) {
      chars <- barcodes_matrix[children_indices, i, drop = TRUE]
      # If all children have the same character, assign that value; otherwise, assign "0".
      if (length(unique(chars)) == 1) {
        barcodes_matrix[node, i] <- unique(chars)
      } else {
        barcodes_matrix[node, i] <- "0"
      }
    }
  }

  return(barcodes_matrix)
}


