#' State lineage information
#'
#' @param state_lineages A named list of state lineage paths (e.g., "L1", "L2", ...).
#' @param ncells Number of cells.
#' @param state A data frame containing cell states (with columns "cluster" and "cell_id").
#' @param barcodes A matrix of lineage barcodes for each cell.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{state_lineage}{A list indicating which state lineage(s) each cell belongs to.}
#'   \item{cell_lineages}{A list mapping each state lineage to the corresponding cell IDs.}
#'   \item{barcodes_lineages}{A list of barcode matrices for cells in each state lineage.}
#'   \item{state_labels_lineages}{A list of state labels (clusters) for cells in each state lineage.}
#' }
#' @export
state_lineage_info <- function(state_lineages, ncells, state, barcodes) {
  # Determine the state lineage for each observed cell.
  state_lineage <- vector("list", ncells)
  for (i in 1:ncells) {
    cell_lineages <- character(0)
    for (j in seq_along(state_lineages)) {
      lineage_label <- paste0("L", j)
      # If the cell's cluster is in the current state lineage, add the label.
      if (state$cluster[i] %in% state_lineages[[lineage_label]]) {
        cell_lineages <- c(cell_lineages, lineage_label)
      }
    }
    state_lineage[[i]] <- cell_lineages
  }
  names(state_lineage) <- state$cell_id

  # Determine which cells belong to each state lineage.
  cell_lineages_list <- list()
  for (j in seq_along(state_lineages)) {
    lineage_label <- paste0("L", j)
    # Create a logical vector indicating whether each cell belongs to the current lineage.
    cells_in_lineage <- sapply(state_lineage, function(x) lineage_label %in% x)
    cell_key <- paste0("cell_", lineage_label)
    cell_lineages_list[[cell_key]] <- names(state_lineage)[cells_in_lineage]
  }

  # Extract the barcode matrices for cells in each state lineage.
  barcodes_lineages <- list()
  for (j in seq_along(state_lineages)) {
    lineage_label <- paste0("L", j)
    cell_key <- paste0("cell_", lineage_label)
    barcodes_lineages[[lineage_label]] <- barcodes[cell_lineages_list[[cell_key]], ]
  }

  # Extract the state labels for cells in each state lineage.
  state_labels_lineages <- list()
  for (j in seq_along(state_lineages)) {
    lineage_label <- paste0("L", j)
    cell_key <- paste0("cell_", lineage_label)
    state_labels_lineages[[lineage_label]] <- state$cluster[state$cell_id %in% cell_lineages_list[[cell_key]]]
  }

  return(list(state_lineage = state_lineage,
              cell_lineages = cell_lineages_list,
              barcodes_lineages = barcodes_lineages,
              state_labels_lineages = state_labels_lineages))
}
