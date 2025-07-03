#' Stochastic dropout generation
#'
#' @description Generate stochastic dropouts on lineage barcodes.
#'
#' @param barcodes_initial A matrix of barcodes containing only inheritable dropout.
#' @param p_ds Stochastic dropout rate.
#'
#' @return The initial barcodes with stochastic dropouts applied.
#' @export
stochastic_dropout <- function(barcodes_initial, p_ds) {
  # Determine the number of cells and number of sites from the input matrix
  ncells <- nrow(barcodes_initial)
  N_char <- ncol(barcodes_initial)

  # Loop over each cell and each site
  for (i in 1:ncells) {
    for (j in 1:N_char) {
      # Randomly decide whether to drop out at the current site using the dropout probability
      f <- sample(c(0, 1), 1, prob = c(p_ds, 1 - p_ds))
      if (f == 0) {
        barcodes_initial[i, j] <- "-"
      }
    }
  }

  return(barcodes_initial)
}
