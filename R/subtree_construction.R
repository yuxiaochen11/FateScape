#' Generate a tree topology using Ward's clustering method
#'
#' @param distance_matrix A Hamming distance matrix.
#'
#' @return A phylogenetic tree obtained by clustering using Ward's method, with all edge lengths set to 1.
#' @export
 nodes_clustering<- function(distance_matrix) {
  # Perform hierarchical clustering using Ward's method.
  tree_WD <- hclust(distance_matrix, method = "ward.D")

  # Convert the clustering object to a phylogenetic tree.
  tree_WD <- as.phylo(tree_WD)

  # Set all edge lengths to 1.
  tree_WD$edge.length <- rep(1, length(tree_WD$edge.length))

  return(tree_WD)
 }


 #' Reconstruct initial sub-cell division trees based solely on lineage barcodes
 #'
 #' @param state_lineages A list of state lineages.
 #' @param barcodes_lineages A list of barcode matrices corresponding to each state lineage.
 #'
 #' @return A list of initial sub-cell division trees reconstructed only from lineage barcodes.
 #' @export
 initial_tree_construction <- function(state_lineages, barcodes_lineages) {
   trees_initial <- list()

   for (j in seq_along(state_lineages)) {
     lineage_label <- paste0("L", j)
     barcode_lineage <- barcodes_lineages[[lineage_label]]
     dist_matrix <- hamming_distance(barcode_lineage)
     tree_WD <-  nodes_clustering(dist_matrix)
     trees_initial[[lineage_label]] <- tree_WD
   }

   return(trees_initial)
 }









