#' Group label barcodes (For imputation)
#'
#' @param label_barcodes A matrix of label barcodes.
#' @param remain_sites A matrix of remaining sites.
#' @param n If n==1, return grouped label barcodes; otherwise, return grouped remaining sites.
#'
#' @return A list of grouped barcodes.
barcodes_grouping <- function(label_barcodes, remain_sites, n) {
  g <- 0
  group_list <- list()
  remain_list <- list()

  while (nrow(label_barcodes) != 0) {
    g <- g + 1
    # Calculate the Hamming distance between the first row and each row; if the distance is 0, they belong to the same group.
    group_row <- which(apply(label_barcodes, 1, function(x) distance_with_dropout(x, label_barcodes[1,]) == 0))

    if (length(group_row) == 1) {
      group_list[[g]] <- t(as.matrix(label_barcodes[group_row, ]))
      rownames(group_list[[g]]) <- rownames(label_barcodes)[group_row]
      remain_list[[g]] <- t(as.matrix(remain_sites[group_row, ]))
      rownames(remain_list[[g]]) <- rownames(remain_sites)[group_row]
    } else {
      group_list[[g]] <- label_barcodes[group_row, ]
      remain_list[[g]] <- remain_sites[group_row, ]
    }

    label_barcodes <- label_barcodes[-group_row, , drop = FALSE]
    remain_sites <- remain_sites[-group_row, , drop = FALSE]
  }

  if (n == 1) {
    return(group_list)
  } else {
    return(remain_list)
  }
}


#' Impute label barcodes
#'
#' @param label A matrix of label barcodes.
#'
#' @return Label barcodes after imputing stochastic dropout.
lable_imputation <- function(label) {
  # For each column, if there is a unique non-dropout value, impute it across the entire column.
  for (i in 1:ncol(label)) {
    ref <- unique(label[, i])
    ref <- ref[ref != "-"]
    if (length(ref) == 0) next
    label[, i] <- rep(ref, nrow(label))
  }
  return(label)
}


#' Impute remaining sites
#'
#' @param remain_s A matrix of remaining sites.
#' @param alpha Threshold for dropout ratio to decide if imputation should occur.
#'
#' @return Remaining sites after imputing stochastic dropout.
remaining_sites_imputation <- function(remain_s, alpha) {
  for (i in 1:ncol(remain_s)) {
    missing_prop <- sum(remain_s[, i] == "-") / nrow(remain_s)
    if (missing_prop > alpha) next

    ref <- unique(remain_s[, i])
    ref <- ref[ref != "-"]

    # If there is no valid reference or no dropout exists, skip imputation.
    if (length(ref) == 0 || !("-" %in% unique(remain_s[, i]))) next

    if (length(ref) == 1) {
      remain_s[, i] <- rep(ref, nrow(remain_s))
    } else if (length(ref) > 1) {
      if ("0" %in% ref) {
        remain_s[remain_s[, i] == "-", i] <- "0"
      } else {
        drop_pos <- which(remain_s[, i] == "-")
        # Use non-dropout rows as reference barcodes.
        ref_barcode <- remain_s[-drop_pos, , drop = FALSE]
        for (n_bar in drop_pos) {
          bar <- remain_s[n_bar, ]
          # Calculate the Hamming distance between the barcode and all reference barcodes.
          HD <- sapply(1:nrow(ref_barcode), function(j) distance_with_dropout(bar, ref_barcode[j, ]))
          best_idx <- which.min(HD)
          ref_site <- ref_barcode[best_idx, i]
          remain_s[n_bar, i] <- ref_site
        }
      }
    }
  }
  return(remain_s)
}


#' Impute stochastic dropout of lineage barcodes
#'
#' @param barcodes_dropout Lineage barcodes with stochastic dropout.
#' @param N_char Total number of sites.
#' @param r_n Number of sites used as label barcodes.
#'
#' @return Lineage barcodes after imputation.
#' @export
dropout_imputation <- function(barcodes_dropout, N_char, r_n) {
  # Reorder target sites by setting column names.
  colnames(barcodes_dropout) <- paste0("s", 1:N_char)
  order_site <- order(colSums(barcodes_dropout != "-"), decreasing = TRUE)

  label_barcodes <- barcodes_dropout[, order_site[1:r_n]]
  remain_sites <- barcodes_dropout[, order_site[(r_n + 1):N_char]]

  # Group barcodes based on the label barcode patterns.
  lbar_group_list <- barcodes_grouping(label_barcodes, remain_sites, 1)
  rbar_group_list <- barcodes_grouping(label_barcodes, remain_sites, 2)
  n_groups <- length(lbar_group_list)

  # Impute label barcodes within each group (skip if the group contains only one cell).
  for (i in 1:n_groups) {
    if (nrow(lbar_group_list[[i]]) > 1) {
      lbar_group_list[[i]] <- lable_imputation(lbar_group_list[[i]])
    }
  }

  # Impute remaining sites within each group.
  for (i in 1:n_groups) {
    if (nrow(rbar_group_list[[i]]) > 1) {
      # Note: the imputation result of the remaining sites is assigned to the corresponding label group.
      lbar_group_list[[i]] <- remaining_sites_imputation(rbar_group_list[[i]])
    }
  }

  barcodes_imputation <- cbind(lbar_group_list, rbar_group_list)
  return(barcodes_imputation)
}


#' Impute stochastic dropout in lineage barcodes
#'
#' @param barcodes_dropout A matrix of lineage barcodes with dropout values.
#' @param N_char Number of characters (sites) in the barcode.
#' @param ncells Number of cells.
#' @param r_n Number of selected sites for reference.
#'
#' @return A matrix of lineage barcodes after imputation.
#' @export
dropout_imputation_alter <- function(barcodes_dropout, N_char, ncells, r_n) {
  # Set column names for all sites
  colnames(barcodes_dropout) <- paste0("s", 1:N_char)

  # Order sites by number of non-dropout entries (in descending order)
  order_site <- order(colSums(barcodes_dropout != "-"), decreasing = TRUE)

  # Select top r_n sites as reference sites
  first_rn_sites <- barcodes_dropout[, order_site[1:r_n]]

  # Impute dropout among reference sites between cell pairs
  for (i in 1:(ncells - 1)) {
    ref <- first_rn_sites[i, ]
    # Mimic original behavior: if all sites are dropout, ref becomes character(0)
    ref <- ref[sum(ref == "-") != ncol(first_rn_sites)]

    for (j in (i + 1):ncells) {
      # Identify indices where both cells have non-dropout values
      common_idx <- which(ref != "-" & first_rn_sites[j, ] != "-")
      # Condition equivalent to (s+1)/(t+1)==1 in original code:
      # If no common non-dropout values or all such positions match exactly
      if (length(common_idx) == 0 || all(ref[common_idx] == first_rn_sites[j, common_idx])) {
        for (k in 1:r_n) {
          if (ref[k] != "-" && first_rn_sites[j, k] == "-") {
            first_rn_sites[j, k] <- ref[k]
          }
          if (ref[k] == "-" && first_rn_sites[j, k] != "-") {
            ref[k] <- first_rn_sites[j, k]
          }
        }
      }
    }
    first_rn_sites[i, ] <- ref
  }

  # Combine the imputed reference sites with the remaining sites
  last_rn_sites <- barcodes_dropout[, -order_site[1:r_n]]
  all_sites <- cbind(first_rn_sites, last_rn_sites)

  # Group cells by identical reference barcode patterns
  unique_refs <- first_rn_sites[!duplicated(first_rn_sites), , drop = FALSE]
  n_groups <- nrow(unique_refs)
  groups <- list()
  for (i in 1:n_groups) {
    rows_match <- apply(first_rn_sites, 1, function(x) all(x == unique_refs[i, ]))
    g <- all_sites[rows_match, , drop = FALSE]
    if (nrow(g) < 1) {
      g <- t(as.matrix(g))
      rownames(g) <- paste("cell", which(rows_match), sep = "_")
    }
    groups[[i]] <- g
  }
  for (i in 1:n_groups) {
    if (nrow(groups[[i]]) > 1) {
      last_rn_sites_i <- last_rn_sites[rownames(groups[[i]]), , drop = FALSE]
      ref_new_n <- 0

      for (k in 1:ncol(last_rn_sites_i)) {
        non_drop_vals <- last_rn_sites_i[last_rn_sites_i[, k] != "-", k]
        tbl <- table(non_drop_vals)
        if (length(tbl) > 1 && all(tbl > 1)) {
          ref_name <- colnames(last_rn_sites_i)[k]
          left_names <- colnames(last_rn_sites_i)[-k]

          last_rn_sites_i <- last_rn_sites_i[, c(ref_name, left_names), drop = FALSE]
          ref_new_n <- ref_new_n + 1
        }
      }
      if ((ref_new_n + 1) > ncol(last_rn_sites_i)) {
        ref_new_n <- ref_new_n - 1
      }

      if (ref_new_n > 0) {
        ref_new_i <- last_rn_sites_i[, 1:ref_new_n, drop = FALSE]
        if (is.null(dim(ref_new_i))) {
          ref_new_i <- as.matrix(ref_new_i)
          colnames(ref_new_i) <- ref_name
        }

        last_new_i <- last_rn_sites_i[rowSums(ref_new_i == "-") != ncol(ref_new_i), (ref_new_n + 1):ncol(last_rn_sites_i), drop = FALSE]
        ref_new_i <- ref_new_i[rowSums(ref_new_i == "-") != ncol(ref_new_i), , drop = FALSE]
        if (is.null(dim(ref_new_i))) {
          ref_new_i <- as.matrix(ref_new_i)
          colnames(ref_new_i) <- ref_name
        }

        for (q in 1:(nrow(ref_new_i) - 1)) {
          ref2 <- ref_new_i[q, ]
          for (j in (q + 1):nrow(ref_new_i)) {
            common_idx <- which(ref2 != "-" & ref_new_i[j, ] != "-")
            if (length(common_idx) == 0 || all(ref2[common_idx] == ref_new_i[j, common_idx])) {
              for (k in 1:ref_new_n) {
                if (ref2[k] != "-" && ref_new_i[j, k] == "-") {
                  ref_new_i[j, k] <- ref2[k]
                }
                if (ref2[k] == "-" && ref_new_i[j, k] != "-") {
                  ref2[k] <- ref_new_i[j, k]
                }
              }
            }
          }
          ref_new_i[q, ] <- ref2
        }
        if (is.null(dim(last_new_i))) {
          last_new_i <- as.matrix(last_new_i)
          colnames(last_new_i) <- setdiff(colnames(last_rn_sites_i), colnames(ref_new_i))
        }
        last_rn_sites_i <- cbind(ref_new_i, last_new_i)

        if (ncol(ref_new_i) == 1) {
          n_groups_new <- length(ref_new_i[!duplicated(ref_new_i)])
          groups_new <- list()
          for (p in 1:n_groups_new) {
            rows_match <- apply(ref_new_i, 1, function(x) all(x == ref_new_i[!duplicated(ref_new_i), ][p]))
            g <- last_rn_sites_i[rows_match, , drop = FALSE]
            if (nrow(g) < 1) {
              g <- t(as.matrix(g))
              rownames(g) <- paste("cell", which(rows_match), sep = "_")
            }
            groups_new[[p]] <- g
          }
        } else {
          unique_ref_new <- ref_new_i[!duplicated(ref_new_i), , drop = FALSE]
          n_groups_new <- nrow(unique_ref_new)
          groups_new <- list()
          for (p in 1:n_groups_new) {
            rows_match <- apply(ref_new_i, 1, function(x) all(x == unique_ref_new[p, ]))
            g <- last_rn_sites_i[rows_match, , drop = FALSE]
            if (nrow(g) < 1) {
              g <- t(as.matrix(g))
              rownames(g) <- paste("cell", which(rows_match), sep = "_")
            }
            groups_new[[p]] <- g
          }
        }

        for (j in seq_along(groups_new)) {
          obj <- groups_new[[j]]
          for (k in (ncol(ref_new_i) + 1):ncol(obj)) {
            unique_vals <- unique(obj[obj[, k] != "-", k])
            if (length(unique_vals) == 1) {
              obj[obj[, k] == "-", k] <- unique_vals
            } else {
              obj[obj[, k] == "-", k] <- 0
            }
          }
          groups_new[[j]] <- obj
        }
        l <- do.call(rbind, groups_new)
        l <- l[, colnames(last_rn_sites), drop = FALSE]
        a <- last_rn_sites[rownames(groups[[i]]), , drop = FALSE]
        b <- l
        f <- which(!(rownames(a) %in% rownames(b)))
        name_f <- rownames(a)[f]
        last_rn_sites_i <- rbind(l, last_rn_sites[rownames(groups[[i]]), , drop = FALSE][f, , drop = FALSE])
        rownames(last_rn_sites_i) <- c(rownames(b), name_f)
        aa <- rownames(a)
        bb <- rownames(last_rn_sites_i)
        last_rn_sites_i <- last_rn_sites_i[match(aa, bb), , drop = FALSE]
      } else if (ref_new_n == 0) {
        for (k in 1:ncol(last_rn_sites_i)) {
          unique_vals <- unique(last_rn_sites_i[last_rn_sites_i[, k] != "-", k])
          if (length(unique_vals) == 1) {
            last_rn_sites_i[last_rn_sites_i[, k] == "-", k] <- unique_vals
          } else {
            last_rn_sites_i[last_rn_sites_i[, k] == "-", k] <- 0
          }
        }
      }
      groups[[i]] <- cbind(groups[[i]][, 1:r_n, drop = FALSE], last_rn_sites_i)
    }
  }
  barcodes <- do.call(rbind, groups)
  barcodes <- barcodes[paste("cell", 1:ncells, sep = "_"), paste0("s", 1:N_char)]
  return(barcodes)
}
