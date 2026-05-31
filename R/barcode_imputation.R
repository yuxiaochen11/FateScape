#' Missing-safe barcode imputation utilities for FateScape
#'
#' This file is a drop-in replacement for FateScape/R/barcode_imputation.R.
#' It preserves the original public function names while normalizing all dropout
#' encodings to "-" before any logical comparison. This prevents errors such as
#' `if (x == "-") : missing value where TRUE/FALSE needed` when barcode entries
#' contain NA, "?", "", "NA", "NaN", or -1.

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.fs_is_missing_barcode <- function(x) {
  is.na(x) | as.character(x) %in% c("-", "?", "", "NA", "NaN", "nan", "NULL", "None", "-1")
}

.fs_normalize_barcode_matrix <- function(x, N_char = NULL, ncells = NULL) {
  rn <- rownames(x)
  cn <- colnames(x)

  x <- as.matrix(x)
  storage.mode(x) <- "character"
  x[.fs_is_missing_barcode(x)] <- "-"

  if (!is.null(N_char)) {
    N_char <- as.integer(N_char)
    if (ncol(x) < N_char) {
      stop("Barcode matrix has fewer columns than N_char: ", ncol(x), " < ", N_char)
    }
    x <- x[, seq_len(N_char), drop = FALSE]
  }

  if (!is.null(ncells)) {
    ncells <- as.integer(ncells)
    if (nrow(x) < ncells) {
      stop("Barcode matrix has fewer rows than ncells: ", nrow(x), " < ", ncells)
    }
    x <- x[seq_len(ncells), , drop = FALSE]
  }

  if (is.null(rn) || any(is.na(rn)) || any(!nzchar(rn))) {
    rownames(x) <- paste("cell", seq_len(nrow(x)), sep = "_")
  } else {
    rownames(x) <- rn[seq_len(nrow(x))]
  }

  if (!is.null(N_char)) {
    colnames(x) <- paste0("s", seq_len(ncol(x)))
  } else if (is.null(cn) || any(is.na(cn)) || any(!nzchar(cn))) {
    colnames(x) <- paste0("s", seq_len(ncol(x)))
  } else {
    colnames(x) <- cn[seq_len(ncol(x))]
  }

  x
}

.fs_distance_with_dropout <- function(barcode_1, barcode_2, N_char = NULL) {
  if (length(barcode_1) != length(barcode_2)) {
    stop("The two barcodes are not of the same length")
  }
  if (is.null(N_char)) N_char <- length(barcode_1)
  barcode_1 <- as.character(barcode_1[seq_len(N_char)])
  barcode_2 <- as.character(barcode_2[seq_len(N_char)])
  barcode_1[.fs_is_missing_barcode(barcode_1)] <- "-"
  barcode_2[.fs_is_missing_barcode(barcode_2)] <- "-"
  sum((barcode_1 != barcode_2) & (barcode_1 != "-") & (barcode_2 != "-"), na.rm = TRUE)
}

.fs_row_equal <- function(x, y) {
  x <- as.character(x)
  y <- as.character(y)
  x[.fs_is_missing_barcode(x)] <- "-"
  y[.fs_is_missing_barcode(y)] <- "-"
  length(x) == length(y) && all(x == y)
}

.fs_non_dropout <- function(x) {
  x <- as.character(x)
  x[! .fs_is_missing_barcode(x)]
}

.fs_seq <- function(from, to) {
  if (is.na(from) || is.na(to) || from > to) integer(0) else seq.int(from, to)
}

# -------------------------------------------------------------------------
# Original public functions, made missing-safe
# -------------------------------------------------------------------------

#' Group label barcodes (For imputation)
#'
#' @param label_barcodes A matrix of label barcodes.
#' @param remain_sites A matrix of remaining sites.
#' @param n If n==1, return grouped label barcodes; otherwise, return grouped remaining sites.
#'
#' @return A list of grouped barcodes.
barcodes_grouping <- function(label_barcodes, remain_sites, n) {
  label_barcodes <- .fs_normalize_barcode_matrix(label_barcodes)
  remain_sites <- .fs_normalize_barcode_matrix(remain_sites)

  # Keep row names synchronized.
  rownames(remain_sites) <- rownames(label_barcodes)

  g <- 0
  group_list <- list()
  remain_list <- list()

  while (nrow(label_barcodes) != 0) {
    g <- g + 1
    ref <- label_barcodes[1, , drop = TRUE]
    group_row <- which(apply(label_barcodes, 1, function(x) .fs_distance_with_dropout(x, ref) == 0))

    group_list[[g]] <- label_barcodes[group_row, , drop = FALSE]
    remain_list[[g]] <- remain_sites[group_row, , drop = FALSE]

    label_barcodes <- label_barcodes[-group_row, , drop = FALSE]
    remain_sites <- remain_sites[-group_row, , drop = FALSE]
  }

  if (n == 1) group_list else remain_list
}


#' Impute label barcodes
#'
#' @param label A matrix of label barcodes.
#'
#' @return Label barcodes after imputing stochastic dropout.
lable_imputation <- function(label) {
  label <- .fs_normalize_barcode_matrix(label)
  if (nrow(label) == 0 || ncol(label) == 0) return(label)

  for (i in seq_len(ncol(label))) {
    ref <- unique(.fs_non_dropout(label[, i]))
    if (length(ref) == 1) {
      label[, i] <- rep(ref, nrow(label))
    }
  }
  label
}


#' Impute remaining sites
#'
#' @param remain_s A matrix of remaining sites.
#' @param alpha Threshold for dropout ratio to decide if imputation should occur.
#'
#' @return Remaining sites after imputing stochastic dropout.
remaining_sites_imputation <- function(remain_s, alpha) {
  remain_s <- .fs_normalize_barcode_matrix(remain_s)
  if (nrow(remain_s) == 0 || ncol(remain_s) == 0) return(remain_s)

  for (i in seq_len(ncol(remain_s))) {
    missing_prop <- sum(remain_s[, i] == "-", na.rm = TRUE) / nrow(remain_s)
    if (is.na(missing_prop) || missing_prop > alpha) next

    ref <- unique(.fs_non_dropout(remain_s[, i]))
    if (length(ref) == 0 || !("-" %in% remain_s[, i])) next

    if (length(ref) == 1) {
      remain_s[, i] <- rep(ref, nrow(remain_s))
    } else if (length(ref) > 1) {
      if ("0" %in% ref) {
        remain_s[remain_s[, i] == "-", i] <- "0"
      } else {
        drop_pos <- which(remain_s[, i] == "-")
        ref_barcode <- remain_s[-drop_pos, , drop = FALSE]
        if (nrow(ref_barcode) == 0) next
        for (n_bar in drop_pos) {
          bar <- remain_s[n_bar, ]
          HD <- sapply(seq_len(nrow(ref_barcode)), function(j) .fs_distance_with_dropout(bar, ref_barcode[j, ]))
          best_idx <- which.min(HD)
          ref_site <- ref_barcode[best_idx, i]
          if (!.fs_is_missing_barcode(ref_site)) remain_s[n_bar, i] <- ref_site
        }
      }
    }
  }
  remain_s
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
  barcodes_dropout <- .fs_normalize_barcode_matrix(barcodes_dropout, N_char = N_char)
  N_char <- ncol(barcodes_dropout)
  r_n <- max(0L, min(as.integer(r_n), N_char))

  if (r_n == 0L || N_char == 0L) return(barcodes_dropout)

  order_site <- order(colSums(barcodes_dropout != "-", na.rm = TRUE), decreasing = TRUE)
  label_barcodes <- barcodes_dropout[, order_site[seq_len(r_n)], drop = FALSE]
  remain_cols <- setdiff(seq_len(N_char), order_site[seq_len(r_n)])
  remain_sites <- barcodes_dropout[, remain_cols, drop = FALSE]

  lbar_group_list <- barcodes_grouping(label_barcodes, remain_sites, 1)
  rbar_group_list <- barcodes_grouping(label_barcodes, remain_sites, 2)
  n_groups <- length(lbar_group_list)

  out_list <- vector("list", n_groups)
  for (i in seq_len(n_groups)) {
    lab <- lbar_group_list[[i]]
    rem <- rbar_group_list[[i]]
    if (nrow(lab) > 1) lab <- lable_imputation(lab)
    if (nrow(rem) > 1) rem <- remaining_sites_imputation(rem, alpha = 0.5)
    out_list[[i]] <- cbind(lab, rem)
  }

  barcodes_imputation <- do.call(rbind, out_list)
  # Restore original site order.
  barcodes_imputation <- barcodes_imputation[rownames(barcodes_dropout), colnames(barcodes_dropout), drop = FALSE]
  barcodes_imputation
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
  barcodes_dropout <- .fs_normalize_barcode_matrix(barcodes_dropout, N_char = N_char, ncells = ncells)
  N_char <- ncol(barcodes_dropout)
  ncells <- nrow(barcodes_dropout)
  r_n <- max(0L, min(as.integer(r_n), N_char))

  if (N_char == 0L || ncells == 0L || r_n == 0L) return(barcodes_dropout)

  # Order sites by number of non-dropout entries (descending).
  order_site <- order(colSums(barcodes_dropout != "-", na.rm = TRUE), decreasing = TRUE)
  ref_cols <- order_site[seq_len(r_n)]
  remain_cols <- setdiff(seq_len(N_char), ref_cols)

  first_rn_sites <- barcodes_dropout[, ref_cols, drop = FALSE]

  # Impute dropout among reference sites between compatible cell pairs.
  if (ncells >= 2L) {
    for (i in seq_len(ncells - 1L)) {
      ref <- as.character(first_rn_sites[i, ])
      ref[.fs_is_missing_barcode(ref)] <- "-"

      for (j in seq.int(i + 1L, ncells)) {
        row_j <- as.character(first_rn_sites[j, ])
        row_j[.fs_is_missing_barcode(row_j)] <- "-"

        common_idx <- which(ref != "-" & row_j != "-")
        compatible <- length(common_idx) == 0L || all(ref[common_idx] == row_j[common_idx])

        if (isTRUE(compatible)) {
          for (k in seq_len(r_n)) {
            a <- ref[k]
            b <- row_j[k]
            if (a != "-" && b == "-") {
              row_j[k] <- a
            } else if (a == "-" && b != "-") {
              ref[k] <- b
            }
          }
          first_rn_sites[j, ] <- row_j
        }
      }
      first_rn_sites[i, ] <- ref
    }
  }

  last_rn_sites <- if (length(remain_cols) > 0L) {
    barcodes_dropout[, remain_cols, drop = FALSE]
  } else {
    matrix(nrow = ncells, ncol = 0, dimnames = list(rownames(barcodes_dropout), character(0)))
  }

  all_sites <- cbind(first_rn_sites, last_rn_sites)

  # Group cells by identical reference barcode patterns.
  unique_refs <- first_rn_sites[!duplicated(as.data.frame(first_rn_sites, stringsAsFactors = FALSE)), , drop = FALSE]
  n_groups <- nrow(unique_refs)
  groups <- vector("list", n_groups)

  for (i in seq_len(n_groups)) {
    rows_match <- apply(first_rn_sites, 1, function(x) .fs_row_equal(x, unique_refs[i, ]))
    g <- all_sites[rows_match, , drop = FALSE]
    groups[[i]] <- g
  }

  for (i in seq_len(n_groups)) {
    if (nrow(groups[[i]]) > 1L && ncol(last_rn_sites) > 0L) {
      last_rn_sites_i <- last_rn_sites[rownames(groups[[i]]), , drop = FALSE]
      ref_new_n <- 0L

      # Move columns with repeated non-dropout values to the front as additional references.
      if (ncol(last_rn_sites_i) > 0L) {
        for (k in seq_len(ncol(last_rn_sites_i))) {
          vals <- .fs_non_dropout(last_rn_sites_i[, k])
          tbl <- table(vals)
          if (length(tbl) > 1L && all(tbl > 1L)) {
            ref_name <- colnames(last_rn_sites_i)[k]
            left_names <- setdiff(colnames(last_rn_sites_i), ref_name)
            last_rn_sites_i <- last_rn_sites_i[, c(ref_name, left_names), drop = FALSE]
            ref_new_n <- ref_new_n + 1L
          }
        }
      }

      if ((ref_new_n + 1L) > ncol(last_rn_sites_i)) ref_new_n <- max(0L, ref_new_n - 1L)

      if (ref_new_n > 0L) {
        ref_new_i <- last_rn_sites_i[, seq_len(ref_new_n), drop = FALSE]
        keep_ref_rows <- rowSums(ref_new_i == "-", na.rm = TRUE) != ncol(ref_new_i)
        if (!any(keep_ref_rows)) {
          ref_new_n <- 0L
        } else {
          last_new_cols <- .fs_seq(ref_new_n + 1L, ncol(last_rn_sites_i))
          last_new_i <- last_rn_sites_i[keep_ref_rows, last_new_cols, drop = FALSE]
          ref_new_i <- ref_new_i[keep_ref_rows, , drop = FALSE]

          # Pairwise imputation among additional reference sites.
          if (nrow(ref_new_i) >= 2L) {
            for (q in seq_len(nrow(ref_new_i) - 1L)) {
              ref2 <- as.character(ref_new_i[q, ])
              ref2[.fs_is_missing_barcode(ref2)] <- "-"
              for (j in seq.int(q + 1L, nrow(ref_new_i))) {
                row_j <- as.character(ref_new_i[j, ])
                row_j[.fs_is_missing_barcode(row_j)] <- "-"
                common_idx <- which(ref2 != "-" & row_j != "-")
                compatible <- length(common_idx) == 0L || all(ref2[common_idx] == row_j[common_idx])
                if (isTRUE(compatible)) {
                  for (k in seq_len(ref_new_n)) {
                    a <- ref2[k]
                    b <- row_j[k]
                    if (a != "-" && b == "-") {
                      row_j[k] <- a
                    } else if (a == "-" && b != "-") {
                      ref2[k] <- b
                    }
                  }
                  ref_new_i[j, ] <- row_j
                }
              }
              ref_new_i[q, ] <- ref2
            }
          }

          last_rn_sites_i2 <- cbind(ref_new_i, last_new_i)

          # Group by additional references and impute the remaining columns.
          if (ncol(ref_new_i) == 1L) {
            ref_key <- as.character(ref_new_i[, 1])
            unique_key <- unique(ref_key)
            groups_new <- vector("list", length(unique_key))
            for (p in seq_along(unique_key)) {
              rows_match <- ref_key == unique_key[p]
              groups_new[[p]] <- last_rn_sites_i2[rows_match, , drop = FALSE]
            }
          } else {
            unique_ref_new <- ref_new_i[!duplicated(as.data.frame(ref_new_i, stringsAsFactors = FALSE)), , drop = FALSE]
            groups_new <- vector("list", nrow(unique_ref_new))
            for (p in seq_len(nrow(unique_ref_new))) {
              rows_match <- apply(ref_new_i, 1, function(x) .fs_row_equal(x, unique_ref_new[p, ]))
              groups_new[[p]] <- last_rn_sites_i2[rows_match, , drop = FALSE]
            }
          }

          for (j in seq_along(groups_new)) {
            obj <- groups_new[[j]]
            if (ncol(obj) > ncol(ref_new_i)) {
              for (k in seq.int(ncol(ref_new_i) + 1L, ncol(obj))) {
                unique_vals <- unique(.fs_non_dropout(obj[, k]))
                miss <- obj[, k] == "-"
                if (length(unique_vals) == 1L) {
                  obj[miss, k] <- unique_vals
                } else {
                  obj[miss, k] <- "0"
                }
              }
            }
            groups_new[[j]] <- obj
          }

          l <- do.call(rbind, groups_new)
          if (!is.null(l) && nrow(l) > 0L) {
            l <- l[, colnames(last_rn_sites), drop = FALSE]
            a <- last_rn_sites[rownames(groups[[i]]), , drop = FALSE]
            missing_rows <- setdiff(rownames(a), rownames(l))
            if (length(missing_rows) > 0L) {
              last_rn_sites_i <- rbind(l, last_rn_sites[missing_rows, , drop = FALSE])
            } else {
              last_rn_sites_i <- l
            }
            last_rn_sites_i <- last_rn_sites_i[rownames(a), , drop = FALSE]
          }
        }
      }

      if (ref_new_n == 0L) {
        for (k in seq_len(ncol(last_rn_sites_i))) {
          unique_vals <- unique(.fs_non_dropout(last_rn_sites_i[, k]))
          miss <- last_rn_sites_i[, k] == "-"
          if (length(unique_vals) == 1L) {
            last_rn_sites_i[miss, k] <- unique_vals
          } else {
            last_rn_sites_i[miss, k] <- "0"
          }
        }
      }

      groups[[i]] <- cbind(groups[[i]][, seq_len(r_n), drop = FALSE], last_rn_sites_i)
    }
  }

  barcodes <- do.call(rbind, groups)
  if (is.null(barcodes) || nrow(barcodes) == 0L) {
    barcodes <- cbind(first_rn_sites, last_rn_sites)
  }

  # Restore cell and site order.
  desired_rows <- rownames(barcodes_dropout)
  desired_cols <- paste0("s", seq_len(N_char))
  barcodes <- barcodes[desired_rows, desired_cols, drop = FALSE]

  # Final guard: no R NA values should leave this function.
  barcodes[.fs_is_missing_barcode(barcodes)] <- "-"
  barcodes
}
