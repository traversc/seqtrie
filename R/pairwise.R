#' @title Compute distances between all combinations of two sets of sequences
#' @description Compute distances between all combinations of query and target sequences
#' @param query `r rdoc("query")`
#' @param target `r rdoc("target")`
#' @param mode `r rdoc("mode")`
#' @param cost_matrix `r rdoc("cost_matrix")`
#' @param gap_cost `r rdoc("gap_cost")`
#' @param gap_open_cost `r rdoc("gap_open_cost")`
#' @param nthreads `r rdoc("nthreads")`
#' @param show_progress `r rdoc("show_progress")`
#' @details This function calculates all combinations of pairwise distances based on Hamming, Levenshtein or Anchored algorithms.
#' The output is a NxM matrix where N = length(query) and M = length(target).
#' Note: this can take a _really_ long time; be careful with input size.
#' 
#' `r rdoc("details")`
#' @return The output is a distance matrix between all query (rows) and target (columns) sequences.
#' For anchored searches, the output also includes attributes "query_size" and "target_size" which are matrices containing the lengths of the query and target sequences that are aligned.
#' @examples
#' dist_matrix(c("ACGT", "AAAA"), c("ACG", "ACGT"), mode = "global")
#' @name dist_matrix
dist_matrix <- function(query, target, mode, cost_matrix = NULL, gap_cost = NA_integer_, gap_open_cost = NA_integer_, nthreads = 1, show_progress = FALSE) {
  charset <- unique(c(get_charset(query), get_charset(target)))
  check_alignment_params(mode, cost_matrix, gap_cost, gap_open_cost, charset, diag_must_be_zero = FALSE)
  mode <- normalize_mode_parameter(mode)
  if (!is.na(gap_open_cost) && !is.na(gap_cost) && gap_open_cost > 0L) {
    gap_open_cost <- gap_open_cost + gap_cost
  }
  c_dist_matrix(query, target, mode, cost_matrix, gap_cost, gap_open_cost, nthreads, show_progress)
}

#' @title Pairwise distance between two sets of sequences
#' @description Compute the pairwise distance between two sets of sequences
#' @param query `r rdoc("query")`
#' @param target `r rdoc("target")`. Must be the same length as query.
#' @param mode `r rdoc("mode")`
#' @param cost_matrix `r rdoc("cost_matrix")`
#' @param gap_cost `r rdoc("gap_cost")`
#' @param gap_open_cost `r rdoc("gap_open_cost")`
#' @param nthreads `r rdoc("nthreads")`
#' @param show_progress `r rdoc("show_progress")`
#' @details
#' This function calculates pairwise distances based on Hamming, Levenshtein or Anchored algorithms. _query_ and _target_ must be the same length.
#' 
#' `r rdoc("details")`
#'
#' @return The output of this function is a vector of distances. If mode == "anchored" then the output also includes attributes "query_size" and
#' "target_size" which are vectors containing the lengths of the query and target sequences that are aligned.
#' @examples
#' dist_pairwise(c("ACGT", "AAAA"), c("ACG", "ACGT"), mode = "global")
#' @name dist_pairwise
dist_pairwise <- function(query, target, mode, cost_matrix = NULL, gap_cost = NA_integer_, gap_open_cost = NA_integer_, nthreads = 1, show_progress = FALSE) {
  if (length(query) != length(target)) {
    stop("query and target must be the same length")
  }
  charset <- unique(c(get_charset(query), get_charset(target)))
  check_alignment_params(mode, cost_matrix, gap_cost, gap_open_cost, charset, diag_must_be_zero = FALSE)
  mode <- normalize_mode_parameter(mode)
  if (!is.na(gap_open_cost) && !is.na(gap_cost) && gap_open_cost > 0L) {
    gap_open_cost <- gap_open_cost + gap_cost
  }
  c_dist_pairwise(query, target, mode, cost_matrix, gap_cost, gap_open_cost, nthreads, show_progress)
}
