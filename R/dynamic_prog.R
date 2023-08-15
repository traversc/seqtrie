#' @title Compute distances between all combinations of two sets of sequences
#' @description Compute distances between all combinations of query and target sequences
#' @param query A character vector of query sequences
#' @param target A character vector of target sequences
#' @param mode The distance metric to use. One of hamming (hm), levenshtein (lv) or anchored (an)
#' @param cost_matrix A custom cost matrix for use with the "levenshtein" or "anchored" distance metrics
#' @param gap_cost The cost of a gap for use with the "levenshtein" or "anchored" distance metrics
#' @param nthreads The number of threads to use for parallel computation
#' @param show_progress Whether to show a progress bar
#' @return The output of this function is a distance matrix between all query (rows) and target (columns) sequences.
#' For anchored searches, the output also includes attributes "query_size" and "target_size" which are matrices containing the lengths of the query and target sequences that are aligned. 
#' @examples
#' dist_matrix(c("ACGT", "AAAA"), c("ACG", "ACGT"), mode = "levenshtein")
#' @name dist_matrix
dist_matrix <- function(query, target, mode, cost_matrix = NULL, gap_cost = NULL, nthreads = 1, show_progress = FALSE) {
    charset <- unique(c(get_charset(query), get_charset(target)))
    check_alignment_params(mode, cost_matrix, gap_cost, charset)
    cost_matrix <- append_gap_cost(cost_matrix, gap_cost)
    mode <- normalize_mode_parameter(mode)
    c_dist_matrix(query, target, mode, cost_matrix, nthreads, show_progress)
}

#' @title Pairwise distance between two sets of sequences
#' @description Compute the pairwise distance between two sets of sequences
#' @param query A character vector of query sequences
#' @param target A character vector of target sequences
#' @param mode The distance metric to use. One of hamming (hm), levenshtein (lv) or anchored (an)
#' @param cost_matrix A custom cost matrix for use with the "levenshtein" or "anchored" distance metrics
#' @param gap_cost The cost of a gap for use with the "levenshtein" or "anchored" distance metrics
#' @param nthreads The number of threads to use for parallel computation
#' @param show_progress Whether to show a progress bar
#' @details
#' The "hamming" distance is the number of positions at which the corresponding symbols are different.
#' @return The output of this function is a data.frame of all matches with columns "query" (the sequences input to the search function), 
#' "target" (the sequences inserted into the tree) and "distance" the absolute distance between query and target sequences. 
#' For anchored searches, the output also includes "query_size" and "target_size" which are the lengths of the query and target sequences that are aligned. 
#' @examples
#' dist_pairwise(c("ACGT", "AAAA"), c("ACG", "ACGT"), mode = "levenshtein")
#' @name dist_pairwise
dist_pairwise <- function(query, target, mode, cost_matrix = NULL, gap_cost = NULL, nthreads = 1, show_progress = FALSE) {
    charset <- unique(c(get_charset(query), get_charset(target)))
    check_alignment_params(mode, cost_matrix, gap_cost, charset)
    cost_matrix <- append_gap_cost(cost_matrix, gap_cost)
    mode <- normalize_mode_parameter(mode)
    c_dist_pairwise(query, target, mode, cost_matrix, nthreads, show_progress)
}
