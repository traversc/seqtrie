# Internal function for print Roxygen documentation, since a lot of it is repeated between functions
rdoc <- function(what) {
    if(what == "details") {
        cat('Three types of distance metrics are supported, based on the form of alignment performed. These are: Hamming, Global (Levenshtein) and Anchored.

An anchored alignment is a form of semi-global alignment, where the query sequence is "anchored" (global) to the beginning of both the query and target sequences, 
but is semi-global in that the end of the either the query sequence or target sequence (but not both) can be unaligned. This type of alignment is sometimes called an "extension" alignment in literature. 

In contrast a global alignment must align the entire query and target sequences. When mismatch and indel costs are equal to 1, this is also known as the Levenshtein distance. 

By default, if mode == "global" or "anchored", all mismatches and indels are given a cost of 1. However, you can define your own distance metric by setting the cost_matrix and gap parameters. 
The cost_matrix is a strictly positive square integer matrix and should include all characters in query and target as column- and rownames. 
To set the cost of a gap (insertion or deletion) you can include a row and column named "gap" in the cost_matrix _OR_ set the gap_cost parameter (a single positive integer). 
Similarly, the affine gap alignment can be set by including a row and column named "gap_open" in the cost_matrix _OR_ setting the gap_open_cost parameter (a single positive integer).
If affine alignment is used, the cost of a gap is defined as:
TOTAL_GAP_COST = gap_open_cost + (gap_cost * gap_length).

If mode == "hamming" all alignment parameters are ignored; z mismatch is given a distance of 1 and gaps are not allowed.
')
    } else if(what == "query") {
        cat('A character vector of query sequences.')
    } else if(what == "sequences") {
        cat('A character vector of sequences.')
    } else if(what == "target") {
        cat('A character vector of target sequences.')
    } else if(what == "mode") {
        cat('The distance metric to use. One of hamming (hm), global (gb) or anchored (an).')
    } else if(what == "cost_matrix") {
        cat('A custom cost matrix for use with the "global" or "anchored" distance metrics. See details.')
    } else if(what == "gap_cost") {
        cat('The cost of a gap for use with the "global" or "anchored" distance metrics. See details.')
    } else if(what == "gap_open_cost") {
        cat('The cost of a gap opening. See details.')
    } else if(what == "nthreads") {
        cat('The number of threads to use for parallel computation.')
    } else if(what == "show_progress") {
        cat('Whether to show a progress bar.')
    }
}

#' @title Compute distances between all combinations of two sets of sequences
#' @description Compute distances between all combinations of query and target sequences
#' @param query `r rdoc("query")`
#' @param target `r rdoc("target")`
#' @param mode `r rdoc("mode")`
#' @param cost_matrix `r rdoc("cost_matrix")`
#' @param gap_cost `r rdoc("gap_cost")`
#' @param nthreads `r rdoc("nthreads")`
#' @param show_progress `r rdoc("show_progress")`
#' @details This function calculates all combinations of pairwise distances based on Hamming, Levenshtein or Anchored algorithms. The output is a NxM matrix where N = length(query) and M = length(target).  
#' `r rdoc("details")`
#' @return The output is a distance matrix between all query (rows) and target (columns) sequences.
#' For anchored searches, the output also includes attributes "query_size" and "target_size" which are matrices containing the lengths of the query and target sequences that are aligned. 
#' @examples
#' dist_matrix(c("ACGT", "AAAA"), c("ACG", "ACGT"), mode = "global")
#' @name dist_matrix
dist_matrix <- function(query, target, mode, cost_matrix = NULL, gap_cost = NULL, gap_open_cost = NULL, nthreads = 1, show_progress = FALSE) {
    charset <- unique(c(get_charset(query), get_charset(target)))
    check_alignment_params(mode, cost_matrix, gap_cost, gap_open_cost, charset, diag_must_be_zero = FALSE)
    finalized_cost_matrix <- finalize_cost_matrix(cost_matrix, gap_cost, gap_open_cost)
    mode <- normalize_mode_parameter(mode)
    gap_type <- get_gap_type(finalized_cost_matrix)
    c_dist_matrix(query, target, mode, gap_type, finalized_cost_matrix, nthreads, show_progress)
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
#' `r rdoc("details")`
#'
#' @return The output of this function is a vector of distances. If mode == "anchored" then the output also includes attributes "query_size" and 
#' "target_size" which are vectors containing the lengths of the query and target sequences that are aligned.
#' @examples
#' dist_pairwise(c("ACGT", "AAAA"), c("ACG", "ACGT"), mode = "global")
#' @name dist_pairwise
dist_pairwise <- function(query, target, mode, cost_matrix = NULL, gap_cost = NULL, gap_open_cost = NULL, nthreads = 1, show_progress = FALSE) {
    if(length(query) != length(target)) {
        stop("query and target must be the same length")
    }
    charset <- unique(c(get_charset(query), get_charset(target)))
    check_alignment_params(mode, cost_matrix, gap_cost, gap_open_cost, charset, diag_must_be_zero = FALSE)
    finalized_cost_matrix <- finalize_cost_matrix(cost_matrix, gap_cost, gap_open_cost)
    mode <- normalize_mode_parameter(mode)
    gap_type <- get_gap_type(finalized_cost_matrix)
    c_dist_pairwise(query, target, mode, gap_type, finalized_cost_matrix, nthreads, show_progress)
}
