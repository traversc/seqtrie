# A cost matrix must
# 1) Be a square matrix _OR_ NULL
# 2) If gap_cost is NULL, cost_matrix must have a "gap" entry instead
# 3) Have all values be non-negative integers
# 4) cost_matrix (if not null) should have an entry for every element in charset (a CharacterVector of single char strings)
# Not exported, internal function only
# This function does not return anything but throws an error if the cost matrix is invalid
check_cost_matrix <- function(cost_matrix, gap_cost, charset) {
    if(is.null(cost_matrix)) return()
    if(!is.matrix(cost_matrix)) {
        stop("Cost matrix must be a matrix or NULL")
    }
    if(is.null(gap_cost)) {
        if(! "gap" %in% rownames(cost_matrix)) {
            stop("Cost matrix must have a 'gap' entry if gap_cost is NULL")
        }
    }
    if(!all(charset %in% rownames(cost_matrix))) {
        stop("Cost matrix must have an entry for every character in query or target sequences")
    }
    if(nrow(cost_matrix) != ncol(cost_matrix)) {
        stop("Cost matrix must be square")
    }
    rnames <- rownames(cost_matrix)
    cnames <- colnames(cost_matrix)
    if(is.null(rnames) || is.null(cnames)) {
        stop("Cost matrix must have row and column names")
    }
    if(!all(rnames == cnames)) {
        stop("Cost matrix must have row and column names that are the same")
    }
    if(!all(cost_matrix == as.integer(cost_matrix))) {
        stop("Cost matrix must have all integer values")
    }
    if(!all(cost_matrix >= 0)) {
        stop("Cost matrix must have all non-negative values")
    }
}

# Function for checking parameters to alignment functions
# This function does not return anything but throws an error if the parameters are invalid
# 1) The alignment "mode" is one of "hamming", "levenshtein" or "anchored"
# 2) check_cost_matrix checks cost_matrix and gap_cost (see above)
# 3) If mode == "hamming" throw a warning that cost_matrix and gap_cost is ignored IF not NULL
# This check applies to RadixTree$search, dist_pairwise, dist_matrix
# Internal function only, not exported
check_alignment_params <- function(mode, cost_matrix, gap_cost, charset) {
    mode <- tolower(mode)
    if(!mode %in% c("hamming", "levenshtein", "anchored", "hm", "lv", "an")) {
        stop("mode must be one of hamming (hm), levenshtein (lv) or anchored (an)")
    }
    check_cost_matrix(cost_matrix, gap_cost, charset)
    if( (mode %in% c("hamming", "hm")) && (!is.null(cost_matrix) || !is.null(gap_cost)) ) {
        warning("cost_matrix and gap_cost are ignored when mode is 'hamming'")
    }
}

# Function to normalize mode parameter for input into Rcpp functions
# Internal function only, not exported
normalize_mode_parameter <- function(mode) {
    mode <- tolower(mode)
    if(mode %in% c("hamming", "hm")) {
        mode <- "hamming"
    } else if(mode %in% c("levenshtein", "lv")) {
        mode <- "levenshtein"
    } else if(mode %in% c("anchored", "an")) {
        mode <- "anchored"
    }
    mode
}

# Append gap cost to cost matrix
# The cost_matrix must have a "gap" entry as input to the internal C++ functions
# Assume parameters were properly checked by check_alignment_params
# Internal function only, not exported
append_gap_cost <- function(cost_matrix, gap_cost) {
    if(is.null(cost_matrix)) {
        # do nothing
    } else if(is.null(gap_cost)) {
        # do nothing, cost_matrix already has a "gap" entry and there is no replacement
    } else if("gap" %in% rownames(cost_matrix)) {
        # replace gap entry in cost matrix with gap_cost
        cost_matrix[,"gap"] <- cost_matrix["gap",] <- gap_cost
    } else {
        # append gap_cost to cost matrix
        cost_matrix_with_gap <- matrix(gap_cost, nrow=nrow(cost_matrix)+1, ncol=nrow(cost_matrix)+1)
        colnames(cost_matrix_with_gap) <- rownames(cost_matrix_with_gap) <- c(rownames(cost_matrix), "gap")
        cost_matrix_with_gap[1:nrow(cost_matrix), 1:nrow(cost_matrix)] <- cost_matrix
        cost_matrix <- cost_matrix_with_gap
    }
    cost_matrix
}

#' @title Generate a simple cost matrix
#' @description Generate a cost matrix for use with the \code{search} method
#' @param alphabet A string representing all possible characters in both query and target sequences (e.g. "ACGT")
#' @param match The cost of a match
#' @param mismatch The cost of a mismatch
#' @param gap The cost of a gap, ignored if `include_gap` is FALSE (default: 1)
#' @param include_gap Whether to encode the cost of a gap in the matrix directly (default: FALSE)
#' @return A cost matrix
#' @examples
#' generate_cost_matrix("ACGT", match = 0, mismatch = 1)
#' @export
generate_cost_matrix <- function(alphabet, match = 0, mismatch = 1, gap = 1, include_gap = FALSE) {
    alphabet <- strsplit(alphabet, "")[[1]]
    if(include_gap) {
        alphabet <- c(alphabet, "gap")
    }
    n <- length(alphabet)
    x <- matrix(nrow = n, ncol = n)
    rownames(x) <- alphabet
    colnames(x) <- alphabet
    x[lower.tri(x)] <- mismatch
    x[upper.tri(x)] <- mismatch
    if(include_gap) {
        x[,"gap"] <- gap
        x["gap",] <- gap
    }
    diag(x) <- match
    x
}