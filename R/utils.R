# Internal function to print Roxygen documentation, since a lot of it is repeated between functions
rdoc <- function(what) {
  if (what == "details") {
    return('Three types of distance metrics are supported, based on the form of alignment performed. These are: Hamming, Global (Levenshtein) and Anchored.

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
  } else if (what == "query") {
    return("A character vector of query sequences.")
  } else if (what == "sequences") {
    return("A character vector of sequences.")
  } else if (what == "target") {
    return("A character vector of target sequences.")
  } else if (what == "mode") {
    return("The distance metric to use. One of hamming (hm), global (gb) or anchored (an).")
  } else if (what == "cost_matrix") {
    return('A custom cost matrix for use with the "global" or "anchored" distance metrics. See details.')
  } else if (what == "gap_cost") {
    return('The cost of a gap for use with the "global" or "anchored" distance metrics. See details.')
  } else if (what == "gap_open_cost") {
    return("The cost of a gap opening. See details.")
  } else if (what == "nthreads") {
    return("The number of threads to use for parallel computation.")
  } else if (what == "show_progress") {
    return("Whether to show a progress bar.")
  } else if (what == "max_distance") {
    return("how far to search in units of absolute distance. Can be a single value or a vector. Mutually exclusive with max_fraction.")
  } else if (what == "max_fraction") {
    return("how far to search in units of relative distance to each query sequence length. Can be a single value or a vector. Mutually exclusive with max_distance.")
  }
}

# Internal function for testing if a numeric or integer object of arbitrary shape contains integer-like values
is_integerlike <- function(x) {
  if (is.integer(x)) {
    return(TRUE)
  }
  if (!is.numeric(x)) {
    return(FALSE)
  }
  same_val <- as.integer(x) == x
  is_na <- is.na(x)
  if (!all(same_val | is_na)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

# A cost matrix must
# 1) Be a square matrix _OR_ NULL
# 2) If gap_cost is NULL, cost_matrix must have a "gap" entry instead
# 3) Have all values be non-negative integers
# 4) cost_matrix (if not null) should have an entry for every element in charset (a CharacterVector of single char strings)
# Handling of gap parameters:
# IF cost_matrix is NULL, gap parameters are ignored
# IF cost_matrix is defined -> either: gap_cost is defined (a single positive integer) OR cost_matrix must have a "gap" entry instead
# IF cost_matrix is defined -> any of: gap_open_cost is defined (a single integer) OR cost_matrix must have a "gap_open" entry OR neither (implies non-affine alignment)
# IF gap_open_cost is defined or in cost_matrix, gap_cost must also be defined or in cost_matrix
# Not exported, internal function only
# This function does not return anything but throws an error if invalid parameters are passed
check_cost_matrix <- function(cost_matrix, gap_cost, gap_open_cost, charset, diag_must_be_zero) {
  if (is.null(cost_matrix)) {
    return()
  }
  if (!is.matrix(cost_matrix)) {
    stop("Cost matrix must be a matrix or NULL")
  }
  if (!all(charset %in% rownames(cost_matrix))) {
    stop("Cost matrix must have an entry for every character in query or target sequences")
  }
  if (nrow(cost_matrix) != ncol(cost_matrix)) {
    stop("Cost matrix must be square")
  }
  rnames <- rownames(cost_matrix)
  cnames <- colnames(cost_matrix)
  if (is.null(rnames) || is.null(cnames)) {
    stop("Cost matrix must have row and column names")
  }
  if (!all(rnames == cnames)) {
    stop("Cost matrix must have row and column names that are the same")
  }
  if (!is_integerlike(cost_matrix)) {
    stop("Cost matrix must have all integer values")
  }
  if (!all(cost_matrix >= 0)) {
    stop("Cost matrix must have all non-negative values")
  }
  if (diag_must_be_zero && any(diag(cost_matrix) != 0)) {
    stop("Cost matrix must have zeros on the diagonal")
  }

  # gap params, cost_matrix is defined
  gap_is_valid <- is_integerlike(gap_cost) && length(gap_cost) == 1 && gap_cost > 0
  cost_matrix_has_gap <- "gap" %in% rnames
  if (!gap_is_valid && !cost_matrix_has_gap) {
    stop("Cost matrix must have a 'gap' entry or gap_cost must be a single positive integer")
  }
  gap_open_is_null <- is.null(gap_open_cost)
  gap_open_is_valid <- is_integerlike(gap_open_cost) && length(gap_open_cost) == 1 && gap_open_cost > 0
  cost_matrix_has_gap_open <- "gap_open" %in% rnames
  if (!gap_open_is_null && !gap_open_is_valid && !cost_matrix_has_gap_open) {
    stop("If gap_open_cost is defined, it must be a single positive integer")
  }

  if ((gap_open_is_valid || cost_matrix_has_gap_open) && !(gap_is_valid || cost_matrix_has_gap)) {
    stop("If gap_open is defined, gap must also be defined")
  }
}

# Function for checking parameters to alignment functions
# This function does not return anything but throws an error if the parameters are invalid
# 1) The alignment "mode" is one of "hamming", "global" or "anchored" (or their abbreviations/synonyms)
# 2) check_cost_matrix checks cost_matrix and gap_cost (see above)
# 3) If mode == "hamming" throw a warning that cost_matrix and gap_cost is ignored IF not NULL
# This check applies to RadixTree$search, dist_pairwise, dist_matrix
# Internal function only, not exported
check_alignment_params <- function(mode, cost_matrix, gap_cost, gap_open_cost, charset, diag_must_be_zero) {
  mode <- tolower(mode)
  if (!mode %in% c("hamming", "global", "anchored", "gb", "lv", "an", "levenshtein", "lv")) {
    stop("mode must be one of hamming (hm), global (gb, lv, levenshtein) or anchored (an, en, extension)")
  }
  check_cost_matrix(cost_matrix, gap_cost, gap_open_cost, charset, diag_must_be_zero)
  if ((mode %in% c("hamming", "hm")) && (!is.null(cost_matrix) || !is.null(gap_cost) || !is.null(gap_open_cost))) {
    warning("cost_matrix and gap parameters are ignored when mode is 'hamming'")
  }
}

# Function to normalize mode parameter for input into Rcpp functions
# Internal function only, not exported
normalize_mode_parameter <- function(mode) {
  mode <- tolower(mode)
  if (mode %in% c("hamming", "hm")) {
    mode <- "hamming"
  } else if (mode %in% c("gb", "global", "levenshtein", "lv")) {
    mode <- "global"
  } else if (mode %in% c("anchored", "an", "extension", "en")) {
    mode <- "anchored"
  }
  mode
}

# This function finalizes the cost matrix for input into the C++ functions
# The input to the C++ function is a single cost_matrix parameter and does not include gap_cost and gap_open_cost parameters
# The finalized cost matrix can be any of the following:
# 1) NULL -- implies linear alignment (mismatches = 1, gaps = 1)
# 2) A square matrix with "gap" column -- implies linear alignment with a custom cost for mismatches and gaps
# 3) A square matrix with "gap" and "gap_open" columns -- implies affine alignment with a custom cost for mismatches, gaps and gap openings
# In this third case, the cost of a gap is added to the gap_open cost to get the total cost of the initial gap.
# TOTAL_GAP_COST = gap_open_cost + (gap_cost * [gap_length-1])
# Note this is different than how it is defined for the R interface, which is more intuitive
# Assume parameters were properly checked by check_alignment_params
# Internal function only, not exported
finalize_cost_matrix <- function(cost_matrix, gap_cost, gap_open_cost) {
  if (is.null(cost_matrix)) {
    return(NULL)
  }

  mode(cost_matrix) <- "integer"
  # Add gap to finalized matrix
  if (is.null(gap_cost)) {
    # do nothing, cost_matrix must already has a "gap" entry
  } else if ("gap" %in% rownames(cost_matrix)) {
    # replace gap entry in cost matrix with gap_cost
    cost_matrix[, "gap"] <- cost_matrix["gap", ] <- gap_cost
  } else {
    # append gap_cost to cost matrix
    cost_matrix_with_gap <- matrix(gap_cost, nrow = nrow(cost_matrix) + 1, ncol = nrow(cost_matrix) + 1)
    colnames(cost_matrix_with_gap) <- rownames(cost_matrix_with_gap) <- c(rownames(cost_matrix), "gap")
    cost_matrix_with_gap[seq_len(nrow(cost_matrix)), seq_len(nrow(cost_matrix))] <- cost_matrix
    cost_matrix <- cost_matrix_with_gap
  }

  # Add gap_open to finalized matrix
  if (is.null(gap_open_cost)) {
    # do nothing, gap_open is optional
  } else if ("gap_open" %in% rownames(cost_matrix)) {
    # replace gap_open entry in cost matrix with gap_open_cost
    cost_matrix[, "gap_open"] <- cost_matrix["gap_open", ] <- gap_open_cost
  } else {
    # append gap_open_cost to cost matrix
    cost_matrix_with_gap_open <- matrix(gap_open_cost, nrow = nrow(cost_matrix) + 1, ncol = nrow(cost_matrix) + 1)
    colnames(cost_matrix_with_gap_open) <- rownames(cost_matrix_with_gap_open) <- c(rownames(cost_matrix), "gap_open")
    cost_matrix_with_gap_open[seq_len(nrow(cost_matrix)), seq_len(nrow(cost_matrix))] <- cost_matrix
    cost_matrix <- cost_matrix_with_gap_open
  }

  # Add the cost of gap to gap_open in finalized matrix
  if ("gap_open" %in% rownames(cost_matrix)) {
    cost_matrix["gap_open", ] <- cost_matrix["gap_open", ] + cost_matrix["gap", ]
    cost_matrix[, "gap_open"] <- cost_matrix[, "gap_open"] + cost_matrix[, "gap"]
  }

  # Set gap to gap entries as NA_integer_ since they can never be used in an alignment
  cost_matrix["gap", "gap"] <- NA_integer_
  if ("gap_open" %in% rownames(cost_matrix)) {
    cost_matrix["gap_open", "gap"] <- NA_integer_
    cost_matrix["gap", "gap_open"] <- NA_integer_
    cost_matrix["gap_open", "gap_open"] <- NA_integer_
  }
  cost_matrix
}

# This function returns either "linear" or "affine" depending on the finalized cost matrix
# Internal function only, not exported
get_gap_type <- function(finalized_cost_matrix) {
  if (is.null(finalized_cost_matrix)) {
    return("linear")
  } else if ("gap_open" %in% rownames(finalized_cost_matrix) && "gap" %in% rownames(finalized_cost_matrix)) {
    return("affine")
  } else {
    return("linear")
  }
}

#' @title Generate a simple cost matrix
#' @description Generate a cost matrix for use with the \code{search} method
#' @param charset A string representing all possible characters in both query and target sequences (e.g. "ACGT")
#' @param match The cost of a match
#' @param mismatch The cost of a mismatch
#' @param gap The cost of a gap or NULL if this parameter will be set later.
#' @param gap_open The cost of a gap opening or NULL. If this parameter is set, gap must also be set.
#' @return A cost matrix
#' @examples
#' generate_cost_matrix("ACGT", match = 0, mismatch = 1)
#' @export
generate_cost_matrix <- function(charset, match = 0L, mismatch = 1L, gap = NULL, gap_open = NULL) {
  charset <- strsplit(charset, "")[[1]]
  gap_is_defined <- !is.null(gap)
  gap_open_is_defined <- !is.null(gap_open)

  if (!is_integerlike(match) || !is_integerlike(mismatch)) {
    stop("Cost parameters must have integer values")
  }
  if (gap_is_defined && !is_integerlike(gap)) {
    stop("Cost parameters must have integer values")
  }
  if (gap_open_is_defined && !is_integerlike(gap_open)) {
    stop("Cost parameters must have integer values")
  }

  if (gap_open_is_defined && gap_is_defined) {
    charset <- c(charset, "gap", "gap_open")
  } else if (gap_is_defined) {
    charset <- c(charset, "gap")
  } else if (gap_open_is_defined) {
    stop("If gap_open is defined, gap must also be defined")
  }
  n <- length(charset)
  x <- matrix(nrow = n, ncol = n)
  rownames(x) <- charset
  colnames(x) <- charset
  x[lower.tri(x)] <- mismatch
  x[upper.tri(x)] <- mismatch
  diag(x) <- match
  if (gap_open_is_defined) {
    x["gap", ] <- gap
    x[, "gap"] <- gap
    x["gap_open", ] <- gap_open
    x[, "gap_open"] <- gap_open
    x["gap", "gap"] <- 0L
    x["gap_open", "gap"] <- 0L
    x["gap", "gap_open"] <- 0L
    x["gap_open", "gap_open"] <- 0L
  } else if (gap_is_defined) {
    x["gap", ] <- gap
    x[, "gap"] <- gap
    x["gap", "gap"] <- 0L
  }
  x
}
