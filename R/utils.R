# Internal function to print Roxygen documentation, since a lot of it is repeated between functions
rdoc <- function(what) {
  if (what == "details") {
    return('Three types of distance metrics are supported, based on the form of alignment performed. These are: Hamming, Global (Levenshtein) and Anchored.

An anchored alignment is a form of semi-global alignment, where the query sequence is "anchored" (global) to the beginning of both the query and target sequences,
but is semi-global in that the end of the either the query sequence or target sequence (but not both) can be unaligned. This type of alignment is sometimes called an "extension" alignment in literature.

In contrast a global alignment must align the entire query and target sequences. When mismatch and indel costs are equal to 1, this is also known as the Levenshtein distance.

By default, if mode == "global" or "anchored", all mismatches and indels are given a cost of 1. However, you can define your own distance metric by setting the substitution cost_matrix and separate gap parameters.
The cost_matrix is a strictly positive square integer matrix of substitution costs and should include all characters in query and target as column- and rownames. Any rows/columns named "gap" or "gap_open" are ignored.
To set the cost of a gap (insertion or deletion), use the gap_cost parameter (a single positive integer). To enable affine gaps, provide the gap_open_cost parameter (a single positive integer) in addition to gap_cost.
If affine alignment is used, the total cost of a gap of length L is defined as:
TOTAL_GAP_COST = gap_open_cost + (gap_cost * gap_length).

If mode == "hamming" all alignment parameters are ignored; mismatch is given a distance of 1 and gaps are not allowed.
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

  # gap parameters are provided separately in the R APIs now
  # validate gap_cost / gap_open_cost here, independent of matrix contents
  if (!is.null(gap_cost)) {
    if (!(is_integerlike(gap_cost) && length(gap_cost) == 1 && gap_cost > 0)) {
      stop("gap_cost must be a single positive integer when provided")
    }
  }
  if (!is.null(gap_open_cost)) {
    if (!(is_integerlike(gap_open_cost) && length(gap_open_cost) == 1 && gap_open_cost >= 0)) {
      stop("gap_open_cost must be a single non-negative integer when provided")
    }
    if (is.null(gap_cost)) {
      stop("If gap_open_cost is defined, gap_cost must also be defined")
    }
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
  if (!mode %in% c("hamming", "hm",
                   "global", "gb", "levenshtein", "lv",
                   "anchored", "an", "extension", "en")) {
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

#' @title Generate a simple cost matrix
#' @description Generate a cost matrix for use with the \code{search} method.
#' @param charset A string of all allowed characters in both query and target sequences (e.g. \code{"ACGT"}).
#' @param ambiguity_base A single character (e.g. \code{"N"}) that will match any character in \code{charset} at the cost of \code{match}.  Defaults to \code{NULL}.
#' @param match Integer cost of a match.
#' @param mismatch Integer cost of a mismatch.
#' @return A square cost matrix with row- and column-names given by \code{charset}, plus the optional \code{ambiguity_base}. Gap costs are no longer included here; pass \code{gap_cost} and \code{gap_open_cost} to distance/search functions.
#' @examples
#' generate_cost_matrix("ACGT", match = 0, mismatch = 1)
#' generate_cost_matrix("ACGT", ambiguity_base = "N", match = 0, mismatch = 1)
#' @export
generate_cost_matrix <- function(charset,
                                 ambiguity_base = NULL,
                                 match       = 0L,
                                 mismatch    = 1L) 
{
  # split into vector of single chars
  chars <- strsplit(charset, "")[[1]]
  
  # handle single ambiguity base
  if (!is.null(ambiguity_base)) {
    if (!is.character(ambiguity_base) || nchar(ambiguity_base) != 1L) {
      stop("`ambiguity_base` must be a single character, e.g. 'N'")
    }
    amb <- ambiguity_base
    # add to charset if not already present
    if (!amb %in% chars) {
      chars <- c(chars, amb)
    }
  }

  # check integer-like costs
  if (!is_integerlike(match) || !is_integerlike(mismatch)) {
    stop("match and mismatch must be integer-like")
  }

  # initialize matrix
  n <- length(chars)
  x <- matrix(0L, nrow = n, ncol = n, dimnames = list(chars, chars))

  # set mismatch off-diagonal, match on-diagonal
  x[lower.tri(x)] <- mismatch
  x[upper.tri(x)] <- mismatch
  diag(x)         <- match

  # override ambiguity_base interactions to always be 'match'
  if (!is.null(ambiguity_base)) {
    others <- setdiff(chars, amb)
    x[amb, others] <- match
    x[others, amb] <- match
    x[amb, amb]    <- match
  }

  x
}
