# StarTree helpers --------------------------------------------------------------

seqtrie_check_startree_sequences <- function(sequences, name, allow_na = FALSE) {
  check_character_vector(sequences, name)
  if (!allow_na && any(is.na(sequences))) {
    stop(name, " must not contain missing values")
  }
  invisible(sequences)
}

seqtrie_check_startree_integer <- function(x, name, positive = FALSE) {
  if (!is_integerlike(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be a single integer")
  }
  x <- as.integer(x)
  if (positive && x <= 0L) {
    stop(name, " must be a single positive integer")
  }
  if (!positive && x < 0L) {
    stop(name, " must be a single non-negative integer")
  }
  x
}

seqtrie_check_startree_params <- function(max_distance, mismatch_cost, gap_cost) {
  max_distance <- seqtrie_check_startree_integer(max_distance, "max_distance")
  mismatch_cost <- seqtrie_check_startree_integer(mismatch_cost, "mismatch_cost", positive = TRUE)
  gap_cost <- seqtrie_check_startree_integer(gap_cost, "gap_cost", positive = TRUE)
  if (max_distance %/% min(mismatch_cost, gap_cost) > 8L) {
    stop("StarTree requires max_distance / min(mismatch_cost, gap_cost) <= 8")
  }
  list(max_distance = max_distance, mismatch_cost = mismatch_cost, gap_cost = gap_cost)
}

seqtrie_check_startree_mode <- function(mode) {
  mode <- normalize_mode_parameter(mode)
  if (mode != "global") {
    stop("StarTree only supports global/Levenshtein search")
  }
  invisible(mode)
}

seqtrie_check_startree_fixed_search <- function(x,
                                                max_distance,
                                                max_fraction,
                                                mode,
                                                cost_matrix,
                                                gap_cost,
                                                gap_open_cost,
                                                lower_triangle,
                                                match_mode,
                                                nthreads,
                                                show_progress) {
  if (!is.null(max_distance)) {
    max_distance <- seqtrie_check_startree_integer(max_distance, "max_distance")
    if (!identical(max_distance, x@max_distance)) {
      stop("StarTree max_distance is fixed at construction")
    }
  }
  if (!is.null(max_fraction)) {
    stop("StarTree does not support max_fraction")
  }
  seqtrie_check_startree_mode(mode)
  if (!is.null(cost_matrix)) {
    stop("StarTree does not support custom cost_matrix values; use mismatch_cost and gap_cost at construction")
  }
  if (!is_missing_arg(gap_cost)) {
    gap_cost <- seqtrie_check_startree_integer(gap_cost, "gap_cost", positive = TRUE)
    if (!identical(gap_cost, x@gap_cost)) {
      stop("StarTree gap_cost is fixed at construction")
    }
  }
  if (!is_missing_arg(gap_open_cost)) {
    stop("StarTree does not support affine gap penalties")
  }
  lower_triangle <- check_flag(lower_triangle, "lower_triangle")
  if (lower_triangle) {
    stop("StarTree align_search does not support lower_triangle")
  }
  match_mode <- match.arg(match_mode, c("all", "best"))
  if (match_mode != "all") {
    stop("StarTree align_search does not support match_mode = 'best'")
  }
  # nthreads and show_progress are per-call for the query search (the query
  # search is independent of how the construction-time self-join was threaded),
  # so they are validated and used directly rather than fixed at construction.
  nthreads <- check_threads(nthreads)
  show_progress <- check_flag(show_progress, "show_progress")
  list(nthreads = nthreads, show_progress = show_progress)
}

# S7 generic --------------------------------------------------------------------

#' Return a fixed-tree self-similarity result
#'
#' @param x A seqtrie S7 object.
#' @return A data frame with columns `query`, `target`, and `distance`.
#' @export
result <- S7::new_generic("result", "x", function(x) S7::S7_dispatch())

# S7 class ----------------------------------------------------------------------

#' Starcode-style fixed tree
#'
#' `star_tree()` constructs a fixed DNA-only tree using a modified version of
#' the Starcode all-pairs search strategy, adapted to operate over a radix trie. The input
#' sequences, `max_distance`, `mismatch_cost`, and `gap_cost` are fixed at
#' construction, and the self-similarity join runs immediately. Use `result()`
#' to retrieve that self-join, and
#' [align_search()] to search additional query sequences against the fixed
#' target set.
#'
#' StarTree supports global/Levenshtein-style DNA alignment only. It accepts
#' `A`, `C`, `G`, `T`, and `N` in either case; sequences are stored and returned
#' in uppercase. `N` is treated as a regular ambiguous base with mismatch cost,
#' not as a wildcard. Custom substitution matrices, affine gaps, anchored
#' search, Hamming search, insertion, and deletion are not supported.
#'
#' The algorithm is based on Starcode (Zorita, Cuscó, and Filion 2015)
#' \doi{10.1093/bioinformatics/btv053}.
#'
#' @param sequences A required character vector of DNA sequences.
#' @param max_distance A single non-negative integer distance threshold.
#' @param mismatch_cost A single positive integer mismatch cost.
#' @param gap_cost A single positive integer gap cost.
#' @param nthreads `r rdoc("nthreads")`
#' @param show_progress `r rdoc("show_progress")`
#' @return A `star_tree` object.
#' @seealso [align_search()]
#' @examples
#' tree <- star_tree(c("ACGT", "ACGA", "AAAA"), max_distance = 1)
#' result(tree)
#' align_search(tree, c("ACGT", "AAAT"))
#' @export
star_tree <- S7::new_class(
  "star_tree",
  properties = list(
    tree_pointer = S7::class_any,
    max_distance = S7::class_integer,
    mismatch_cost = S7::class_integer,
    gap_cost = S7::class_integer,
    nthreads = S7::class_integer,
    show_progress = S7::class_logical
  ),
  constructor = function(sequences,
                         max_distance,
                         mismatch_cost = 1L,
                         gap_cost = 1L,
                         nthreads = 1L,
                         show_progress = FALSE) {
    seqtrie_check_startree_sequences(sequences, "sequences")
    params <- seqtrie_check_startree_params(max_distance, mismatch_cost, gap_cost)
    nthreads <- check_threads(nthreads)
    show_progress <- check_flag(show_progress, "show_progress")
    S7::new_object(
      star_tree,
      tree_pointer = StarTree_create(
        sequences,
        params$max_distance,
        params$mismatch_cost,
        params$gap_cost,
        nthreads,
        show_progress
      ),
      max_distance = params$max_distance,
      mismatch_cost = params$mismatch_cost,
      gap_cost = params$gap_cost,
      nthreads = nthreads,
      show_progress = show_progress
    )
  }
)

S7::method(result, star_tree) <- function(x) {
  StarTree_result(x@tree_pointer)
}

S7::method(to_vector, star_tree) <- function(x) {
  StarTree_to_vector(x@tree_pointer)
}

S7::method(size, star_tree) <- function(x) {
  StarTree_size(x@tree_pointer)
}

S7::method(align_search, star_tree) <- function(x,
                                                query,
                                                max_distance = NULL,
                                                max_fraction = NULL,
                                                mode = "levenshtein",
                                                cost_matrix = NULL,
                                                gap_cost = NA_integer_,
                                                gap_open_cost = NA_integer_,
                                                lower_triangle = FALSE,
                                                match_mode = c("all", "best"),
                                                nthreads = 1L,
                                                show_progress = FALSE) {
  opts <- seqtrie_check_startree_fixed_search(
    x,
    max_distance = max_distance,
    max_fraction = max_fraction,
    mode = mode,
    cost_matrix = cost_matrix,
    gap_cost = gap_cost,
    gap_open_cost = gap_open_cost,
    lower_triangle = lower_triangle,
    match_mode = match_mode,
    nthreads = nthreads,
    show_progress = show_progress
  )

  seqtrie_check_startree_sequences(query, "query", allow_na = TRUE)
  keep <- not_na_character(query)
  if (!any(keep)) {
    return(seqtrie_empty_match_result())
  }

  StarTree_search(
    x@tree_pointer,
    query[keep],
    opts$nthreads,
    opts$show_progress
  )
}
