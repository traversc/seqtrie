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
  if (!mode %in% c("global", "anchored", "hamming")) {
    stop("StarTree mode must be one of global (gb, lv, levenshtein), anchored (an, en, extension), or hamming (hm)")
  }
  mode
}

seqtrie_dna_scalar_cost_matrix <- function(query, target, mismatch_cost) {
  chars <- sort(unique(unlist(strsplit(c(query, target), "", fixed = TRUE))))
  cost_matrix <- generate_cost_matrix(paste0(chars, collapse = ""), match = 0L, mismatch = mismatch_cost)
  if ("N" %in% chars) {
    cost_matrix["N", "N"] <- mismatch_cost
  }
  cost_matrix
}

seqtrie_add_anchored_startree_sizes <- function(result, mismatch_cost, gap_cost, nthreads) {
  if (nrow(result) == 0L) {
    return(result)
  }
  cost_matrix <- seqtrie_dna_scalar_cost_matrix(result$query, result$target, mismatch_cost)
  result2 <- c_dist_pairwise(
    result$query,
    result$target,
    mode = "anchored",
    cost_matrix = cost_matrix,
    gap_cost = gap_cost,
    gap_open_cost = NA_integer_,
    nthreads = nthreads,
    show_progress = FALSE
  )
  if (any(result$distance != result2)) {
    stop("Internal error: anchored StarTree search results do not match pairwise results")
  }
  result$query_size <- attr(result2, "query_size")
  result$target_size <- attr(result2, "target_size")
  result
}

seqtrie_stop_startree_query_param <- function(class_name, param) {
  stop(class_name, " ", param, " is fixed at construction or unsupported; omit ",
       param, " in align_search()")
}

seqtrie_check_startree_fixed_search <- function(class_name,
                                                max_distance_supplied,
                                                max_fraction_supplied,
                                                mode_supplied,
                                                cost_matrix_supplied,
                                                gap_cost_supplied,
                                                gap_open_cost_supplied,
                                                lower_triangle_supplied,
                                                match_mode_supplied,
                                                nthreads,
                                                show_progress) {
  if (max_distance_supplied) {
    seqtrie_stop_startree_query_param(class_name, "max_distance")
  }
  if (max_fraction_supplied) {
    seqtrie_stop_startree_query_param(class_name, "max_fraction")
  }
  if (mode_supplied) {
    seqtrie_stop_startree_query_param(class_name, "mode")
  }
  if (cost_matrix_supplied) {
    seqtrie_stop_startree_query_param(class_name, "cost_matrix")
  }
  if (gap_cost_supplied) {
    seqtrie_stop_startree_query_param(class_name, "gap_cost")
  }
  if (gap_open_cost_supplied) {
    seqtrie_stop_startree_query_param(class_name, "gap_open_cost")
  }
  if (lower_triangle_supplied) {
    seqtrie_stop_startree_query_param(class_name, "lower_triangle")
  }
  if (match_mode_supplied) {
    seqtrie_stop_startree_query_param(class_name, "match_mode")
  }
  nthreads <- check_threads(nthreads)
  show_progress <- check_flag(show_progress, "show_progress")
  list(nthreads = nthreads, show_progress = show_progress)
}

# S7 generic --------------------------------------------------------------------

#' Return a fixed-tree self-similarity result
#'
#' @param x A seqtrie S7 object.
#' @return A data frame with columns `query`, `target`, and `distance`.
#' Anchored mode also includes `query_size` and `target_size`.
#' @export
result <- S7::new_generic("result", "x", function(x) S7::S7_dispatch())

seqtrie_is_s7_dispatch_call <- function(call) {
  if (!is.call(call)) {
    return(FALSE)
  }
  head <- call[[1L]]
  if (identical(head, quote(S7_dispatch))) {
    return(TRUE)
  }
  is.call(head) &&
    identical(head[[1L]], quote(`::`)) &&
    identical(head[[2L]], quote(S7)) &&
    identical(head[[3L]], quote(S7_dispatch))
}

seqtrie_align_search_supplied_args <- function() {
  calls <- sys.calls()
  dispatch <- which(vapply(calls, seqtrie_is_s7_dispatch_call, logical(1)))
  if (length(dispatch) == 0L) {
    return(character())
  }
  generic_call_index <- dispatch[[length(dispatch)]] - 1L
  if (generic_call_index < 1L) {
    return(character())
  }
  call <- tryCatch(
    match.call(
      definition = align_search,
      call = calls[[generic_call_index]],
      expand.dots = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(call)) {
    return(character())
  }
  setdiff(names(as.list(call))[-1L], c("", "..."))
}

# S7 class ----------------------------------------------------------------------

#' Starcode-style fixed tree
#'
#' `star_tree()` constructs a fixed DNA-only tree using a modified version of
#' the Starcode all-pairs search strategy, adapted to operate over a radix trie.
#' The input sequences, alignment `mode`, `max_distance`, `mismatch_cost`, and
#' `gap_cost` are fixed at construction, and the self-similarity join runs
#' immediately. Use `result()` to retrieve that self-join, and [align_search()]
#' to search additional query sequences against the fixed target set.
#'
#' StarTree supports global/Levenshtein-style, anchored, and Hamming DNA
#' alignment. It accepts `A`, `C`, `G`, `T`, and `N` in either case; sequences
#' are stored and returned in uppercase. `N` is treated as a regular ambiguous
#' base with mismatch cost, not as a wildcard. Custom substitution matrices,
#' affine gaps, insertion, and deletion are not supported.
#'
#' Hamming mode (`mode = "hamming"`) is substitution-only: only equal-length
#' sequences can match, and `max_distance` is the maximum number of mismatching
#' positions (unit substitution cost; `mismatch_cost` and `gap_cost` do not
#' apply). It is typically much faster than global mode for the same data.
#'
#' For `star_tree` objects, [align_search()] only accepts `query`, `nthreads`,
#' and `show_progress`; all alignment parameters are fixed here at construction.
#' Anchored-mode results also include `query_size` and `target_size`.
#'
#' The algorithm is based on Starcode (Zorita, Cuscó, and Filion 2015)
#' \doi{10.1093/bioinformatics/btv053}.
#'
#' @param sequences A required character vector of DNA sequences.
#' @param max_distance A single non-negative integer distance threshold.
#' @param mode Alignment mode: global/Levenshtein, anchored, or hamming.
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
#'
#' anchored <- star_tree(c("ACGT", "ACG", "AAAA"), max_distance = 1,
#'                       mode = "anchored")
#' result(anchored)
#'
#' hamming <- star_tree(c("ACGT", "ACGA", "TCGT"), max_distance = 1,
#'                      mode = "hamming")
#' result(hamming)
#' @export
star_tree <- S7::new_class(
  "star_tree",
  properties = list(
    tree_pointer = S7::class_any,
    mode = S7::class_character,
    max_distance = S7::class_integer,
    mismatch_cost = S7::class_integer,
    gap_cost = S7::class_integer,
    nthreads = S7::class_integer,
    show_progress = S7::class_logical
  ),
  constructor = function(sequences,
                         max_distance,
                         mode = "levenshtein",
                         mismatch_cost = 1L,
                         gap_cost = 1L,
                         nthreads = 1L,
                         show_progress = FALSE) {
    seqtrie_check_startree_sequences(sequences, "sequences")
    mode <- seqtrie_check_startree_mode(mode)
    params <- seqtrie_check_startree_params(max_distance, mismatch_cost, gap_cost)
    nthreads <- check_threads(nthreads)
    show_progress <- check_flag(show_progress, "show_progress")
    S7::new_object(
      star_tree,
      tree_pointer = if (mode == "anchored") AnchoredStarTree_create(
        sequences,
        params$max_distance,
        params$mismatch_cost,
        params$gap_cost,
        nthreads,
        show_progress
      ) else StarTree_create(
        sequences,
        params$max_distance,
        params$mismatch_cost,
        params$gap_cost,
        nthreads,
        show_progress,
        hamming = identical(mode, "hamming")
      ),
      mode = mode,
      max_distance = params$max_distance,
      mismatch_cost = params$mismatch_cost,
      gap_cost = params$gap_cost,
      nthreads = nthreads,
      show_progress = show_progress
    )
  }
)

S7::method(result, star_tree) <- function(x) {
  if (identical(S7::prop(x, "mode"), "anchored")) {
    result <- AnchoredStarTree_result(S7::prop(x, "tree_pointer"))
    return(seqtrie_add_anchored_startree_sizes(
      result,
      S7::prop(x, "mismatch_cost"),
      S7::prop(x, "gap_cost"),
      S7::prop(x, "nthreads")
    ))
  }
  StarTree_result(S7::prop(x, "tree_pointer"))
}

S7::method(to_vector, star_tree) <- function(x) {
  if (identical(S7::prop(x, "mode"), "anchored")) {
    return(AnchoredStarTree_to_vector(S7::prop(x, "tree_pointer")))
  }
  StarTree_to_vector(S7::prop(x, "tree_pointer"))
}

S7::method(size, star_tree) <- function(x) {
  if (identical(S7::prop(x, "mode"), "anchored")) {
    return(AnchoredStarTree_size(S7::prop(x, "tree_pointer")))
  }
  StarTree_size(S7::prop(x, "tree_pointer"))
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
  class_name <- "StarTree"
  supplied <- seqtrie_align_search_supplied_args()
  opts <- seqtrie_check_startree_fixed_search(
    class_name = class_name,
    max_distance_supplied = "max_distance" %in% supplied,
    max_fraction_supplied = "max_fraction" %in% supplied,
    mode_supplied = "mode" %in% supplied,
    cost_matrix_supplied = "cost_matrix" %in% supplied,
    gap_cost_supplied = "gap_cost" %in% supplied,
    gap_open_cost_supplied = "gap_open_cost" %in% supplied,
    lower_triangle_supplied = "lower_triangle" %in% supplied,
    match_mode_supplied = "match_mode" %in% supplied,
    nthreads = nthreads,
    show_progress = show_progress
  )

  seqtrie_check_startree_sequences(query, "query", allow_na = TRUE)
  keep <- not_na_character(query)
  if (!any(keep)) {
    return(seqtrie_empty_match_result())
  }

  if (identical(S7::prop(x, "mode"), "anchored")) {
    result <- AnchoredStarTree_search(
      S7::prop(x, "tree_pointer"),
      query[keep],
      opts$nthreads,
      opts$show_progress
    )
    return(seqtrie_add_anchored_startree_sizes(
      result,
      S7::prop(x, "mismatch_cost"),
      S7::prop(x, "gap_cost"),
      opts$nthreads
    ))
  }

  StarTree_search(
    S7::prop(x, "tree_pointer"),
    query[keep],
    opts$nthreads,
    opts$show_progress
  )
}
