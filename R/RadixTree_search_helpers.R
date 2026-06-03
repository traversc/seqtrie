#' @title Distance search for similar sequences
#' @description Find similar sequences within a distance threshold
#' @param query `r rdoc("query")`
#' @param target `r rdoc("target")` If `NULL`, `query` is searched against itself and only pairs where the query index is strictly greater than the target terminal index are returned.
#' @param max_distance `r rdoc("max_distance")`
#' @param max_fraction `r rdoc("max_fraction")`
#' @param mode `r rdoc("mode")`
#' @param cost_matrix `r rdoc("cost_matrix")`
#' @param gap_cost `r rdoc("gap_cost")`
#' @param gap_open_cost `r rdoc("gap_open_cost")`
#' @param tree_class Which tree implementation to use. One of RadixTree, RadixForest, StarTree, radix_tree, radix_forest, or star_tree (default: RadixTree)
#' @param nthreads `r rdoc("nthreads")`
#' @param show_progress `r rdoc("show_progress")`
#' @param mismatch_cost A single positive integer mismatch cost for fixed StarTree classes.
#' @details
#' This function finds all sequences in _target_ that are within a distance threshold of any sequence in _query_.
#' If `target = NULL`, the tree is built from `query` and each query is searched against that tree while requiring the query index to be strictly greater than the target terminal index. This returns lower-triangle self-pairs. Duplicate query strings are not given special handling; because the underlying tree stores one terminal index per unique sequence, duplicates naturally collapse to the first inserted occurrence.
#' This function uses a radix_tree/RadixTree, radix_forest/RadixForest, or fixed star_tree/StarTree to store _target_ sequences. Use `tree_class = "StarTree"` with `mode = "anchored"` for fixed anchored DNA joins.
#' 
#' `r rdoc("details")`
#' @return The output is a data frame of all matches with columns "query", "target", and "distance".
#' For anchored searches, the output also includes columns "query_size" and "target_size"
#' containing the portion of the query and target sequences that are aligned.
#' @examples
#' dist_search(c("ACGT", "AAAA"), c("ACG", "ACGT"), max_distance = 1, mode = "levenshtein")
#' @name dist_search
dist_search <- function(query, target = NULL, max_distance = NULL, max_fraction = NULL, mode = "levenshtein",
                        cost_matrix = NULL, gap_cost = NA_integer_, gap_open_cost = NA_integer_, tree_class = "RadixTree",
                        nthreads = 1, show_progress = FALSE, mismatch_cost = 1L) {
  if (!tree_class %in% c("RadixTree", "RadixForest", "StarTree",
                         "radix_tree", "radix_forest", "star_tree")) {
    stop("tree_class must be one of RadixTree, RadixForest, StarTree, radix_tree, radix_forest, or star_tree")
  }
  cpp_tree_class <- switch(
    tree_class,
    radix_tree = "RadixTree",
    radix_forest = "RadixForest",
    star_tree = "StarTree",
    tree_class
  )

  if (cpp_tree_class == "StarTree") {
    nthreads <- check_threads(nthreads)
    show_progress <- check_flag(show_progress, "show_progress")
    if (!is.null(max_fraction)) {
      stop(cpp_tree_class, " does not support max_fraction")
    }
    star_gap_cost <- if (is_missing_arg(gap_cost)) 1L else gap_cost
    star_mode <- seqtrie_check_startree_mode(mode)
    params <- seqtrie_check_startree_params(max_distance, mismatch_cost, star_gap_cost)
    if (!is.null(cost_matrix)) {
      stop(cpp_tree_class, " does not support custom cost_matrix values; use mismatch_cost and gap_cost")
    }
    if (!is_missing_arg(gap_open_cost)) {
      stop(cpp_tree_class, " does not support affine gap penalties")
    }
    if (is.null(target)) {
      seqtrie_check_startree_sequences(query, "query")
      if (star_mode == "anchored") {
        result <- AnchoredStarTree_self_search(
          query,
          max_distance = params$max_distance,
          mismatch_cost = params$mismatch_cost,
          gap_cost = params$gap_cost,
          nthreads = nthreads,
          show_progress = show_progress
        )
        return(seqtrie_add_anchored_startree_sizes(
          result, params$mismatch_cost, params$gap_cost, nthreads
        ))
      }
      return(StarTree_self_search(
        query,
        max_distance = params$max_distance,
        mismatch_cost = params$mismatch_cost,
        gap_cost = params$gap_cost,
        nthreads = nthreads,
        show_progress = show_progress,
        hamming = identical(star_mode, "hamming")
      ))
    }
    obj <- star_tree(
      target,
      max_distance = params$max_distance,
      mode = star_mode,
      mismatch_cost = params$mismatch_cost,
      gap_cost = params$gap_cost,
      nthreads = nthreads,
      show_progress = show_progress
    )
    return(align_search(obj, query = query, nthreads = nthreads, show_progress = show_progress))
  }

  if (is.null(target)) {
    obj <- if (cpp_tree_class == "RadixTree") radix_tree(query) else radix_forest(query)
    return(align_search(
      obj,
      query = query,
      max_distance = max_distance,
      max_fraction = max_fraction,
      mode = mode,
      cost_matrix = cost_matrix,
      gap_cost = gap_cost,
      gap_open_cost = gap_open_cost,
      lower_triangle = TRUE,
      nthreads = nthreads,
      show_progress = show_progress
    ))
  }

  if (cpp_tree_class == "RadixTree") {
    obj <- radix_tree(target)
    align_search(
      obj,
      query = query,
      max_distance = max_distance,
      max_fraction = max_fraction,
      mode = mode,
      cost_matrix = cost_matrix,
      gap_cost = gap_cost,
      gap_open_cost = gap_open_cost,
      nthreads = nthreads,
      show_progress = show_progress
    )
  } else if (cpp_tree_class == "RadixForest") {
    obj <- radix_forest(target)
    align_search(
      obj,
      query = query,
      max_distance = max_distance,
      max_fraction = max_fraction,
      mode = mode,
      cost_matrix = cost_matrix,
      gap_cost = gap_cost,
      gap_open_cost = gap_open_cost,
      nthreads = nthreads,
      show_progress = show_progress
    )
  }
}

#' @title split_search
#' @description Search for similar sequences based on splitting sequences into left and right sides
#' and searching for matches on each side using a bidirectional anchored alignment.
#' @param query `r rdoc("query")`
#' @param target `r rdoc("target")`
#' @param query_split index to split query sequence. Should be within (edge_trim, nchar(query)-edge_trim] or -1 to indicate no split.
#' @param target_split index to split target sequence. Should be within (edge_trim, nchar(target)-edge_trim] or -1 to indicate no split.
#' @param edge_trim number of bases to trim from each side of the sequence (default value: 0).
#' @param max_distance `r rdoc("max_distance")`
#' @param max_fraction `r rdoc("max_fraction")`
#' @param cost_matrix `r rdoc("cost_matrix")`
#' @param gap_cost `r rdoc("gap_cost")`
#' @param gap_open_cost `r rdoc("gap_open_cost")`
#' @param nthreads `r rdoc("nthreads")`
#' @param show_progress `r rdoc("show_progress")`
#' @return A data frame with columns query, target, and distance.
#' @details
#' This function is useful for searching for similar sequences that may have variable sequencing windows (e.g. different 5' and 3' primers)
#' but contain the same core sequence or position.
#' The two split parameters partition the query and target sequences into left and right sides, 
#' where left = rev(substr(sequence, edge_trim+1, split)) and right = substr(sequence, split+1, nchar(sequence)-edge_trim).
#' @examples
#' # Consider two sets of sequences
#' # query1   AGACCTAA CCC
#' # target1 AAGACCTAA CC
#' # query2   GGGTGTAA CCACCC
#' # target2   GGTGTAA CCAC
#' # Despite having different frames, query1 and query2 can clearly
#' # match to target1 and target2, respectively.
#' # One could consider splitting based on a common core sequence, 
#' # e.g. a common TAA stop codon. 
#' split_search(query=c(  "AGACCTAACCC", "GGGTGTAACCACCC"),
#'              target=c("AAGACCTAACC",   "GGTGTAACCAC"),
#'              query_split=c(8, 8),
#'              target_split=c(9, 7),
#'              edge_trim=0,
#'              max_distance=0)
split_search <- function(query, target, query_split, target_split, edge_trim = 0L,
                         max_distance = 0L, max_fraction = NULL,
                         cost_matrix = NULL, gap_cost = NA_integer_, gap_open_cost = NA_integer_,
                         nthreads = 1, show_progress = FALSE) {
  if (!is.character(query) || !is.character(target)) {
    stop("query and target must be character vectors")
  }
  nthreads <- check_threads(nthreads)
  if (!is_integerlike(edge_trim) || length(edge_trim) != 1L || is.na(edge_trim) || edge_trim < 0L) {
    stop("edge_trim must be a single non-negative integer")
  }
  edge_trim <- as.integer(edge_trim)

  query_split <- recycle_arg(query_split, query)
  target_split <- recycle_arg(target_split, target)

  if (!is.null(max_distance)) {
    max_distance <- recycle_arg(max_distance, query)
  } else if (!is.null(max_fraction)) {
    max_distance <- as.integer(nchar(query) * recycle_arg(max_fraction, query))
  } else {
    stop("Either max_distance or max_fraction must be non-null")
  }

  keep_query <- not_na_character(query)
  keep_target <- not_na_character(target)
  query <- query[keep_query]
  query_split <- query_split[keep_query]
  max_distance <- max_distance[keep_query]
  target <- target[keep_target]
  target_split <- target_split[keep_target]

  if (length(query) == 0L || length(target) == 0L) {
    return(data.frame(query = character(), target = character(), distance = integer()))
  }

  normalize_split <- function(split, sequences, label) {
    if (!is_integerlike(split) || any(is.na(split))) {
      stop(label, " must contain integer-like split positions")
    }
    split <- as.integer(split)
    if (any(split < 0L & split != -1L)) {
      stop(label, " values must be -1 or non-negative split positions")
    }
    seq_len <- nchar(sequences)
    split <- ifelse(split == -1L, seq_len - edge_trim, split)
    invalid <- split <= edge_trim | split > (seq_len - edge_trim)
    if (any(invalid)) {
      stop(label, " values must be -1 or within (edge_trim, nchar(sequence) - edge_trim]")
    }
    split
  }

  query_split <- normalize_split(query_split, query, "query_split")
  target_split <- normalize_split(target_split, target, "target_split")

  charset <- unique(c(get_charset(query), get_charset(target)))
  check_alignment_params("anchored", cost_matrix, gap_cost, gap_open_cost, charset, diag_must_be_zero = FALSE)

  if (!is_integerlike(max_distance) || any(is.na(max_distance)) || any(max_distance < 0)) {
    stop("max_distance/max_fraction must be non-negative integer-like values")
  }
  max_distance <- as.integer(max_distance)

  if (!is.na(gap_open_cost) && !is.na(gap_cost) && gap_open_cost > 0L) {
    gap_open_cost <- gap_open_cost + gap_cost
  }

  c_split_search(query, target, query_split, target_split, edge_trim, max_distance,
                 cost_matrix, gap_cost, gap_open_cost, nthreads, show_progress)
}
