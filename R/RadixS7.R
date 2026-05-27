# Shared implementation helpers -------------------------------------------------

seqtrie_empty_match_result <- function() {
  data.frame(query = character(0), target = character(0), distance = integer(0), stringsAsFactors = FALSE)
}

seqtrie_empty_prefix_result <- function() {
  data.frame(query = character(0), target = character(0), stringsAsFactors = FALSE)
}

seqtrie_plot_graph <- function(result, root_label, plot) {
  if (is.null(result)) {
    return(data.frame(parent = character(0), child = character(0), stringsAsFactors = FALSE))
  }
  if (!plot) {
    return(result)
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph package is required to plot the tree.") # nocov
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required to plot the tree.") # nocov
  }

  result$parent <- ifelse(result$parent == "", root_label, result$parent)
  gr <- igraph::graph_from_data_frame(result)
  fr <- igraph::layout.fruchterman.reingold(gr)
  fr <- as.data.frame(fr)
  fr$node <- names(igraph::V(gr))
  fr$fill <- ifelse(fr$node == root_label, "white", "skyblue")
  fr$size <- ifelse(fr$node == root_label, 16, 12)

  result$parent_x <- fr$V1[match(result$parent, fr$node)]
  result$parent_y <- fr$V2[match(result$parent, fr$node)]
  result$child_x <- fr$V1[match(result$child, fr$node)]
  result$child_y <- fr$V2[match(result$child, fr$node)]

  ggplot2::ggplot() +
    ggplot2::geom_segment(data = result, ggplot2::aes(x = parent_x, xend = child_x, y = parent_y, yend = child_y)) +
    ggplot2::geom_point(data = fr, ggplot2::aes(x = V1, y = V2, fill = fill, size = size), shape = 21, color = "black") +
    ggplot2::geom_text(data = fr, ggplot2::aes(x = V1, y = V2, label = node)) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_size_identity() +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title = ggplot2::element_blank())
}

seqtrie_prepare_search <- function(query,
                                   max_distance,
                                   max_fraction,
                                   mode,
                                   cost_matrix,
                                   gap_cost,
                                   gap_open_cost,
                                   lower_triangle,
                                   match_mode,
                                   charset,
                                   allow_anchored = TRUE,
                                   nthreads) {
  check_character_vector(query, "query")
  nthreads <- check_threads(nthreads)
  lower_triangle <- check_flag(lower_triangle, "lower_triangle")
  match_mode <- match.arg(match_mode, c("all", "best"))
  keep <- not_na_character(query)
  query_search <- query[keep]

  charset <- unique(c(charset, get_charset(query_search)))
  check_alignment_params(mode, cost_matrix, gap_cost, gap_open_cost, charset, diag_must_be_zero = FALSE)
  mode <- normalize_mode_parameter(mode)
  if (!allow_anchored && !mode %in% c("hamming", "global")) {
    stop("mode must be one of hamming (hm) or global (gb, lv, levenshtein)")
  }

  if (!is.null(max_distance)) {
    max_distance <- recycle_arg(max_distance, query)
  } else if (!is.null(max_fraction)) {
    max_distance <- as.integer(nchar(query) * max_fraction)
  } else {
    stop("Either max_distance or max_fraction must be non-null")
  }
  if (any(is.na(max_distance[keep])) || any(max_distance[keep] < 0)) {
    stop("max_distance/max_fraction must be non-negative")
  }

  if (!is.na(gap_open_cost) && !is.na(gap_cost) && gap_open_cost > 0L) {
    gap_open_cost <- gap_open_cost + gap_cost
  }

  list(
    query = query,
    keep = keep,
    query_search = query_search,
    query_index = which(keep),
    max_distance = max_distance,
    mode = mode,
    gap_cost = gap_cost,
    gap_open_cost = gap_open_cost,
    lower_triangle = lower_triangle,
    match_mode = match_mode,
    nthreads = nthreads
  )
}

seqtrie_add_anchored_sizes <- function(result, mode, cost_matrix, gap_cost, gap_open_cost, nthreads) {
  if (mode != "anchored") {
    return(result)
  }
  result2 <- c_dist_pairwise(result$query, result$target, mode, cost_matrix, gap_cost, gap_open_cost, nthreads, show_progress = FALSE)
  if (any(result$distance != result2)) {
    stop("Internal error: anchored search results do not match pairwise results")
  }
  result$query_size <- attr(result2, "query_size")
  result$target_size <- attr(result2, "target_size")
  result
}

# S7 generics -------------------------------------------------------------------

#' Insert sequences
#'
#' @param x A seqtrie S7 object.
#' @param sequences `r rdoc("sequences")`
#' @return A logical vector indicating whether each sequence was inserted.
#' @export
insert <- S7::new_generic("insert", "x", function(x, sequences) S7::S7_dispatch())

#' Erase sequences
#'
#' @param x A seqtrie S7 object.
#' @param sequences `r rdoc("sequences")`
#' @return A logical vector indicating whether each sequence was erased.
#' @export
erase <- S7::new_generic("erase", "x", function(x, sequences) S7::S7_dispatch())

#' Test sequence membership
#'
#' @param x A seqtrie S7 object.
#' @param query `r rdoc("query")`
#' @return A logical vector indicating whether each query is present.
#' @export
has_sequence <- S7::new_generic("has_sequence", "x", function(x, query) S7::S7_dispatch())

#' Prefix search
#'
#' @param x A seqtrie S7 object.
#' @param query `r rdoc("query")`
#' @return A data frame with columns `query` and `target`.
#' @export
prefix_search <- S7::new_generic("prefix_search", "x", function(x, query) S7::S7_dispatch())

#' Alignment search
#'
#' @param x A seqtrie S7 object.
#' @param query `r rdoc("query")`
#' @param max_distance `r rdoc("max_distance")`
#' @param max_fraction `r rdoc("max_fraction")`
#' @param mode `r rdoc("mode")`
#' @param cost_matrix `r rdoc("cost_matrix")`
#' @param gap_cost `r rdoc("gap_cost")`
#' @param gap_open_cost `r rdoc("gap_open_cost")`
#' @param lower_triangle `r rdoc("lower_triangle")`
#' @param match_mode `r rdoc("match_mode")`
#' @param nthreads `r rdoc("nthreads")`
#' @param show_progress `r rdoc("show_progress")`
#' @param ... Additional method-specific arguments.
#' @return A data frame with columns `query`, `target`, and `distance`.
#' @export
align_search <- S7::new_generic(
  "align_search", "x",
  function(x, query, max_distance = NULL, max_fraction = NULL,
           mode = "levenshtein", cost_matrix = NULL, gap_cost = NA_integer_,
           gap_open_cost = NA_integer_, lower_triangle = FALSE,
           match_mode = c("all", "best"), nthreads = 1L,
           show_progress = FALSE, ...) S7::S7_dispatch()
)

#' Single-gap alignment search
#'
#' @param x A `radix_tree` object.
#' @param query `r rdoc("query")`
#' @param max_distance `r rdoc("max_distance")`
#' @param gap_cost `r rdoc("gap_cost")`
#' @param nthreads `r rdoc("nthreads")`
#' @param show_progress `r rdoc("show_progress")`
#' @return A data frame with columns `query`, `target`, and `distance`.
#' @export
single_gap_search <- S7::new_generic(
  "single_gap_search", "x",
  function(x, query, max_distance, gap_cost = 1L, nthreads = 1L,
           show_progress = FALSE) S7::S7_dispatch()
)

#' Convert sequences to a character vector
#'
#' @param x A seqtrie S7 object.
#' @return A character vector of stored sequences.
#' @export
to_vector <- S7::new_generic("to_vector", "x", function(x) S7::S7_dispatch())

#' Convert a trie or forest to a string
#'
#' @param x A seqtrie S7 object.
#' @return A string representation of the trie or forest.
#' @export
to_string <- S7::new_generic("to_string", "x", function(x) S7::S7_dispatch())

#' Count stored sequences
#'
#' @param x A seqtrie S7 object.
#' @return The number of stored sequences.
#' @export
size <- S7::new_generic("size", "x", function(x) S7::S7_dispatch())

#' Plot trie or forest structure
#'
#' @param x A seqtrie S7 object.
#' @param depth The tree depth to plot. If -1, plot the entire tree.
#' @param root_label The label of the root node in the plot.
#' @param plot Whether to create a plot or return the graph data.
#' @return A data frame of parent-child relationships or a ggplot2 object.
#' @export
plot_tree <- S7::new_generic(
  "plot_tree", "x",
  function(x, depth = -1, root_label = "root", plot = TRUE) S7::S7_dispatch()
)

#' Validate trie or forest structure
#'
#' @param x A seqtrie S7 object.
#' @return A logical value indicating whether the object is valid.
#' @export
is_valid <- S7::new_generic("is_valid", "x", function(x) S7::S7_dispatch())

# S7 classes --------------------------------------------------------------------

#' Radix tree
#'
#' `radix_tree()` constructs a mutable S7 wrapper around the seqtrie C++ radix
#' tree implementation. It supports hamming, global/Levenshtein, anchored, and
#' single-gap searches.
#'
#' @param sequences Optional character vector of sequences to insert.
#' @return A `radix_tree` object.
#' @seealso [RadixTree] for the R6 compatibility wrapper.
#' @examples
#' tree <- radix_tree(c("ACGT", "AAAA"))
#' align_search(tree, "ACG", max_distance = 1, mode = "levenshtein")
#' @export
radix_tree <- S7::new_class(
  "radix_tree",
  properties = list(
    root_pointer = S7::class_any,
    char_counter_pointer = S7::class_any
  ),
  constructor = function(sequences = NULL) {
    x <- S7::new_object(
      radix_tree,
      root_pointer = RadixTree_create(),
      char_counter_pointer = CharCounter_create()
    )
    if (!is.null(sequences)) {
      insert(x, sequences)
    }
    x
  }
)

#' Radix forest
#'
#' `radix_forest()` constructs a mutable S7 wrapper around the seqtrie C++ radix
#' forest implementation. It partitions sequences by length and supports hamming
#' and global/Levenshtein searches.
#'
#' @param sequences Optional character vector of sequences to insert.
#' @return A `radix_forest` object.
#' @seealso [RadixForest] for the R6 compatibility wrapper.
#' @examples
#' forest <- radix_forest(c("ACGT", "AAAA"))
#' align_search(forest, "ACG", max_distance = 1, mode = "levenshtein")
#' @export
radix_forest <- S7::new_class(
  "radix_forest",
  properties = list(
    forest_pointer = S7::class_any,
    char_counter_pointer = S7::class_any
  ),
  constructor = function(sequences = NULL) {
    x <- S7::new_object(
      radix_forest,
      forest_pointer = RadixForest_create(),
      char_counter_pointer = CharCounter_create()
    )
    if (!is.null(sequences)) {
      insert(x, sequences)
    }
    x
  }
)

# Radix tree methods ------------------------------------------------------------

S7::method(insert, radix_tree) <- function(x, sequences) {
  check_character_vector(sequences, "sequences")
  keep <- not_na_character(sequences)
  result <- rep(FALSE, length(sequences))
  if (any(keep)) {
    result[keep] <- RadixTree_insert(S7::prop(x, "root_pointer"), sequences[keep])
    CharCounter_add(S7::prop(x, "char_counter_pointer"), sequences[keep][result[keep]])
  }
  invisible(result)
}

S7::method(erase, radix_tree) <- function(x, sequences) {
  check_character_vector(sequences, "sequences")
  keep <- not_na_character(sequences)
  result <- rep(FALSE, length(sequences))
  if (any(keep)) {
    result[keep] <- RadixTree_erase(S7::prop(x, "root_pointer"), sequences[keep])
    CharCounter_subtract(S7::prop(x, "char_counter_pointer"), sequences[keep][result[keep]])
  }
  invisible(result)
}

S7::method(has_sequence, radix_tree) <- function(x, query) {
  check_character_vector(query, "query")
  keep <- not_na_character(query)
  result <- rep(FALSE, length(query))
  if (any(keep)) {
    result[keep] <- RadixTree_find(S7::prop(x, "root_pointer"), query[keep])
  }
  result
}

S7::method(prefix_search, radix_tree) <- function(x, query) {
  check_character_vector(query, "query")
  keep <- not_na_character(query)
  if (!any(keep)) {
    return(seqtrie_empty_prefix_result())
  }
  result <- RadixTree_prefix_search(S7::prop(x, "root_pointer"), query[keep])
  if (is.null(result)) {
    seqtrie_empty_prefix_result()
  } else {
    result
  }
}

S7::method(align_search, radix_tree) <- function(x,
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
  args <- seqtrie_prepare_search(
    query = query,
    max_distance = max_distance,
    max_fraction = max_fraction,
    mode = mode,
    cost_matrix = cost_matrix,
    gap_cost = gap_cost,
    gap_open_cost = gap_open_cost,
    lower_triangle = lower_triangle,
    match_mode = match_mode,
    charset = CharCounter_keys(S7::prop(x, "char_counter_pointer")),
    allow_anchored = TRUE,
    nthreads = nthreads
  )
  if (!any(args$keep)) {
    return(seqtrie_empty_match_result())
  }

  result <- RadixTree_search(
    S7::prop(x, "root_pointer"),
    args$query_search,
    args$max_distance[args$keep],
    args$mode,
    cost_matrix,
    args$gap_cost,
    args$gap_open_cost,
    args$lower_triangle,
    args$match_mode,
    args$query_index,
    args$nthreads,
    show_progress
  )
  seqtrie_add_anchored_sizes(result, args$mode, cost_matrix, args$gap_cost, args$gap_open_cost, args$nthreads)
}

S7::method(single_gap_search, radix_tree) <- function(x,
                                                      query,
                                                      max_distance,
                                                      gap_cost = 1L,
                                                      nthreads = 1L,
                                                      show_progress = FALSE) {
  check_character_vector(query, "query")
  nthreads <- check_threads(nthreads)
  if (!is_integerlike(gap_cost) || length(gap_cost) != 1L || is.na(gap_cost) || gap_cost <= 0L) {
    stop("gap_cost must be a single positive integer")
  }
  gap_cost <- as.integer(gap_cost)

  keep <- not_na_character(query)
  if (!is.null(max_distance)) {
    max_distance <- recycle_arg(max_distance, query)
  } else {
    stop("Either max_distance must be non-null")
  }
  if (!is_integerlike(max_distance) || any(is.na(max_distance[keep])) || any(max_distance[keep] < 0)) {
    stop("max_distance/max_fraction must be non-negative integer-like values")
  }
  max_distance <- as.integer(max_distance)

  if (!any(keep)) {
    return(seqtrie_empty_match_result())
  }

  RadixTree_single_gap_search(S7::prop(x, "root_pointer"), query[keep], max_distance[keep], gap_cost, nthreads, show_progress)
}

S7::method(to_vector, radix_tree) <- function(x) {
  RadixTree_to_vector(S7::prop(x, "root_pointer"))
}

S7::method(to_string, radix_tree) <- function(x) {
  RadixTree_print(S7::prop(x, "root_pointer"))
}

S7::method(size, radix_tree) <- function(x) {
  RadixTree_size(S7::prop(x, "root_pointer"))
}

S7::method(plot_tree, radix_tree) <- function(x, depth = -1, root_label = "root", plot = TRUE) {
  seqtrie_plot_graph(RadixTree_graph(S7::prop(x, "root_pointer"), depth), root_label, plot)
}

S7::method(is_valid, radix_tree) <- function(x) {
  RadixTree_validate(S7::prop(x, "root_pointer"))
}

# Radix forest methods ----------------------------------------------------------

S7::method(insert, radix_forest) <- function(x, sequences) {
  check_character_vector(sequences, "sequences")
  keep <- not_na_character(sequences)
  result <- rep(FALSE, length(sequences))
  if (any(keep)) {
    result[keep] <- RadixForest_insert(S7::prop(x, "forest_pointer"), sequences[keep])
    CharCounter_add(S7::prop(x, "char_counter_pointer"), sequences[keep][result[keep]])
  }
  invisible(result)
}

S7::method(erase, radix_forest) <- function(x, sequences) {
  check_character_vector(sequences, "sequences")
  keep <- not_na_character(sequences)
  result <- rep(FALSE, length(sequences))
  if (any(keep)) {
    result[keep] <- RadixForest_erase(S7::prop(x, "forest_pointer"), sequences[keep])
    CharCounter_subtract(S7::prop(x, "char_counter_pointer"), sequences[keep][result[keep]])
  }
  invisible(result)
}

S7::method(has_sequence, radix_forest) <- function(x, query) {
  check_character_vector(query, "query")
  keep <- not_na_character(query)
  result <- rep(FALSE, length(query))
  if (any(keep)) {
    result[keep] <- RadixForest_find(S7::prop(x, "forest_pointer"), query[keep])
  }
  result
}

S7::method(prefix_search, radix_forest) <- function(x, query) {
  check_character_vector(query, "query")
  keep <- not_na_character(query)
  if (!any(keep)) {
    return(seqtrie_empty_prefix_result())
  }
  result <- RadixForest_prefix_search(S7::prop(x, "forest_pointer"), query[keep])
  if (is.null(result)) {
    seqtrie_empty_prefix_result()
  } else {
    result
  }
}

S7::method(align_search, radix_forest) <- function(x,
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
  args <- seqtrie_prepare_search(
    query = query,
    max_distance = max_distance,
    max_fraction = max_fraction,
    mode = mode,
    cost_matrix = cost_matrix,
    gap_cost = gap_cost,
    gap_open_cost = gap_open_cost,
    lower_triangle = lower_triangle,
    match_mode = match_mode,
    charset = CharCounter_keys(S7::prop(x, "char_counter_pointer")),
    allow_anchored = FALSE,
    nthreads = nthreads
  )
  if (!any(args$keep)) {
    return(seqtrie_empty_match_result())
  }

  RadixForest_search(
    S7::prop(x, "forest_pointer"),
    args$query_search,
    args$max_distance[args$keep],
    args$mode,
    cost_matrix,
    args$gap_cost,
    args$gap_open_cost,
    args$lower_triangle,
    args$match_mode,
    args$query_index,
    args$nthreads,
    show_progress
  )
}

S7::method(to_vector, radix_forest) <- function(x) {
  RadixForest_to_vector(S7::prop(x, "forest_pointer"))
}

S7::method(to_string, radix_forest) <- function(x) {
  RadixForest_print(S7::prop(x, "forest_pointer"))
}

S7::method(size, radix_forest) <- function(x) {
  RadixForest_size(S7::prop(x, "forest_pointer"))
}

S7::method(plot_tree, radix_forest) <- function(x, depth = -1, root_label = "root", plot = TRUE) {
  seqtrie_plot_graph(RadixForest_graph(S7::prop(x, "forest_pointer"), depth), root_label, plot)
}

S7::method(is_valid, radix_forest) <- function(x) {
  RadixForest_validate(S7::prop(x, "forest_pointer"))
}
