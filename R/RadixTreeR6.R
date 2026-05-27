#' @title RadixTree
#' @description R6 compatibility wrapper for `radix_tree`
#'
#' @details
#' `RadixTree` preserves the original R6 API while delegating implementation to
#' the S7 [radix_tree] class. New code can use `radix_tree()` with the exported
#' S7 generics directly; existing code using `$insert()`, `$erase()`, `$find()`,
#' `$prefix_search()`, `$search()`, `$single_gap_search()`, `$to_vector()`,
#' `$to_string()`, `$size()`, `$graph()`, and `$validate()` remains supported.
#'
#' The RadixTree implementation stores sequences in a trie and searches for
#' similar sequences with hamming, global/Levenshtein, anchored, or single-gap
#' alignment metrics.
#'
#' `r rdoc("details")`
#' @examples
#' tree <- RadixTree$new()
#' tree$insert(c("ACGT", "AAAA"))
#' tree$erase("AAAA")
#' tree$search("ACG", max_distance = 1, mode = "levenshtein")
#' #   query target distance
#' # 1   ACG   ACGT        1
#'
#' tree$search("ACG", max_distance = 1, mode = "hamming")
#' # query    target   distance
#' # <0 rows> (or 0-length row.names)
#' @seealso [radix_tree], <https://en.wikipedia.org/wiki/Radix_tree>
#' @export
RadixTree <- R6::R6Class("RadixTree", public = list(
  #' @field root_pointer Root of the RadixTree.
  root_pointer = NULL,
  #' @field char_counter_pointer Character count data for validating input.
  char_counter_pointer = NULL,

  #' @description Create a new RadixTree object.
  #' @param sequences A character vector of sequences to insert into the tree.
  initialize = function(sequences = NULL) {
    private$s7_object <- radix_tree(sequences)
    self$root_pointer <- S7::prop(private$s7_object, "root_pointer")
    self$char_counter_pointer <- S7::prop(private$s7_object, "char_counter_pointer")
  },

  #' @description Print the tree to screen.
  show = function() {
    cat(to_string(private$s7_object))
  },

  #' @description Print the tree to a string.
  #' @return A string representation of the tree.
  to_string = function() {
    to_string(private$s7_object)
  },

  #' @description Plot the tree using igraph and ggplot2.
  #' @param depth The tree depth to plot. If -1, plot the entire tree.
  #' @param root_label The label of the root node in the plot.
  #' @param plot Whether to create a plot or return the graph data.
  #' @return A data frame of parent-child relationships or a ggplot2 object.
  graph = function(depth = -1, root_label = "root", plot = TRUE) {
    plot_tree(private$s7_object, depth = depth, root_label = root_label, plot = plot)
  },

  #' @description Output all stored sequences as a character vector.
  #' @return A character vector of all sequences contained in the tree.
  to_vector = function() {
    to_vector(private$s7_object)
  },

  #' @description Output the size of the tree.
  #' @return The number of stored sequences.
  size = function() {
    size(private$s7_object)
  },

  #' @description Insert new sequences into the tree.
  #' @param sequences A character vector of sequences to insert into the tree.
  #' @return A logical vector indicating whether each sequence was inserted.
  insert = function(sequences) {
    insert(private$s7_object, sequences)
  },

  #' @description Erase sequences from the tree.
  #' @param sequences A character vector of sequences to erase from the tree.
  #' @return A logical vector indicating whether each sequence was erased.
  erase = function(sequences) {
    erase(private$s7_object, sequences)
  },

  #' @description Find sequences in the tree.
  #' @param query A character vector of sequences to find in the tree.
  #' @return A logical vector indicating whether each sequence was found.
  find = function(query) {
    has_sequence(private$s7_object, query)
  },

  #' @description Find sequences in the tree.
  #' @param query A character vector of sequences to find in the tree.
  #' @return A logical vector indicating whether each sequence was found.
  has_sequence = function(query) {
    has_sequence(private$s7_object, query)
  },

  #' @description Search for sequences in the tree that start with a prefix.
  #' @param query A character vector of prefixes.
  #' @return A data frame with columns `query` and `target`.
  prefix_search = function(query) {
    prefix_search(private$s7_object, query)
  },

  #' @description Alignment search.
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
  #' @return A data frame with columns `query`, `target`, and `distance`.
  search = function(query,
                    max_distance = NULL,
                    max_fraction = NULL,
                    mode = "levenshtein",
                    cost_matrix = NULL,
                    gap_cost = NA_integer_,
                    gap_open_cost = NA_integer_,
                    lower_triangle = FALSE,
                    match_mode = c("all", "best"),
                    nthreads = 1,
                    show_progress = FALSE) {
    align_search(
      private$s7_object,
      query = query,
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
  },

  #' @description Alignment search.
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
  #' @return A data frame with columns `query`, `target`, and `distance`.
  align_search = function(query,
                          max_distance = NULL,
                          max_fraction = NULL,
                          mode = "levenshtein",
                          cost_matrix = NULL,
                          gap_cost = NA_integer_,
                          gap_open_cost = NA_integer_,
                          lower_triangle = FALSE,
                          match_mode = c("all", "best"),
                          nthreads = 1,
                          show_progress = FALSE) {
    align_search(
      private$s7_object,
      query = query,
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
  },

  #' @description A specialized anchored search allowing at most one internal gap.
  #' @param query `r rdoc("query")`
  #' @param max_distance `r rdoc("max_distance")`
  #' @param gap_cost `r rdoc("gap_cost")`
  #' @param nthreads `r rdoc("nthreads")`
  #' @param show_progress `r rdoc("show_progress")`
  #' @return A data frame with columns `query`, `target`, and `distance`.
  single_gap_search = function(query,
                               max_distance,
                               gap_cost = 1L,
                               nthreads = 1,
                               show_progress = FALSE) {
    single_gap_search(
      private$s7_object,
      query = query,
      max_distance = max_distance,
      gap_cost = gap_cost,
      nthreads = nthreads,
      show_progress = show_progress
    )
  },

  #' @description Validate the tree.
  #' @return A logical value indicating whether the tree is valid.
  validate = function() {
    is_valid(private$s7_object)
  },

  #' @description Validate the tree.
  #' @return A logical value indicating whether the tree is valid.
  is_valid = function() {
    is_valid(private$s7_object)
  }
),
private = list(
  s7_object = NULL
),
cloneable = FALSE)
