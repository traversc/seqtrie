#' @title StarTree
#' @description R6 compatibility wrapper for `star_tree`
#'
#' @details
#' `StarTree` is a fixed DNA-only tree using a modified version of the Starcode
#' all-pairs search strategy, adapted to operate over a radix trie rather than
#' being a direct reimplementation. It supports global/Levenshtein, anchored,
#' and Hamming alignment modes. Unlike [RadixTree], all sequences and alignment parameters
#' are supplied at construction time and the self-similarity join runs
#' immediately. The tree does not support insertion or deletion after
#' construction.
#'
#' Use `$result()` to retrieve the construction-time self-similarity join, and
#' `$align_search()` or `$search()` to search additional query sequences against
#' the fixed target set using the same `mode`, `max_distance`, `mismatch_cost`,
#' and `gap_cost`.
#'
#' The algorithm is based on Starcode (Zorita, Cuscó, and Filion 2015)
#' \doi{10.1093/bioinformatics/btv053}.
#'
#' @examples
#' tree <- StarTree$new(c("ACGT", "ACGA", "AAAA"), max_distance = 1)
#' tree$result()
#' tree$search(c("ACGT", "AAAT"))
#'
#' anchored <- StarTree$new(c("ACGT", "ACG", "AAAA"), max_distance = 1,
#'                          mode = "anchored")
#' anchored$result()
#' @seealso [align_search()]
#' @export
StarTree <- R6::R6Class("StarTree", public = list(
  #' @field tree_pointer External pointer to the fixed-tree C++ object.
  tree_pointer = NULL,
  #' @field mode Fixed alignment mode.
  mode = NULL,
  #' @field max_distance Fixed distance threshold.
  max_distance = NULL,
  #' @field mismatch_cost Fixed mismatch cost.
  mismatch_cost = NULL,
  #' @field gap_cost Fixed gap cost.
  gap_cost = NULL,
  #' @field nthreads Fixed thread count.
  nthreads = NULL,
  #' @field show_progress Fixed progress flag.
  show_progress = NULL,

  #' @description Create a new StarTree object.
  #' @param sequences A required character vector of DNA sequences.
  #' @param max_distance A single non-negative integer distance threshold.
  #' @param mode Alignment mode: global/Levenshtein, anchored, or hamming.
  #' @param mismatch_cost A single positive integer mismatch cost.
  #' @param gap_cost A single positive integer gap cost.
  #' @param nthreads `r rdoc("nthreads")`
  #' @param show_progress `r rdoc("show_progress")`
  initialize = function(sequences,
                        max_distance,
                        mode = "levenshtein",
                        mismatch_cost = 1L,
                        gap_cost = 1L,
                        nthreads = 1L,
                        show_progress = FALSE) {
    private$s7_object <- star_tree(
      sequences = sequences,
      max_distance = max_distance,
      mode = mode,
      mismatch_cost = mismatch_cost,
      gap_cost = gap_cost,
      nthreads = nthreads,
      show_progress = show_progress
    )
    self$tree_pointer <- S7::prop(private$s7_object, "tree_pointer")
    self$mode <- S7::prop(private$s7_object, "mode")
    self$max_distance <- S7::prop(private$s7_object, "max_distance")
    self$mismatch_cost <- S7::prop(private$s7_object, "mismatch_cost")
    self$gap_cost <- S7::prop(private$s7_object, "gap_cost")
    self$nthreads <- S7::prop(private$s7_object, "nthreads")
    self$show_progress <- S7::prop(private$s7_object, "show_progress")
  },

  #' @description Output all stored unique sequences as a character vector.
  #' @return A character vector of all sequences contained in the tree.
  to_vector = function() {
    to_vector(private$s7_object)
  },

  #' @description Output the size of the tree.
  #' @return The number of stored unique sequences.
  size = function() {
    size(private$s7_object)
  },

  #' @description Return the construction-time self-similarity join.
  #' @return A data frame with columns `query`, `target`, and `distance`.
  #' Anchored mode also includes `query_size` and `target_size`.
  result = function() {
    result(private$s7_object)
  },

  #' @description Search additional query sequences against the fixed tree.
  #' @param query `r rdoc("query")`
  #' @param nthreads `r rdoc("nthreads")`
  #' @param show_progress `r rdoc("show_progress")`
  #' @return A data frame with columns `query`, `target`, and `distance`.
  #' Anchored mode also includes `query_size` and `target_size`.
  search = function(query,
                    nthreads = self$nthreads,
                    show_progress = self$show_progress) {
    align_search(
      private$s7_object,
      query = query,
      nthreads = nthreads,
      show_progress = show_progress
    )
  },

  #' @description Search additional query sequences against the fixed tree.
  #' @param query `r rdoc("query")`
  #' @param nthreads `r rdoc("nthreads")`
  #' @param show_progress `r rdoc("show_progress")`
  #' @return A data frame with columns `query`, `target`, and `distance`.
  #' Anchored mode also includes `query_size` and `target_size`.
  align_search = function(query,
                          nthreads = self$nthreads,
                          show_progress = self$show_progress) {
    align_search(
      private$s7_object,
      query = query,
      nthreads = nthreads,
      show_progress = show_progress
    )
  }
),
private = list(
  s7_object = NULL
),
cloneable = FALSE)
