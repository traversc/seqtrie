#' @title RadixTree
#' @description Radix Tree (trie) class implementation
#'
#' @details
#' The RadixTree class is a trie implementation. The primary usage is to be able to search of similar sequences based on a dynamic programming framework.
#' This can be done using the _search_ method which searches for similar sequences based on the Global, Anchored or Hamming distance metrics.
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
#' @seealso
#' https://en.wikipedia.org/wiki/Radix_tree
RadixTree <- R6::R6Class("RadixTree", public = list(
  #' @field root_pointer Root of the RadixTree
  root_pointer = NULL,
  #' @field char_counter_pointer Character count data for the purpose of validating input
  char_counter_pointer = NULL,
  #' @description Create a new RadixTree object
  #' @param sequences A character vector of sequences to insert into the tree
  initialize = function(sequences = NULL) {
    self$root_pointer <- RadixTree_create()
    self$char_counter_pointer <- CharCounter_create()
    if (!is.null(sequences)) {
      self$insert(sequences)
    }
  },
  #' @description Print the tree to screen
  show = function() {
    cat(RadixTree_print(self$root_pointer))
  },
  #' @description Print the tree to a string
  #' @return A string representation of the tree
  to_string = function() {
    RadixTree_print(self$root_pointer)
  },
  #' @description Plot of the tree using igraph (needs to be installed separately)
  #' @param depth The tree depth to plot. If -1 (default), plot the entire tree.
  #' @param root_label The label of the root node in the plot.
  #' @param plot Whether to create a plot or return the data used to generate the plot.
  #' @return A data frame of parent-child relationships used to generate the igraph plot OR a ggplot2 object
  graph = function(depth = -1, root_label = "root", plot = TRUE) {
    result <- RadixTree_graph(self$root_pointer, depth)
    if (is.null(result)) {
      result <- data.frame(parent = character(0), child = character(0), stringsAsFactors = FALSE)
      return(result)
    } else if (plot) {
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

      g <- ggplot2::ggplot() +
        ggplot2::geom_segment(data = result, ggplot2::aes(x = parent_x, xend = child_x, y = parent_y, yend = child_y)) +
        ggplot2::geom_point(data = fr, ggplot2::aes(x = V1, y = V2, fill = fill, size = size), shape = 21, color = "black") +
        ggplot2::geom_text(data = fr, ggplot2::aes(x = V1, y = V2, label = node)) +
        ggplot2::scale_fill_identity() +
        ggplot2::scale_size_identity() +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.title = ggplot2::element_blank())
      return(g)
    } else {
      return(result)
    }
  },
  #' @description Output all sequences held by the tree as a character vector
  #' @return A character vector of all sequences contained in the tree. Return order is not guaranteed.
  to_vector = function() {
    RadixTree_to_vector(self$root_pointer)
  },
  #' @description Output the size of the tree (i.e. how many sequences are contained)
  #' @return The size of the tree
  size = function() {
    RadixTree_size(self$root_pointer)
  },
  #' @description Insert new sequences into the tree
  #' @param sequences A character vector of sequences to insert into the tree
  #' @return A logical vector indicating whether the sequence was inserted (TRUE) or already existing in the tree (FALSE)
  insert = function(sequences) {
    result <- RadixTree_insert(self$root_pointer, sequences)
    CharCounter_add(self$char_counter_pointer, sequences[result])
    invisible(result)
  },
  #' @description Erase sequences from the tree
  #' @param sequences A character vector of sequences to erase from the tree
  #' @return A logical vector indicating whether the sequence was erased (TRUE) or not found in the tree (FALSE)
  erase = function(sequences) {
    result <- RadixTree_erase(self$root_pointer, sequences)
    CharCounter_subtract(self$char_counter_pointer, sequences[result])
    invisible(result)
  },
  #' @description Find sequences in the tree
  #' @param query A character vector of sequences to find in the tree
  #' @return A logical vector indicating whether the sequence was found (TRUE) or not found in the tree (FALSE)
  find = function(query) {
    RadixTree_find(self$root_pointer, query)
  },
  #' @description Search for sequences in the tree that start with a specified prefix.
  #' E.g.: a query of "CAR" will find "CART", "CARBON", "CARROT", etc. but not "CATS".
  #' @param query A character vector of sequences to search for in the tree
  #' @return A data frame of all matches with columns "query" and "target".
  prefix_search = function(query) {
    result <- RadixTree_prefix_search(self$root_pointer, query)
    if (is.null(result)) {
      data.frame(query = character(0), target = character(0), stringsAsFactors = FALSE)
    } else {
      result
    }
  },
  #' @description Search for sequences in the tree that are with a specified distance metric to a specified query.
  #' @param query `r rdoc("query")`
  #' @param max_distance `r rdoc("max_distance")`
  #' @param max_fraction `r rdoc("max_fraction")`
  #' @param mode `r rdoc("mode")`
  #' @param cost_matrix `r rdoc("cost_matrix")`
  #' @param gap_cost `r rdoc("gap_cost")`
  #' @param gap_open_cost `r rdoc("gap_open_cost")`
  #' @param nthreads `r rdoc("nthreads")`
  #' @param show_progress `r rdoc("show_progress")`
  #' @return The output is a data.frame of all matches with columns "query" and "target".
  #' For anchored searches, the output also includes attributes "query_size" and "target_size" which are vectors containing the portion of the query and target sequences that are aligned.
  search = function(query, max_distance = NULL, max_fraction = NULL, mode = "levenshtein", cost_matrix = NULL, gap_cost = NULL, gap_open_cost = NULL, nthreads = 1, show_progress = FALSE) {
    charset <- unique(c(CharCounter_keys(self$char_counter_pointer), get_charset(query)))
    check_alignment_params(mode, cost_matrix, gap_cost, gap_open_cost, charset, diag_must_be_zero = TRUE)
    mode <- normalize_mode_parameter(mode)

    if (!is.null(max_distance)) {
      if (length(max_distance) == 1) {
        max_distance <- rep(max_distance, length(query))
      }
    } else if (!is.null(max_fraction)) {
      max_distance <- as.integer(nchar(query) * max_fraction)
    } else {
      stop("Either max_distance or max_fraction must be non-null")
    }
    if (any(max_distance < 0)) {
      stop("max_distance/max_fraction must be non-negative")
    }

    # defaults for C++ plain ints
    if (is.null(gap_cost)) gap_cost <- 1L
    if (is.null(gap_open_cost)) gap_open_cost <- 0L
    # Align conventions with pwalign/Biostrings: first gap includes one extension
    if (gap_open_cost > 0L) gap_open_cost <- gap_open_cost + gap_cost
    result <- RadixTree_search(self$root_pointer, query, max_distance, mode, cost_matrix, as.integer(gap_cost), as.integer(gap_open_cost), nthreads, show_progress)
    if (mode == "anchored") { # Append query_size and target_size attributes
      result2 <- c_dist_pairwise(result$query, result$target, mode, cost_matrix, as.integer(gap_cost), as.integer(gap_open_cost), nthreads, show_progress = FALSE)
      if (any(result$distance != result2)) {
        stop("Internal error: anchored search results do not match pairwise results")
      }
      result$query_size <- attr(result2, "query_size")
      result$target_size <- attr(result2, "target_size")
    }
    result
  },
  #' @description Validate the tree
  #' @return A logical indicating whether the tree is valid (TRUE) or not (FALSE). This is mostly an internal function for debugging purposes and should always return TRUE.
  validate = function() {
    RadixTree_validate(self$root_pointer)
  }
),
cloneable=FALSE)
