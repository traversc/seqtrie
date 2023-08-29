#' @title RadixForest
#' @description Radix Forest class implementation
#' 
#' @details
#' The RadixForest class is a specialization of the RadixTree implementation. Instead of putting sequences into a single tree, the RadixForest class puts sequences into a tree based on sequence length. 
#' I.e. *map<sequence_length, RadixTree>. This allows for faster searching of similar sequences based on Hamming or Levenshtein distance metrics.
#' Unlike the RadixTree class, the RadixForest class does NOT support Anchored searches or a custom cost matrix.
#' See *RadixTree* for additional details. 
#' @examples
#' forest <- RadixForest$new()
#' forest$insert(c("ACGT", "AAAA"))
#' forest$erase("AAAA")
#' forest$search("ACG", max_distance = 1, mode = "levenshtein")
#'  #   query target distance
#'  # 1   ACG   ACGT        1
#'  
#' forest$search("ACG", max_distance = 1, mode = "hamming")
#'  # query    target   distance
#'  # <0 rows> (or 0-length row.names)
RadixForest <- R6::R6Class("RadixForest", list(
  #' @field forest_pointer Pointer C++ implementation (holds trie map)
  forest_pointer = NULL,
  #' @field char_counter_pointer Pointer to C++ object (holds character counts for the purpose of validating input)
  char_counter_pointer = NULL,
  #' @description Create a new RadixForest object
  #' @param sequences A character vector of sequences to insert into the forest
  initialize = function(sequences = NULL) {
    self$forest_pointer <- RadixForest_create()
    if(!is.null(sequences)) {
      self$insert(sequences)
    }
  },
  #' @description Print the forest to screen
  show = function() {
    cat(RadixForest_print(self$forest_pointer))
  },
  #' @description Print the forest to a string
  to_string = function() {
    RadixForest_print(self$forest_pointer)
  },
  #' @description Plot of the forest using igraph
  #' @param depth The tree depth to plot for each tree in the forest.
  #' @param root_label The label of the root node(s) in the plot.
  #' @param plot Whether to create a plot or return the data used to generate the plot.
  #' @return A data frame of parent-child relationships used to generate the igraph plot OR a ggplot2 object
  graph = function(depth = -1, root_label = "root", plot = TRUE) {
    result <- RadixForest_graph(self$forest_pointer, depth)
    if(is.null(result)) {
      result <- data.frame(parent = character(0), child = character(0), stringsAsFactors=FALSE)
      return(result)
    } else if(plot) {
      if (!requireNamespace("igraph", quietly = TRUE)) {
        stopf("igraph package is required to plot the tree.") # nocov
      }
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stopf("ggplot2 package is required to plot the tree.") # nocov
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
        ggplot2::geom_segment(data=result, ggplot2::aes(x=parent_x, xend=child_x, y=parent_y, yend=child_y)) + 
        ggplot2::geom_point(data=fr, ggplot2::aes(x=V1, y=V2, fill=fill, size=size), shape = 21, color = "black") + 
        ggplot2::geom_text(data=fr, ggplot2::aes(x=V1, y=V2, label=node)) + 
        ggplot2::scale_fill_identity() + 
        ggplot2::scale_size_identity() + 
        ggplot2::theme_bw() + 
        ggplot2::theme(axis.title = ggplot2::element_blank())
      return(g)
    } else {
      return(result)
    }
  },
  #' @description Output all sequences held by the forest as a character vector
  #' @return A character vector of all sequences contained in the forest.
  to_vector = function() {
    RadixForest_to_vector(self$forest_pointer)
  },
  #' @description Output the size of the forest (i.e. how many sequences are contained)
  #' @return The size of the forest
  size = function() {
    RadixForest_size(self$forest_pointer)
  },
  #' @description Insert new sequences into the forest
  #' @param sequences A character vector of sequences to insert into the forest
  #' @return A logical vector indicating whether the sequence was inserted (TRUE) or already existing in the forest (FALSE)
  insert = function(sequences) {
    result <- RadixForest_insert(self$forest_pointer, sequences)
    invisible(result)
  },
  #' @description Erase sequences from the forest
  #' @param sequences A character vector of sequences to erase from the forest
  #' @return A logical vector indicating whether the sequence was erased (TRUE) or not found in the forest (FALSE)
  erase = function(sequences) {
    result <- RadixForest_erase(self$forest_pointer, sequences)
    invisible(result)
  },
  #' @description Find sequences in the forest
  #' @param query A character vector of sequences to find in the forest
  #' @return A logical vector indicating whether the sequence was found (TRUE) or not found in the forest (FALSE)
  find = function(query) {
    RadixForest_find(self$forest_pointer, query)
  },
  #' @description Search for sequences in the forest that start with a specified prefix.
  #' E.g.: a query of "CAR" will find "CART", "CARBON", "CARROT", etc. but not "CATS".
  #' @param query A character vector of sequences to search for in the forest
  #' @return A data frame of all matches with columns "query" and "target".
  prefix_search = function(query) {
    result <- RadixForest_prefix_search(self$forest_pointer, query)
    if(is.null(result)) {
      data.frame(query = character(0), target = character(0), stringsAsFactors=FALSE)
    } else {
      result
    }
  },
  #' @description Search for sequences in the forest that are with a specified distance metric to a specified query.
  #' @param query `r rdoc("query")`
  #' @param max_distance `r rdoc("max_distance")`
  #' @param max_fraction `r rdoc("max_fraction")`
  #' @param mode `r rdoc("mode")`
  #' @param nthreads `r rdoc("nthreads")`
  #' @param show_progress `r rdoc("show_progress")`
  #' @return The output is a data.frame of all matches with columns "query" and "target".
  
  search = function(query, max_distance = NULL, max_fraction = NULL, mode = "levenshtein", nthreads = 1, show_progress = FALSE) {
    check_alignment_params(mode, cost_matrix=NULL, gap_cost=NULL, gap_open_cost=NULL, charset = "", diag_must_be_zero = TRUE)
    mode <- normalize_mode_parameter(mode)
    if(!mode %in% c("hamming", "global")) {
      stop("mode must be one of hamming (hm) or global (gb, lv, levenshtein)")
    }
    
    if(!is.null(max_distance)) {
      if(length(max_distance) == 1) {
        max_distance <- rep(max_distance, length(query))
      }
    } else if(!is.null(max_fraction)) {
      max_distance <- as.integer(nchar(query) * max_fraction)
    } else {
      stop("Either max_distance or max_fraction must be non-null")
    }
    if(any(max_distance < 0)) {
      stop("max_distance/max_fraction must be non-negative")
    }
    RadixForest_search(self$forest_pointer, query, max_distance, mode, nthreads, show_progress)
  },
  #' @description Validate the forest
  #' @return A logical indicating whether the forest is valid (TRUE) or not (FALSE). This is mostly an internal function for debugging purposes and should always return TRUE.
  validate = function() {
    RadixForest_validate(self$forest_pointer)
  }
))

