#' @title RadixTree
#' @description Radix Tree (trie) class implementation
#' 
#' @details
#' The RadixTree class is a trie implementation. The primary usage is to be able to search of similar sequences based on a dynamic programming framework. 
#' This can be done using the _search_ method which searches for similar sequences based on the Global, Anchored or Hamming distance metrics.
#' `r rdoc("details")`
#' @examples
#' tree <- RadixTree$new()
#' tree$insert(c("ACGT", "AAAA"))
#' tree$erase("AAAA")
#' tree$search("ACG", max_distance = 1, mode = "levenshtein")
#'  #   query target distance
#'  # 1   ACG   ACGT        1
#'  
#' tree$search("ACG", max_distance = 1, mode = "hamming")
#'  # query    target   distance
#'  # <0 rows> (or 0-length row.names)
#' @seealso 
#' https://en.wikipedia.org/wiki/Radix_tree
RadixTree <- R6::R6Class("RadixTree", public=list(
  #' @field root_pointer Pointer C++ implementation (holds trie map)
  root_pointer = NULL,
  #' @field char_counter_pointer Pointer to C++ object (holds character counts for the purpose of validating input)
  char_counter_pointer = NULL,
  #' @description Create a new RadixTree object
  #' @param sequences A character vector of sequences to insert into the tree
  initialize = function(sequences = NULL) {
    self$root_pointer <- RadixTree_create()
    self$char_counter_pointer <- CharCounter_create()
    if(!is.null(sequences)) {
      self$insert(sequences)
    }
  },
  #' @description Print the tree to screen
  print = function() {
    cat(RadixTree_print(self$root_pointer))
  },
  #' @description Print the tree to a string
  #' @return A string representation of the tree
  to_string = function() {
    RadixTree_print(self$root_pointer)
  },
  #' @description Plot of the tree using igraph
  #' @param depth The tree depth to plot. If -1 (default), plot the entire tree.
  #' @param root_label The label of the root node in the plot. 
  #' @param plot Whether to create a plot or return the data used to generate the plot.
  #' @return a data frame of parent-child relationships used to generate the igraph plot
  graph = function(depth = -1, root_label = "root", plot = TRUE) {
    result <- RadixTree_graph(self$root_pointer, depth)
    if(is.null(result)) {
      result <- data.frame(parent = character(0), child = character(0), stringsAsFactors=F)
    } else if(plot) {
      result$parent <- ifelse(result$parent == "", root_label, result$parent)
      gr <- igraph::graph_from_data_frame(result)
      igraph::V(gr)$color <- ifelse(names(igraph::V(gr)) == root_label, "white", "skyblue2")
      igraph::V(gr)$size <- ifelse(names(igraph::V(gr)) == root_label, 21, 15)
      plot(gr, layout=igraph::layout.fruchterman.reingold, vertex.color=igraph::V(gr)$color, vertex.label.family = "sans", margin = 0) 
    }
    invisible(result)
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
    if(is.null(result)) {
      data.frame(query = character(0), target = character(0), stringsAsFactors=FALSE)
    } else {
      result
    }
  },
  #' @description Search for sequences in the tree that are with a specified distance metric to a specified query.
  #' @param query `r rdoc("query")`
  #' @param max_distance how far to search in units of absolute distance. Can be a single value or a vector. Mutually exclusive with max_fraction.
  #' @param max_fraction how far to search in units of relative distance to each query sequence length. Can be a single value or a vector. Mutually exclusive with max_distance.
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
    finalized_cost_matrix <- finalize_cost_matrix(cost_matrix, gap_cost, gap_open_cost)
    mode <- normalize_mode_parameter(mode)
    gap_type <- get_gap_type(finalized_cost_matrix)
    
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

    result <- RadixTree_search(self$root_pointer, query, max_distance, mode, gap_type, finalized_cost_matrix, nthreads, show_progress)
    if(mode == "anchored") { # Append query_size and target_size attributes
      result2 <- c_dist_pairwise(result$query, result$target, mode, gap_type = gap_type, finalized_cost_matrix, nthreads, show_progress = FALSE)
      if(any(result$distance != result2)) {
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
))

