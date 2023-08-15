#' RadixForest
#'
#' Radix Forest (trie) class implementation
#'
#' @section Usage:
#' \preformatted{tree <- RadixForest$new()
#' 
#' tree$print()
#' 
#' tree$graph(depth = -1, root_label = "root", plot = TRUE)
#' 
#' tree$to_string()
#' 
#' tree$to_vector()
#' 
#' tree$size()
#' 
#' tree$insert(sequences)
#' 
#' tree$erase(sequences)
#' 
#' tree$find(sequences)
#' 
#' tree$search(sequences, max_distance = NULL, max_fraction = NULL, mode = "levenshtein", cost_matrix = NULL, nthreads = 1, show_progress = FALSE)
#' }
#'
#' @section Arguments:
#' \describe{
#'   \item{depth}{- In \code{$graph()}, the tree depth to plot. Default -1 means plot the entire tree.}
#'   \item{root_label}{- In \code{$graph()}, the label of the root node.}
#'   \item{plot}{- In \code{$graph()}, whether to create a plot or return the data used to generate the plot.}
#'   \item{sequences}{- In \code{$search()}, the sequences to operate on.}
#'   \item{max_distance}{- In \code{$search()}, how far to search for similar sequences, in units of absolute distance. See details.}
#'   \item{max_fraction}{- In \code{$search()}, how far to search for similar sequences, relative to the sequence length. See details.}
#'   \item{mode}{- In \code{$search()}, either hamming (hm) or levenshtein (lv). Levenshtein will allows for insertions and deletions and calculates "edit distance". Hamming does not allow for insertions or deletions.}
#'   \item{nthreads}{- How many threads to use in the search.}
#'   \item{show_progress}{- Display progress.}
#' }
#'
#' @section Details:
#'
#' \code{$new()} creates a new Radix Forest object, which holds the pointer to the underlying C++ implentation. The RadixForest class accepts any strings of single-width characters. 
#'
#' \code{$print()} and \code{$to_string()} prints to screen or outputs the tree to a string representation.
#' 
#' \code{$to_vector()} outputs all sequences held by the tree. 
#' 
#' \code{$size()} outputs the size of the tree (i.e. how many sequences are contained). 
#' 
#' \code{$insert()}, \code{$erase()} and \code{$find()} insert, erase and find sequences in the tree, respectively.
#' 
#' \code{$search()} This function searches for similar sequences within a threshold (given by max_distance or max_fraction) based on Hamming or Levenshtein. Unlike the `RadixTree` class, 
#' `RadixForest` does not support anchored searches or custom cost matrices. 
#' 
#' The output of this function is a data.frame of all matches with columns "query" (the sequences input to the search function), 
#' "target" (the sequences inserted into the tree) and "distance" the absolute distance between query and target sequences. 
#' 
#' @seealso 
#' https://en.wikipedia.org/wiki/Radix_tree
#' 
#' @examples
#' tree <- RadixForest$new()
#' tree$insert(c("ACGT", "AAAA"))
#' tree$erase("AAAA")
#' tree$search("ACG", max_distance = 1, mode = "levenshtein")
#'  #   query target distance
#'  # 1   ACG   ACGT        1
#'  
#' tree$search("ACG", max_distance = 1, mode = "hamming")
#'  # query    target   distance
#'  # <0 rows> (or 0-length row.names)
#' @name RadixForest
NULL


RadixForest <- R6::R6Class("RadixForest", list(
  forest_pointer = NULL,
  initialize = function(sequences = NULL) {
    self$forest_pointer <- RadixForest_create()
    if(!is.null(sequences)) {
      RadixForest_insert(self$forest_pointer, sequences)
    }
  },
  print = function() {
    cat(RadixForest_print(self$forest_pointer))
  },
  to_string = function() {
    RadixForest_print(self$forest_pointer)
  },
  graph = function(depth = -1, root_label = "root", plot = TRUE) {
    result <- RadixForest_graph(self$forest_pointer, depth)
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
  to_vector = function() {
    RadixForest_to_vector(self$forest_pointer)
  },
  size = function() {
    RadixForest_size(self$forest_pointer)
  },
  insert = function(sequences) {
    result <- RadixForest_insert(self$forest_pointer, sequences)
    invisible(result)
  },
  erase = function(sequences) {
    result <- RadixForest_erase(self$forest_pointer, sequences)
    invisible(result)
  },
  find = function(sequences) {
    RadixForest_find(self$forest_pointer, sequences)
  },
  prefix_search = function(sequences) {
    result <- RadixForest_prefix_search(self$forest_pointer, sequences)
    if(is.null(result)) {
      data.frame(query = character(0), target = character(0), stringsAsFactors=F)
    } else {
      result
    }
  },
  search = function(sequences, max_distance = NULL, max_fraction = NULL, mode = "levenshtein", nthreads = 1, show_progress = FALSE) {
    check_alignment_params(mode, cost_matrix=NULL, gap_cost=NULL, charset = "")
    stopifnot(mode %in% c("hamming", "levenshtein","hm", "lv"))
    mode <- normalize_mode_parameter(mode)
    
    if(!is.null(max_distance)) {
      if(length(max_distance) == 1) {
        max_distance <- rep(max_distance, length(sequences))
      }
    } else if(!is.null(max_fraction)) {
      max_distance <- as.integer(nchar(sequences) * max_fraction)
    } else {
      stop("Either max_distance or max_fraction must be non-null")
    }
    if(any(max_distance < 0)) {
      stop("max_distance/max_fraction must be non-negative")
    }
    if(mode == "hamming") {
      RadixForest_hamming_search(self$forest_pointer, sequences, max_distance, nthreads, show_progress)
    } else if(mode == "levenshtein") {
      RadixForest_levenshtein_search(self$forest_pointer, sequences, max_distance, nthreads, show_progress)
    }
  },
  validate = function() {
    RadixForest_validate(self$forest_pointer)
  }
))

