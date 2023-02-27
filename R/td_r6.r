#' RadixTree
#'
#' This class is a generic Radix Tree class
#'
#' @section Usage:
#' \preformatted{tree <- RadixTree$new()
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
#' tree$search(sequences, max_distance = NULL, max_fraction = NULL, mode = "levenshtein", nthreads = 1, show_progress = FALSE)
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
#'   \item{mode}{- In \code{$search()}, either Levenshtein or Hamming. Levenshtein will allows for insertions and deletions and calculates "edit distance". Hamming does not allow for insertions or deletions.}
#'   \item{nthreads}{- How many threads to use in the search.}
#'   \item{show_progress}{- Display progress.}
#' }
#'
#' @section Details:
#' \code{$new()} creates a new Tree object, which holds the pointer to the underlying C++ implentation. The RadixTree class accepts any strings of single-width characters. 
#'
#' \code{$print()} and \code{$to_string()} prints to screen or outputs the tree to a string representation.
#' 
#' \code{$to_vector()} outputs all sequences held by the tree. 
#' 
#' \code{$size()} outputs the size of the tree (i.e. how many sequences are contained). 
#' 
#' \code{$insert()}, \code{$erase()} and \code{$find()} insert, erase and find sequences in the tree, respectively.
#' 
#' \code{$search()} This function searches for similar sequences within a threshold (given by max_distance or max_fraction) based on Levenshtein or Hamming algorithms. 
#' The output of this function is a data.frame of all matches with columns "query" (the sequences input to the search function), 
#' "target" (the sequences inserted into the tree) and "distance" the absolute distance between query and target sequences. 
#' 
#' @seealso 
#' https://en.wikipedia.org/wiki/Radix_tree
#' 
#' @examples
#' # Plot Data: x,y,r
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

#' @name RadixTree
NULL


RadixTree <- R6::R6Class("RadixTree", list(
  xp = NULL,
  initialize = function(sequences = NULL) {
    self$xp <- RadixTree_create()
    if(!is.null(sequences)) {
      RadixTree_insert(self$xp, sequences)
    }
  },
  print = function() {
    cat(RadixTree_print(self$xp))
  },
  to_string = function() {
    RadixTree_print(self$xp)
  },
  graph = function(depth = -1, root_label = "root", plot = TRUE) {
    result <- RadixTree_graph(self$xp, depth)
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
    RadixTree_to_vector(self$xp)
  },
  size = function() {
    RadixTree_size(self$xp)
  },
  insert = function(sequences) {
    invisible(RadixTree_insert(self$xp, sequences))
  },
  erase = function(sequences) {
    invisible(RadixTree_erase(self$xp, sequences))
  },
  find = function(sequences) {
    RadixTree_find(self$xp, sequences)
  },
  prefix_search = function(sequences) {
    result <- RadixTree_prefix_search(self$xp, sequences)
    if(is.null(result)) {
      data.frame(query = character(0), target = character(0), stringsAsFactors=F)
    } else {
      result
    }
  },
  search = function(sequences, max_distance = NULL, max_fraction = NULL, mode = "levenshtein", nthreads = 1, show_progress = FALSE) {
    if(!is.null(max_distance)) {
      if(length(max_distance) == 1) {
        max_distance <- rep(max_distance, length(sequences))
      }
    } else if(!is.null(max_fraction)) {
      max_distance <- as.integer(nchar(sequences) * max_fraction)
    } else {
      stop("Either max_distance or max_fraction must be non-null")
    }
    if(mode == "levenshtein") {
      result <- RadixTree_levenshtein_search(self$xp, sequences, max_distance, nthreads, show_progress)
    } else if(mode == "hamming") {
      result <- RadixTree_hamming_search(self$xp, sequences, max_distance, nthreads, show_progress)
    } else if(mode == "anchored") {
      result <- RadixTree_anchored_search(self$xp, sequences, max_distance, nthreads, show_progress)
    } else {
      stop("mode should be levenshtein, hamming or anchored")
    }
    if(is.null(result)) {
      data.frame(query = character(0), target = character(0), distance = integer(0), stringsAsFactors=F)
    } else {
      result
    }
  },
  validate = function() {
    RadixTree_validate(self$xp)
  }
))

