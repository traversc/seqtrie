#' RadixTree
#'
#' Radix Tree (trie) class implementation
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
#'   \item{mode}{- In \code{$search()}, One of hamming (hm), levenshtein (lv) or anchored (an). Levenshtein will allows for insertions and deletions and calculates "edit distance". Hamming does not allow for insertions or deletions.}
#'   \item{cost_matrix}{In \code{$search()}, a cost matrix for use with Levenshtein or Anchored searches. See details.}
#'   \item{gap_cost}{In \code{$search()}, a cost matrix for use with Levenshtein or Anchored searches. See details.}
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
#' \code{$search()} This function searches for similar sequences within a threshold (given by max_distance or max_fraction) based on Hamming, Levenshtein or Anchored algorithms.
#' An anchored search is a form of semi-global alignment, where the query sequence is "anchored" (global) to the beginning of both the query and target sequences, 
#' but is semi-global in that the end of the either the query sequence or target sequence (but not both) can be unaligned. This contrasts with the Levenshtein distance which is global at 
#' both the starts and ends of the sequences. 
#' 
#' The output of this function is a data.frame of all matches with columns "query" (the sequences input to the search function), 
#' "target" (the sequences inserted into the tree) and "distance" the absolute distance between query and target sequences. 
#' For anchored searches, the output also includes "query_size" and "target_size" which are the partial lengths of the query and target sequences that are aligned. 
#' 
#' @seealso 
#' https://en.wikipedia.org/wiki/Radix_tree
#' 
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
#' @name RadixTree
NULL


RadixTree <- R6::R6Class("RadixTree", list(
  root_pointer = NULL,
  char_counter_pointer = NULL,
  initialize = function(sequences = NULL) {
    self$root_pointer <- RadixTree_create()
    self$char_counter_pointer <- CharCounter_create()
    if(!is.null(sequences)) {
      RadixTree_insert(self$root_pointer, sequences)
    }
  },
  print = function() {
    cat(RadixTree_print(self$root_pointer))
  },
  to_string = function() {
    RadixTree_print(self$root_pointer)
  },
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
  to_vector = function() {
    RadixTree_to_vector(self$root_pointer)
  },
  size = function() {
    RadixTree_size(self$root_pointer)
  },
  insert = function(sequences) {
    result <- RadixTree_insert(self$root_pointer, sequences)
    CharCounter_add(self$char_counter_pointer, sequences[result])
    invisible(result)
  },
  erase = function(sequences) {
    result <- RadixTree_erase(self$root_pointer, sequences)
    CharCounter_subtract(self$char_counter_pointer, sequences[result])
    invisible(result)
  },
  find = function(sequences) {
    RadixTree_find(self$root_pointer, sequences)
  },
  prefix_search = function(sequences) {
    result <- RadixTree_prefix_search(self$root_pointer, sequences)
    if(is.null(result)) {
      data.frame(query = character(0), target = character(0), stringsAsFactors=F)
    } else {
      result
    }
  },
  search = function(sequences, max_distance = NULL, max_fraction = NULL, mode = "levenshtein", cost_matrix = NULL, gap_cost = NULL, nthreads = 1, show_progress = FALSE) {
    charset <- unique(c(CharCounter_keys(self$char_counter_pointer), get_charset(sequences)))
    check_alignment_params(mode, cost_matrix, gap_cost, charset)
    cost_matrix <- append_gap_cost(cost_matrix, gap_cost)
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
      RadixTree_hamming_search(self$root_pointer, sequences, max_distance, cost_matrix, nthreads, show_progress)
    } else if(mode == "levenshtein") {
      RadixTree_levenshtein_search(self$root_pointer, sequences, max_distance, cost_matrix, nthreads, show_progress)
    } else if(mode == "anchored") { # Append query_start and target_start columns from dist_pairwise anchored search
      result <- RadixTree_anchored_search(self$root_pointer, sequences, max_distance, cost_matrix, nthreads, show_progress)
      result2 <- dist_pairwise(result$query, result$target, mode = "anchored", cost_matrix = cost_matrix, nthreads = nthreads, show_progress = FALSE)
      if(any(result$distance != result2)) {
        stop("Internal error: anchored search results do not match pairwise results")
      }
      result$query_size <- attr(result2, "query_size")
      result$target_size <- attr(result2, "target_size")
      result
    }
  },
  validate = function() {
    RadixTree_validate(self$root_pointer)
  }
))

