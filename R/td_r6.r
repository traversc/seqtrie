# export("DNATree_create", "RadixTree_create", "PrefixTree_create", 
#        "DNATree_size", "RadixTree_size", "PrefixTree_size",
#        "DNATree_print", "RadixTree_print", "PrefixTree_print",
#        "DNATree_insert", "RadixTree_insert", "PrefixTree_insert",
#        "DNATree_erase", "RadixTree_erase", "PrefixTree_erase",
#        "DNATree_find", "RadixTree_find", "PrefixTree_find",
#        "DNATree_levenshtein", "RadixTree_levenshtein", "PrefixTree_levenshtein",
#        "DNATree_hamming", "RadixTree_hamming", "PrefixTree_hamming")


# DNATree <- R6::R6Class("DNATree", list(
#   xp = NULL,
#   initialize = function(sequences = NULL) {
#     self$xp <- DNATree_create()
#     if(!is.null(sequences)) {
#       DNATree_insert(self$xp, sequences)
#     }
#   },
#   print = function() {
#     cat(DNATree_print(self$xp))
#   },
#   to_string = function() {
#     DNATree_print(self$xp)
#   },
#   to_dataframe = function() {
#     DNATree_to_dataframe(self$xp)
#   },
#   size = function() {
#     DNATree_size(self$xp)
#   },
#   insert = function(sequences) {
#     invisible(DNATree_insert(self$xp, sequences))
#   },
#   erase = function(sequences) {
#     invisible(DNATree_erase(self$xp, sequences))
#   },
#   find = function(sequences) {
#     !is.na(DNATree_find(self$xp, sequences))
#   },
#   find_prefix = function(sequences) {
#     result <- DNATree_find_prefix(self$xp, sequences)
#     if(is.null(result)) {
#       data.frame(query = character(0), target = character(0), stringsAsFactors=F)
#     } else {
#       result
#     }
#   },
#   search = function(sequences, max_distance = NULL, max_fraction = NULL, mode = "levenshtein", nthreads = 1, show_progress = FALSE) {
#     if(!is.null(max_distance)) {
#       if(length(max_distance) == 1) {
#         max_distance <- rep(max_distance, length(sequences))
#       }
#     } else if(!is.null(max_fraction)) {
#       max_distance <- as.integer(nchar(sequences) * max_fraction)
#     } else {
#       stop("Either max_distance or max_fraction must be non-null")
#     }
#     if(mode == "levenshtein") {
#       result <- DNATree_levenshtein(self$xp, sequences, max_distance, nthreads, show_progress)
#     } else if(mode == "hamming") {
#       result <- DNATree_hamming(self$xp, sequences, max_distance, nthreads, show_progress)
#     } else {
#       stop("mode should be levenshtein or hamming")
#     }
#     if(is.null(result)) {
#       data.frame(query = character(0), target = character(0), distance = integer(0), stringsAsFactors=F)
#     } else {
#       result
#     }
#   }
# ))


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
  graph = function(depth = -1, root_label = "(root)", plot = TRUE) {
    result <- RadixTree_graph(self$xp, depth)
    if(is.null(result)) {
      result <- data.frame(parent = character(0), child = character(0), stringsAsFactors=F)
    } else if(plot) {
      result$parent <- ifelse(result$parent == "", root_label, result$parent)
      gr <- igraph::graph_from_data_frame(result)
      igraph::V(gr)$color <- ifelse(names(igraph::V(gr)) == root_label, "white", "skyblue2")
      igraph::V(gr)$size <- ifelse(names(igraph::V(gr)) == root_label, 21, 15)
      plot(gr, layout=igraph::layout.fruchterman.reingold, vertex.color=igraph::V(gr)$color, vertex.label.family = "sans") 
    }
    invisible(result)
  },
  to_dataframe = function() {
    RadixTree_to_dataframe(self$xp)
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
    !is.na(RadixTree_find(self$xp, sequences))
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
    } else {
      stop("mode should be levenshtein or hamming")
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




# PrefixTree <- R6::R6Class("PrefixTree", list(
#   xp = NULL,
#   initialize = function(sequences = NULL) {
#     self$xp <- PrefixTree_create()
#     if(!is.null(sequences)) {
#       PrefixTree_insert(self$xp, sequences)
#     }
#   },
#   print = function() {
#     cat(PrefixTree_print(self$xp))
#   },
#   to_string = function() {
#     PrefixTree_print(self$xp)
#   },
#   to_dataframe = function() {
#     PrefixTree_to_dataframe(self$xp)
#   },
#   size = function() {
#     PrefixTree_size(self$xp)
#   },
#   insert = function(sequences) {
#     invisible(PrefixTree_insert(self$xp, sequences))
#   },
#   erase = function(sequences) {
#     invisible(PrefixTree_erase(self$xp, sequences))
#   },
#   find = function(sequences) {
#     !is.na(PrefixTree_find(self$xp, sequences))
#   },
#   find_prefix = function(sequences) {
#     result <- PrefixTree_find_prefix(self$xp, sequences)
#     if(is.null(result)) {
#       data.frame(query = character(0), target = character(0), stringsAsFactors=F)
#     } else {
#       result
#     }
#   },
#   search = function(sequences, max_distance = NULL, max_fraction = NULL, mode = "levenshtein", nthreads = 1, show_progress = FALSE) {
#     if(!is.null(max_distance)) {
#       if(length(max_distance) == 1) {
#         max_distance <- rep(max_distance, length(sequences))
#       }
#     } else if(!is.null(max_fraction)) {
#       max_distance <- as.integer(nchar(sequences) * max_fraction)
#     } else {
#       stop("Either max_distance or max_fraction must be non-null")
#     }
#     if(mode == "levenshtein") {
#       result <- PrefixTree_levenshtein(self$xp, sequences, max_distance, nthreads, show_progress)
#     } else if(mode == "hamming") {
#       result <- PrefixTree_hamming(self$xp, sequences, max_distance, nthreads, show_progress)
#     } else {
#       stop("mode should be levenshtein or hamming")
#     }
#     if(is.null(result)) {
#       data.frame(query = character(0), target = character(0), distance = integer(0), stringsAsFactors=F)
#     } else {
#       result
#     }
#   }
# ))


