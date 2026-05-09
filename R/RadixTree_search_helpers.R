#' @title Distance search for similar sequences
#' @description Find similar sequences within a distance threshold
#' @param query `r rdoc("query")`
#' @param target `r rdoc("target")`
#' @param max_distance `r rdoc("max_distance")`
#' @param max_fraction `r rdoc("max_fraction")`
#' @param mode `r rdoc("mode")`
#' @param cost_matrix `r rdoc("cost_matrix")`
#' @param gap_cost `r rdoc("gap_cost")`
#' @param gap_open_cost `r rdoc("gap_open_cost")`
#' @param tree_class Which R6 class to use. Either RadixTree or RadixForest (default: RadixTree)
#' @param nthreads `r rdoc("nthreads")`
#' @param show_progress `r rdoc("show_progress")`
#' @details
#' This function finds all sequences in _target_ that are within a distance threshold of any sequence in _query_.
#' This function uses either a RadixTree or RadixForest to store _target_ sequences. See the R6 class documentation for additional details.
#' 
#' `r rdoc("details")`
#' @return The output is a data.frame of all matches with columns "query" and "target".
#' For anchored searches, the output also includes attributes "query_size" and "target_size"
#' which are vectors containing the portion of the query and target sequences that are aligned.
#' @examples
#' dist_search(c("ACGT", "AAAA"), c("ACG", "ACGT"), max_distance = 1, mode = "levenshtein")
#' @name dist_search
dist_search <- function(query, target, max_distance = NULL, max_fraction = NULL, mode = "levenshtein",
                        cost_matrix = NULL, gap_cost = NA_integer_, gap_open_cost = NA_integer_, tree_class = "RadixTree",
                        nthreads = 1, show_progress = FALSE) {
  if (!tree_class %in% c("RadixTree", "RadixForest")) {
    stop("tree_class must be one of RadixTree or RadixForest")
  }
  if (tree_class == "RadixTree") {
    obj <- RadixTree$new()
    obj$insert(target)
    obj$search(query, max_distance, max_fraction, mode, cost_matrix, gap_cost, gap_open_cost, nthreads, show_progress)
  } else if(tree_class == "RadixForest") {
    gap_cost_provided <- !(length(gap_cost) == 1L && is.na(gap_cost))
    gap_open_provided <- !(length(gap_open_cost) == 1L && is.na(gap_open_cost))
    if(!is.null(cost_matrix) || gap_cost_provided || gap_open_provided) {
      stop("cost_matrix, gap_cost and gap_open_cost are not supported for RadixForest")
    }
    obj <- RadixForest$new()
    obj$insert(target)
    obj$search(query, max_distance, max_fraction, mode, nthreads, show_progress)
  }
}

#' @title split_search
#' @description Search for similar sequences based on splitting sequences into left and right sides
#' and searching for matches in each side using a bi-directional anchored alignment. 
#' @param query `r rdoc("query")`
#' @param target `r rdoc("target")`
#' @param query_split index to split query sequence. Should be within (edge_trim, nchar(query)-edge_trim] or -1 to indicate no split.
#' @param target_split index to split target sequence. Should be within (edge_trim, nchar(target)-edge_trim] or -1 to indicate no split.
#' @param edge_trim number of bases to trim from each side of the sequence (default value: 0).
#' @param max_distance `r rdoc("max_distance")`
#' @param max_fraction `r rdoc("max_fraction")`
#' @param cost_matrix `r rdoc("cost_matrix")`
#' @param gap_cost `r rdoc("gap_cost")`
#' @param gap_open_cost `r rdoc("gap_open_cost")`
#' @param nthreads `r rdoc("nthreads")`
#' @param show_progress `r rdoc("show_progress")`
#' @return data.frame with columns query, target, and distance.
#' @details
#' This function is useful for searching for similar sequences that may have variable windows of sequencing (e.g. different 5' and 3' primers)
#' but contain the same core sequence or position. 
#' The two split parameters partition the query and target sequences into left and right sides, 
#' where left = rev(substr(sequence, edge_trim+1, split)) and right = substr(sequence, split+1, nchar(sequence)-edge_trim).
#' @examples
#' # Consider two sets of sequences
#' # query1   AGACCTAA CCC
#' # target1 AAGACCTAA CC
#' # query2   GGGTGTAA CCACCC
#' # target2   GGTGTAA CCAC
#' # Despite having different frames, query1 and query2 and clearly 
#' # match to target1 and target2, respectively.
#' # One could consider splitting based on a common core sequence, 
#' # e.g. a common TAA stop codon. 
#' split_search(query=c(  "AGACCTAACCC", "GGGTGTAACCACCC"),
#'              target=c("AAGACCTAACC",   "GGTGTAACCAC"),
#'              query_split=c(8, 8),
#'              target_split=c(9, 7),
#'              edge_trim=0,
#'              max_distance=0)
split_search <- function(query, target, query_split, target_split, edge_trim = 0L,
                         max_distance = 0L, max_fraction = NULL,
                         cost_matrix = NULL, gap_cost = NA_integer_, gap_open_cost = NA_integer_,
                         nthreads = 1, show_progress = FALSE) {
  if (!is.character(query) || !is.character(target)) {
    stop("query and target must be character vectors")
  }
  if (!is_integerlike(edge_trim) || length(edge_trim) != 1L || is.na(edge_trim) || edge_trim < 0L) {
    stop("edge_trim must be a single non-negative integer")
  }
  edge_trim <- as.integer(edge_trim)

  normalize_split <- function(split, sequences, label) {
    split <- recycle_arg(split, sequences)
    if (!is_integerlike(split) || any(is.na(split))) {
      stop(label, " must contain integer-like split positions")
    }
    split <- as.integer(split)
    if (any(split < 0L & split != -1L)) {
      stop(label, " values must be -1 or non-negative split positions")
    }
    seq_len <- nchar(sequences)
    split <- ifelse(split == -1L, seq_len - edge_trim, split)
    invalid <- split <= edge_trim | split > (seq_len - edge_trim)
    if (any(invalid)) {
      stop(label, " values must be -1 or within (edge_trim, nchar(sequence) - edge_trim]")
    }
    split
  }

  query_split <- normalize_split(query_split, query, "query_split")
  target_split <- normalize_split(target_split, target, "target_split")

  charset <- unique(c(get_charset(query), get_charset(target)))
  check_alignment_params("anchored", cost_matrix, gap_cost, gap_open_cost, charset, diag_must_be_zero = TRUE)

  if (!is.null(max_distance)) {
    max_distance <- recycle_arg(max_distance, query)
  } else if (!is.null(max_fraction)) {
    max_distance <- as.integer(nchar(query) * recycle_arg(max_fraction, query))
  } else {
    stop("Either max_distance or max_fraction must be non-null")
  }
  if (!is_integerlike(max_distance) || any(is.na(max_distance)) || any(max_distance < 0)) {
    stop("max_distance/max_fraction must be non-negative integer-like values")
  }
  max_distance <- as.integer(max_distance)

  if (!is.na(gap_open_cost) && !is.na(gap_cost) && gap_open_cost > 0L) {
    gap_open_cost <- gap_open_cost + gap_cost
  }

  c_split_search(query, target, query_split, target_split, edge_trim, max_distance,
                 cost_matrix, gap_cost, gap_open_cost, nthreads, show_progress)
}
