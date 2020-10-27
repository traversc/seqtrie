td_suffix_tree <- function(subject, min_length = 1) {
  tree <- c_td_suffix_tree(subject, min_length)
  structure(list(tree = tree), class = "suffix_tree")
}

td_prefix_tree <- function(subject) {
  tree <- c_td_prefix_tree(subject)
  structure(list(tree = tree), class = "prefix_tree")
}

td_partial_hamming <- function(query, subject=NULL, max_distance = .Machine$integer.max, symmetric = FALSE, nthreads = 1L) {
  if(is.null(subject)) {
    tree <- c_td_prefix_tree(query)
    c_td_partial_hamming(tree, query, max_distance, symmetric, nthreads)
  } else if(is.character(subject)) {
    tree <- c_td_prefix_tree(subject)
    c_td_partial_hamming(tree, query, max_distance, symmetric, nthreads)
  } else if("prefix_tree" %in% class(subject)) {
    c_td_partial_hamming(subject$tree, query, max_distance, symmetric, nthreads)
  } else {
    stop("improper subject parameter")
  }
}

td_levenshtein <- function(query, subject=NULL, max_distance = .Machine$integer.max, symmetric = FALSE, nthreads = 1L) {
  if(is.null(subject)) {
    tree <- c_td_prefix_tree(query)
    c_td_levenshtein(tree, query, max_distance, symmetric, nthreads)
  } else if(is.character(subject)) {
    tree <- c_td_prefix_tree(subject)
    c_td_levenshtein(tree, query, max_distance, symmetric, nthreads)
  } else if("prefix_tree" %in% class(subject)) {
    c_td_levenshtein(subject$tree, query, max_distance, symmetric, nthreads)
  } else {
    stop("improper subject parameter")
  }
}

td_hamming <- function(query, subject=NULL, max_distance = .Machine$integer.max, symmetric = FALSE, nthreads = 1L) {
  if(is.null(subject)) {
    tree <- c_prefix_tree(query)
    c_td_hamming(tree, query, max_distance, symmetric, nthreads)
  } else if(is.character(subject)) {
    tree <- c_td_prefix_tree(subject)
    c_td_hamming(tree, query, max_distance, symmetric, nthreads)
  } else if("prefix_tree" %in% class(subject)) {
    c_td_hamming(subject$tree, query, max_distance, symmetric, nthreads)
  } else {
    stop("improper subject parameter")
  }
}

td_partial_levenshtein <- function(query, subject=NULL, anchor = "left", max_distance = .Machine$integer.max, symmetric = FALSE, nthreads = 1L) {
  if(is.null(subject)) {
    tree <- c_td_prefix_tree(query)
    c_td_partial_levenshtein(tree, query, anchor, max_distance, symmetric, nthreads)
  } else if(is.character(subject)) {
    tree <- c_td_prefix_tree(subject)
    c_td_partial_levenshtein(tree, query, anchor, max_distance, symmetric, nthreads)
  } else if("prefix_tree" %in% class(subject)) {
    c_td_partial_levenshtein(subject$tree, query, anchor, max_distance, symmetric, nthreads)
  } else {
    stop("improper subject parameter")
  }
}

