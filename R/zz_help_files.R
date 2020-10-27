#' td_suffix_tree
#' 
#' Creates a suffix tree
#' @usage td_suffix_tree(subject, min_length)
#' @param subject A character vector
#' @param min_length Minimum suffix length in the tree
#' @return A suffix tree: a list of class "suffix_tree" containing an external pointer
#' @details 
#' A suffix tree can be used to perform partial hamming matching. 
#' @examples 
#' x <- td_suffix_tree("hello")
#' @name td_suffix_tree
NULL

#' td_prefix_tree
#' 
#' Creates a prefix tree
#' @usage td_prefix_tree(subject)
#' @param subject A character vector
#' @return A prefix tree: a list of class "prefix_tree" containing an external pointer
#' @details 
#' A prefix tree can be used to perform levenshtein and hamming calculations. 
#' @examples 
#' x <- td_prefix_tree("hello")
#' @name td_prefix_tree
NULL

#' td_partial_hamming
#' 
#' Performs partial hamming matching using a suffix tree
#' @usage td_partial_hamming(query, subject=NULL, max_distance = NA_integer_, symmetric = FALSE, nthreads = 1L)
#' @param query A character vector of strings as query
#' @param subject A character vector of strings as subject, a suffix_tree or NULL. See details
#' @param max_distance The max hamming distance to search for matches.
#' @param symmetric Whether to treat query and subject as the same (i.e. pairwise matching)
#' @param nthreads Number of threads to use in computation
#' @return A data.frame containing query, subject and distance.
#' @details 
#' Subject parameter can be a character vector, a suffix_tree or NULL. If NULL, the function will perform a pairwise search on query. 
#' @examples 
#' x <- td_partial_hamming(query = "ell", subject = "hello")
#' @name td_partial_hamming
NULL

#' td_hamming
#' 
#' Performs hamming matching using a prefix tree
#' @usage td_hamming(query, subject=NULL, max_distance = NA_integer_, symmetric = FALSE, nthreads = 1L)
#' @param query A character vector of strings as query
#' @param subject A character vector of strings as subject, a prefix_tree or NULL. See details
#' @param max_distance The max hamming distance to search for matches.
#' @param symmetric Whether to treat query and subject as the same (i.e. pairwise matching)
#' @param nthreads Number of threads to use in computation
#' @return A data.frame containing query, subject and distance.
#' @details 
#' Subject parameter can be a character vector, a prefix_tree or NULL. If NULL, the function will perform a pairwise search on query. 
#' @examples 
#' x <- td_hamming(query = "heloo", subject = "hello")
#' @name td_hamming
NULL

#' td_levenshtein
#' 
#' Performs levenshtein matching using a prefix tree
#' @usage td_levenshtein(query, subject=NULL, max_distance = NA_integer_, symmetric = FALSE, nthreads = 1L)
#' @param query A character vector of strings as query
#' @param subject A character vector of strings as subject, a prefix_tree or NULL. See details
#' @param max_distance The max levenshtein distance to search for matches.
#' @param symmetric Whether to treat query and subject as the same (i.e. pairwise matching)
#' @param nthreads Number of threads to use in computation
#' @return A data.frame containing query, subject and distance.
#' @details 
#' Subject parameter can be a character vector, a prefix_tree or NULL. If NULL, the function will perform a pairwise search on query. 
#' @examples 
#' x <- td_levenshtein(query = "helo", subject = "hello")
#' @name td_levenshtein
NULL


#' td_partial_levenshtein
#' 
#' Performs partial levenshtein matching using a suffix tree
#' @usage td_partial_levenshtein(query, subject=NULL, anchor = "left", max_distance = NA_integer_, symmetric = FALSE, nthreads = 1L)
#' @param query A character vector of strings as query
#' @param subject A character vector of strings as subject, a suffix_tree or NULL. See details
#' @param anchor "left" or "right". Determines whether the left or right side of an alignment's end gap should be penalized
#' @param max_distance The max levenshtein distance to search for matches.
#' @param symmetric Whether to treat query and subject as the same (i.e. pairwise matching)
#' @param nthreads Number of threads to use in computation
#' @return A data.frame containing query, subject and distance.
#' @details 
#' Subject parameter can be a character vector, a suffix_tree or NULL. If NULL, the function will perform a pairwise search on query. 
#' @examples 
#' x <- td_partial_levenshtein(query = "hell", subject = "hello", anchor = "left")
#' @name td_partial_levenshtein
NULL