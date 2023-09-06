#' @title split_search
#' @description Search for similar sequences based on splitting sequences into left and right sides
#' and searching for matches in each side using a bi-directional anchored alignment. 
#' @param query `r rdoc("query")`
#' @param target `r rdoc("target")`
#' @param query_split index to split query sequence. Should be within (edge_trim, nchar(query)-edge_trim] or -1 to indicate no split.
#' @param target_split index to split target sequence. Should be within (edge_trim, nchar(query)-edge_trim] or -1 to indicate no split.
#' @param edge_trim number of bases to trim from each side of the sequence (default value: 0).
#' @param max_distance `r rdoc("max_distance")`
#' @param ... additional arguments passed to `RadixTree$search`
#' @return data.frame with columns query, target, and distance.
#' @details
#' This function is useful for searching for similar sequences that may have variable windows of sequencing (e.g. different 5' and 3' primers)
#' but contain the same core sequence or position. 
#' The two split parameters partition the query and target sequences into left and right sides, 
#' where left = stri_sub(sequence, edge_trim+1, split) and right = stri_sub(query, split+1, -edge_trim-1).
#' @examples
#' # Consider two sets of sequences
#' # query1   AGACCTAA CCC
#' # target1 AAGACCTAA CC
#' # query2   GGGTGTAA CCACCC
#' # target2   GGTGTAA CCAC
#' # Despite having different frames, query1 and query2 and clearly match to target1 and target2, respectively.
#' # One could consider splitting based on a common core sequence, e.g. a common TAA stop codon. 
#' split_search(query=c(  "AGACCTAACCC", "GGGTGTAACCACCC"),
#'              target=c("AAGACCTAACC",   "GGTGTAACCAC"),
#'              query_split=c(8, 8),
#'              target_split=c(9, 7),
#'              edge_trim=0,
#'              max_distance=0)
split_search <- function(query, target, query_split, target_split, edge_trim = 0L, max_distance = 0L, ...) {
  # special case: split == -1 means left side is full sequence
  query_split <- ifelse(query_split < 0, nchar(query)-edge_trim, query_split)
  target_split <- ifelse(target_split < 0, nchar(target)-edge_trim, target_split)
  
  # split sequences
  query_left <- stringi::stri_reverse(stringi::stri_sub(query, edge_trim+1, query_split))
  query_right <- stringi::stri_sub(query, query_split+1, -edge_trim-1)
  target_left <- stringi::stri_reverse(stringi::stri_sub(target, edge_trim+1, target_split))
  target_right <- stringi::stri_sub(target, target_split+1, -edge_trim-1)
  
  # build trees
  left_tree <- RadixTree$new(unique(target_left))
  right_tree <- RadixTree$new(unique(target_right))
  
  # Search for similar sequences between lefts and rights
  left_matches <- left_tree$search(unique(query_left), max_distance = max_distance, mode = "anchored", ...)
  left_matches <- dplyr::rename(left_matches, query_left=query, target_left=target)
  right_matches <- right_tree$search(unique(query_right), max_distance = max_distance, mode = "anchored", ...)
  right_matches <- dplyr::rename(right_matches, query_right=query, target_right=target)
  
  # If either left or right finds no matches, return empty dataframe
  if(nrow(left_matches) == 0 || nrow(right_matches) == 0) {
    return(data.frame(query=character(0), target=character(0), distance=integer(0)))
  }
  
  # construct map of full sequence to left and right
  # filter in only potential matches, i.e. queries or targets that are in both left_matches and right_matches data.frame
  df_query <- data.frame(query, query_left, query_right)
  df_query <- dplyr::filter(df_query, query_left %in% left_matches$query_left, query_right %in% right_matches$query_right)
  df_target <- data.frame(target, target_left, target_right)
  df_target <- dplyr::filter(df_target, target_left %in% left_matches$target_left, target_right %in% right_matches$target_right)
  
  # Join results together, append full query and target sequences to left and right matches
  left_matches <- dplyr::inner_join(left_matches, df_query, by = "query_left")
  left_matches <- dplyr::inner_join(left_matches, df_target, by = "target_left")
  right_matches <- dplyr::inner_join(right_matches, df_query, by = "query_right")
  right_matches <- dplyr::inner_join(right_matches, df_target, by = "target_right")
  
  results <- dplyr::inner_join(left_matches, right_matches, by = c("query", "target"), suffix=c(".left", ".right"))
  results <- dplyr::mutate(results, distance = distance.left + distance.right)
  results <- dplyr::filter(results, distance <= max_distance)
  results <- dplyr::select(results, query, target, distance)
  as.data.frame(results)
}