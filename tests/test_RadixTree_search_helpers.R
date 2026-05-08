print("Running test_RadixTree_search_helpers.R")

if (requireNamespace("seqtrie", quietly = TRUE) &&
    requireNamespace("dplyr", quietly = TRUE)) {
  library(seqtrie)
  library(dplyr)

  query <- c("ACGT", "ACGA", "TTTA")
  target <- c("ACGT", "ACGG", "ACGA", "TTTT", "GGGG")

  tree <- RadixTree$new(target)
  direct_tree <- tree$search(query, max_distance = 1L, mode = "levenshtein", nthreads = 2L) %>%
    arrange(query, target)
  helper_tree <- dist_search(query, target, max_distance = 1L, mode = "levenshtein", tree_class = "RadixTree", nthreads = 2L) %>%
    arrange(query, target)
  stopifnot(identical(helper_tree, direct_tree))

  forest <- RadixForest$new(target)
  direct_forest <- forest$search(query, max_distance = 1L, mode = "levenshtein", nthreads = 2L) %>%
    arrange(query, target)
  helper_forest <- dist_search(query, target, max_distance = 1L, mode = "levenshtein", tree_class = "RadixForest", nthreads = 2L) %>%
    arrange(query, target)
  stopifnot(identical(helper_forest, direct_forest))

  stopifnot(inherits(
    try(dist_search(query, target, max_distance = 1L, tree_class = "RadixForest", gap_cost = 1L), silent = TRUE),
    "try-error"
  ))

  split_result <- split_search(
    query = c("AGACCTAACCC", "GGGTGTAACCACCC"),
    target = c("AAGACCTAACC", "GGTGTAACCAC"),
    query_split = c(8L, 8L),
    target_split = c(9L, 7L),
    edge_trim = 0L,
    max_distance = 0L,
    nthreads = 2L
  ) %>%
    arrange(query, target)

  expected_split <- data.frame(
    query = c("AGACCTAACCC", "GGGTGTAACCACCC"),
    target = c("AAGACCTAACC", "GGTGTAACCAC"),
    distance = c(0L, 0L),
    stringsAsFactors = FALSE
  )
  stopifnot(identical(split_result, expected_split))
}
