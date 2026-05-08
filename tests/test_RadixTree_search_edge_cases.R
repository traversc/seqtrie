print("Running test_RadixTree_search_edge_cases.R")

if(requireNamespace("seqtrie", quietly=TRUE)) {
  library(seqtrie)

  sort_search_result <- function(result) {
    result <- as.data.frame(result, stringsAsFactors = FALSE)
    columns <- intersect(c("query", "target", "distance", "query_size", "target_size"), names(result))
    result <- result[, columns, drop = FALSE]
    if(nrow(result) > 0L) {
      result <- result[order(result$query, result$target, result$distance), , drop = FALSE]
    }
    rownames(result) <- NULL
    result
  }

  empty_result <- function(include_sizes = FALSE) {
    result <- data.frame(
      query = character(),
      target = character(),
      distance = integer(),
      stringsAsFactors = FALSE
    )
    if(include_sizes) {
      result$query_size <- integer()
      result$target_size <- integer()
    }
    result
  }

  check_identical <- function(label, actual, expected) {
    actual <- sort_search_result(actual)
    expected <- sort_search_result(expected)
    if(!identical(actual, expected)) {
      print(label)
      print(actual)
      print(expected)
    }
    stopifnot(identical(actual, expected))
  }

  expected_from_matrix <- function(query,
                                   target,
                                   max_distance,
                                   mode,
                                   cost_matrix = NULL,
                                   gap_cost = NA_integer_,
                                   gap_open_cost = NA_integer_) {
    distances <- seqtrie::dist_matrix(
      query = query,
      target = target,
      mode = mode,
      cost_matrix = cost_matrix,
      gap_cost = gap_cost,
      gap_open_cost = gap_open_cost,
      nthreads = 2L,
      show_progress = FALSE
    )
    include_sizes <- mode == "anchored"
    query_size <- attr(distances, "query_size")
    target_size <- attr(distances, "target_size")
    rows <- vector("list", length(query))
    for(i in seq_along(query)) {
      idx <- which(!is.na(distances[i, ]) & distances[i, ] <= max_distance[i])
      if(length(idx) == 0L) next
      row <- data.frame(
        query = query[i],
        target = target[idx],
        distance = as.integer(distances[i, idx]),
        stringsAsFactors = FALSE
      )
      if(include_sizes) {
        row$query_size <- as.integer(query_size[i, idx])
        row$target_size <- as.integer(target_size[i, idx])
      }
      rows[[i]] <- row
    }
    result <- do.call(rbind, rows)
    if(is.null(result)) {
      empty_result(include_sizes)
    } else {
      sort_search_result(result)
    }
  }

  check_tree_against_matrix <- function(label,
                                        query,
                                        target,
                                        max_distance,
                                        mode,
                                        cost_matrix = NULL,
                                        gap_cost = NA_integer_,
                                        gap_open_cost = NA_integer_) {
    tree <- RadixTree$new(target)
    actual <- tree$search(
      query = query,
      max_distance = max_distance,
      mode = mode,
      cost_matrix = cost_matrix,
      gap_cost = gap_cost,
      gap_open_cost = gap_open_cost,
      nthreads = 2L,
      show_progress = FALSE
    )
    expected <- expected_from_matrix(
      query = query,
      target = target,
      max_distance = max_distance,
      mode = mode,
      cost_matrix = cost_matrix,
      gap_cost = gap_cost,
      gap_open_cost = gap_open_cost
    )
    check_identical(label, actual, expected)
  }

  subtree_target <- c("", "A", "AA", "AAA", "AC", "CAA", "T")
  subtree <- RadixTree$new(subtree_target)

  expected_root_subtree <- data.frame(
    query = rep("", length(subtree_target)),
    target = subtree_target,
    distance = 0L,
    query_size = 0L,
    target_size = 0L,
    stringsAsFactors = FALSE
  )
  check_identical(
    "anchored empty query should add the whole subtree",
    subtree$search("", max_distance = 0L, mode = "anchored", nthreads = 2L, show_progress = FALSE),
    expected_root_subtree
  )

  expected_child_subtree <- data.frame(
    query = rep("A", 5L),
    target = c("", "A", "AA", "AAA", "AC"),
    distance = 0L,
    query_size = c(0L, 1L, 1L, 1L, 1L),
    target_size = c(0L, 1L, 1L, 1L, 1L),
    stringsAsFactors = FALSE
  )
  check_identical(
    "anchored prefix match should add terminal child descendants",
    subtree$search("A", max_distance = 0L, mode = "anchored", nthreads = 2L, show_progress = FALSE),
    expected_child_subtree
  )

  hamming_target <- c("AAAAC", "AAAAT", "AAAGT", "CCCC", "CCCT")
  hamming_query <- c("AAAAG", "CCCC", "CCCA", "AAAGT")
  expected_hamming <- data.frame(
    query = c("AAAAG", "AAAAG", "AAAGT", "CCCC", "CCCA", "CCCA"),
    target = c("AAAAC", "AAAAT", "AAAGT", "CCCC", "CCCC", "CCCT"),
    distance = c(1L, 1L, 0L, 0L, 1L, 1L),
    stringsAsFactors = FALSE
  )
  check_identical(
    "hamming should accumulate distance across compressed branches",
    RadixTree$new(hamming_target)$search(
      hamming_query,
      max_distance = c(1L, 0L, 1L, 0L),
      mode = "hamming",
      nthreads = 2L,
      show_progress = FALSE
    ),
    expected_hamming
  )

  target <- c("", "A", "AA", "AAA", "AC", "CAA", "T", "TTTT", "ACCC")
  query <- c("", "A", "AA", "AAC", "TT", "CCCC")
  max_distance <- c(0L, 0L, 1L, 2L, 2L, 3L)
  cost_matrix <- matrix(
    c(0L, 3L, 4L,
      2L, 0L, 3L,
      4L, 2L, 0L),
    nrow = 3L,
    byrow = TRUE
  )
  rownames(cost_matrix) <- colnames(cost_matrix) <- c("A", "C", "T")

  check_tree_against_matrix("global unit with empty strings and per-query thresholds", query, target, max_distance, "levenshtein")
  check_tree_against_matrix("anchored unit with empty strings and per-query thresholds", query, target, max_distance, "anchored")
  check_tree_against_matrix("global linear with empty strings and custom costs", query, target, max_distance + 2L, "levenshtein", cost_matrix = cost_matrix, gap_cost = 2L)
  check_tree_against_matrix("anchored linear with empty strings and custom costs", query, target, max_distance + 2L, "anchored", cost_matrix = cost_matrix, gap_cost = 2L)
  check_tree_against_matrix("global affine with empty strings and custom costs", query, target, max_distance + 3L, "levenshtein", cost_matrix = cost_matrix, gap_cost = 1L, gap_open_cost = 2L)
  check_tree_against_matrix("anchored affine with empty strings and custom costs", query, target, max_distance + 3L, "anchored", cost_matrix = cost_matrix, gap_cost = 1L, gap_open_cost = 2L)

  forest <- RadixForest$new(target)
  for(mode in c("hamming", "levenshtein")) {
    actual <- forest$search(query, max_distance = max_distance, mode = mode, nthreads = 2L, show_progress = FALSE)
    expected <- expected_from_matrix(query, target, max_distance, mode)
    check_identical(paste("RadixForest", mode, "edge cases"), actual, expected)
  }
}
