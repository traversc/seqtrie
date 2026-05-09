print("Running test_split_search.R")

if(requireNamespace("seqtrie", quietly = TRUE)) {
  library(seqtrie)

  sort_result <- function(x) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    x <- x[c("query", "target", "distance")]
    if(nrow(x) > 0L) {
      x <- x[order(x$query, x$target, x$distance), , drop = FALSE]
    }
    rownames(x) <- NULL
    x
  }

  empty_result <- function() {
    data.frame(query = character(), target = character(), distance = integer(), stringsAsFactors = FALSE)
  }

  reverse_strings <- function(x) {
    vapply(strsplit(x, "", fixed = TRUE), function(chars) paste0(rev(chars), collapse = ""), character(1))
  }

  split_parts <- function(sequences, split, edge_trim) {
    split <- ifelse(split == -1L, nchar(sequences) - edge_trim, split)
    data.frame(
      left = reverse_strings(substr(sequences, edge_trim + 1L, split)),
      right = substr(sequences, split + 1L, nchar(sequences) - edge_trim),
      stringsAsFactors = FALSE
    )
  }

  split_search_expected <- function(query,
                                    target,
                                    query_split,
                                    target_split,
                                    edge_trim = 0L,
                                    max_distance = 0L,
                                    max_fraction = NULL,
                                    cost_matrix = NULL,
                                    gap_cost = NA_integer_,
                                    gap_open_cost = NA_integer_) {
    if(is.null(max_distance)) {
      max_distance <- as.integer(nchar(query) * rep(max_fraction, length.out = length(query)))
    } else {
      max_distance <- rep(max_distance, length.out = length(query))
    }

    keep_target <- !duplicated(target)
    target <- target[keep_target]
    target_split <- target_split[keep_target]

    query_parts <- split_parts(query, query_split, edge_trim)
    target_parts <- split_parts(target, target_split, edge_trim)

    left_dist <- seqtrie::dist_matrix(
      query_parts$left,
      target_parts$left,
      mode = "anchored",
      cost_matrix = cost_matrix,
      gap_cost = gap_cost,
      gap_open_cost = gap_open_cost,
      nthreads = 2L
    )
    right_dist <- seqtrie::dist_matrix(
      query_parts$right,
      target_parts$right,
      mode = "anchored",
      cost_matrix = cost_matrix,
      gap_cost = gap_cost,
      gap_open_cost = gap_open_cost,
      nthreads = 2L
    )
    distances <- left_dist + right_dist

    rows <- vector("list", length(query))
    for(i in seq_along(query)) {
      idx <- which(!is.na(distances[i, ]) & distances[i, ] <= max_distance[i])
      if(length(idx) == 0L) next
      rows[[i]] <- data.frame(
        query = query[i],
        target = target[idx],
        distance = as.integer(distances[i, idx]),
        stringsAsFactors = FALSE
      )
    }

    out <- do.call(rbind, rows)
    if(is.null(out)) empty_result() else sort_result(out)
  }

  check_split_search <- function(label, query, target, query_split, target_split, ...) {
    actual <- sort_result(split_search(query, target, query_split, target_split, nthreads = 2L, ...))
    expected <- split_search_expected(query, target, query_split, target_split, ...)
    if(!identical(actual, expected)) {
      print(label)
      print(actual)
      print(expected)
    }
    stopifnot(identical(actual, expected))
  }

  random_strings <- function(n, charset = "ACGT", min_len = 4L, max_len = 9L) {
    chars <- strsplit(charset, "", fixed = TRUE)[[1]]
    len <- sample(min_len:max_len, n, replace = TRUE)
    vapply(len, function(x) paste0(sample(chars, x, replace = TRUE), collapse = ""), character(1))
  }

  random_split <- function(sequences, edge_trim) {
    vapply(sequences, function(x) {
      if(runif(1) < 0.2) return(-1L)
      sample(seq.int(edge_trim + 1L, nchar(x) - edge_trim), 1L)
    }, integer(1))
  }

  check_split_search(
    "example with duplicate target collapsed",
    query = c("AGACCTAACCC", "GGGTGTAACCACCC"),
    target = c("AAGACCTAACC", "GGTGTAACCAC", "AAGACCTAACC"),
    query_split = c(8L, 8L),
    target_split = c(9L, 7L, 9L),
    edge_trim = 0L,
    max_distance = 0L
  )

  check_split_search(
    "no matches returns empty result",
    query = "AAAAAA",
    target = c("CCCCCC", "GGGGGG"),
    query_split = 3L,
    target_split = c(3L, 3L),
    max_distance = 0L
  )

  set.seed(12)
  for(iter in seq_len(20L)) {
    edge_trim <- sample(0:1, 1L)
    query <- random_strings(8L, min_len = 2L * edge_trim + 3L)
    query <- c(query, query[1L])
    target <- random_strings(12L, min_len = 2L * edge_trim + 3L)
    target <- c(target, target[1L], target[3L])
    query_split <- random_split(query, edge_trim)
    target_split <- random_split(target, edge_trim)
    max_distance <- sample(0:4, length(query), replace = TRUE)

    check_split_search(
      paste("unit fuzz", iter),
      query = query,
      target = target,
      query_split = query_split,
      target_split = target_split,
      edge_trim = edge_trim,
      max_distance = max_distance
    )
  }

  set.seed(23)
  query <- random_strings(10L, min_len = 5L)
  target <- random_strings(14L, min_len = 5L)
  query_split <- random_split(query, 1L)
  target_split <- random_split(target, 1L)
  check_split_search(
    "max_fraction fuzz",
    query = query,
    target = target,
    query_split = query_split,
    target_split = target_split,
    edge_trim = 1L,
    max_distance = NULL,
    max_fraction = 0.35
  )

  cost_matrix <- seqtrie::generate_cost_matrix("ACGT", match = 0L, mismatch = 2L)
  check_split_search(
    "linear custom cost fuzz",
    query = query,
    target = target,
    query_split = query_split,
    target_split = target_split,
    edge_trim = 1L,
    max_distance = 5L,
    cost_matrix = cost_matrix,
    gap_cost = 2L
  )
  check_split_search(
    "affine custom cost fuzz",
    query = query,
    target = target,
    query_split = query_split,
    target_split = target_split,
    edge_trim = 1L,
    max_distance = 6L,
    cost_matrix = cost_matrix,
    gap_cost = 1L,
    gap_open_cost = 2L
  )

  stopifnot(inherits(
    try(split_search("ACGT", "ACGT", query_split = 0L, target_split = 2L), silent = TRUE),
    "try-error"
  ))
  stopifnot(inherits(
    try(split_search("ACGT", "ACGT", query_split = -2L, target_split = 2L), silent = TRUE),
    "try-error"
  ))
  stopifnot(inherits(
    try(split_search("ACGT", "ACGT", query_split = 2L, target_split = 4L, edge_trim = 1L), silent = TRUE),
    "try-error"
  ))
}
