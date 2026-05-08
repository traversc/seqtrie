print("Running test_RadixTree_affine_fuzz.R")

if(requireNamespace("seqtrie", quietly=TRUE)) {
  library(seqtrie)

  enumerate_strings <- function(charset, max_len) {
    out <- ""
    for(len in seq_len(max_len)) {
      grid <- expand.grid(rep(list(charset), len), stringsAsFactors = FALSE)
      out <- c(out, apply(grid, 1L, paste0, collapse = ""))
    }
    unique(out)
  }

  sort_search_result <- function(result) {
    result <- as.data.frame(result, stringsAsFactors = FALSE)
    columns <- intersect(c("query", "target", "distance", "query_size", "target_size"), names(result))
    result <- result[, columns, drop = FALSE]
    if(nrow(result) > 0L) {
      result <- result[do.call(order, result[, columns, drop = FALSE]), , drop = FALSE]
    }
    rownames(result) <- NULL
    result
  }

  empty_result <- function(include_sizes) {
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

  expected_from_distances <- function(distances, query, target, max_distance, include_sizes) {
    if(length(max_distance) == 1L) {
      max_distance <- rep(max_distance, length(query))
    }
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
    if(is.null(result)) empty_result(include_sizes) else sort_search_result(result)
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

  high_mismatch_cost <- function(charset, mismatch) {
    cost_matrix <- matrix(
      as.integer(mismatch),
      nrow = length(charset),
      ncol = length(charset),
      dimnames = list(charset, charset)
    )
    diag(cost_matrix) <- 0L
    cost_matrix
  }

  asymmetric_cost <- function(charset) {
    n <- length(charset)
    values <- ((seq_len(n * n) * 3L) %% 7L) + 1L
    cost_matrix <- matrix(as.integer(values), nrow = n, ncol = n)
    diag(cost_matrix) <- 0L
    dimnames(cost_matrix) <- list(charset, charset)
    cost_matrix
  }

  run_affine_case <- function(label,
                              query,
                              target,
                              cost_matrix,
                              gap_cost,
                              gap_open_cost,
                              max_distance_sets) {
    query <- unique(query)
    target <- unique(target)
    tree <- RadixTree$new(target)
    stopifnot(tree$validate())

    for(mode in c("levenshtein", "anchored")) {
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

      for(set_name in names(max_distance_sets)) {
        max_distance <- max_distance_sets[[set_name]]
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
        expected <- expected_from_distances(distances, query, target, max_distance, include_sizes)
        check_identical(
          paste(label, mode, set_name, sep = " / "),
          actual,
          expected
        )
      }
    }
  }

  small_charset <- c("A", "C")
  small_strings <- enumerate_strings(small_charset, 5L)
  small_thresholds <- list(
    scalar_tight = 2L,
    scalar_wide = 5L,
    per_query = as.integer((seq_along(small_strings) - 1L) %% 6L)
  )

  run_affine_case(
    label = "small alphabet exhaustive high mismatch",
    query = small_strings,
    target = small_strings,
    cost_matrix = high_mismatch_cost(small_charset, 5L),
    gap_cost = 1L,
    gap_open_cost = 3L,
    max_distance_sets = small_thresholds
  )

  run_affine_case(
    label = "small alphabet exhaustive asymmetric costs",
    query = small_strings,
    target = small_strings,
    cost_matrix = asymmetric_cost(small_charset),
    gap_cost = 2L,
    gap_open_cost = 1L,
    max_distance_sets = small_thresholds
  )

  boundary_charset <- c("A", "C", "G", "T")
  boundary_target <- c(
    "",
    "A", "AA", "AAA", "AAAA", "AAAAAA", "AAAAAAAA",
    "CAAAAAAAA", "AAAAAAAAC", "AAAACAAAA",
    "C", "CC", "CCCC", "CCCCCCCC", "CCCCAAAACCCC",
    "ACACACAC", "ACACACACAC",
    "TAAAAAAAAT", "GCCCCCCCG", "AAAAGGGGCCCC"
  )
  boundary_query <- c(
    "",
    "A", "AAAA", "AAAAAAA", "AAAAAAAAAA",
    "CAAAAAAAAC", "AAAAACAAAA", "AAAACAAA",
    "C", "CCCCC", "CCCCCCCCCC", "CCCCAACCCC",
    "ACACAC", "ACACACACACAC",
    "TAAAAT", "TAAAAAAAAAAT", "GCCCCG", "AAAAGGCCCC"
  )
  boundary_thresholds <- list(
    max_1 = 1L,
    max_3 = 3L,
    max_6 = 6L,
    per_query_boundary = as.integer(pmin(8L, 1L + (nchar(boundary_query) %% 8L)))
  )

  run_affine_case(
    label = "constructed affine band boundary high mismatch",
    query = boundary_query,
    target = boundary_target,
    cost_matrix = high_mismatch_cost(boundary_charset, 6L),
    gap_cost = 1L,
    gap_open_cost = 2L,
    max_distance_sets = boundary_thresholds
  )

  run_affine_case(
    label = "constructed affine band boundary asymmetric costs",
    query = boundary_query,
    target = boundary_target,
    cost_matrix = asymmetric_cost(boundary_charset),
    gap_cost = 2L,
    gap_open_cost = 2L,
    max_distance_sets = boundary_thresholds
  )
}
