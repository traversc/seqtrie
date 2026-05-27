print("Running test_utilities.R")
runtime <- Sys.time()

if(requireNamespace("seqtrie", quietly = TRUE)) {
  library(seqtrie)

  NTHREADS <- 2L

  test_seed <- Sys.getenv("SEQTRIE_TEST_SEED")
  if(nzchar(test_seed)) {
    test_seed <- as.integer(test_seed)
  } else {
    test_seed <- as.integer(as.numeric(Sys.time())) %% .Machine$integer.max
  }
  cat("Test seed:", test_seed, "\n")
  set.seed(test_seed)

  sort_result <- function(x) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    columns <- intersect(c("query", "target", "distance", "query_size", "target_size"), names(x))
    x <- x[, columns, drop = FALSE]
    if(nrow(x) > 0L) {
      x <- unique(x)
      x <- x[do.call(order, x[, columns, drop = FALSE]), , drop = FALSE]
    }
    rownames(x) <- NULL
    x
  }

  empty_search_result <- function(include_sizes = FALSE) {
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

  expected_search <- function(query,
                              target,
                              max_distance,
                              mode,
                              cost_matrix = NULL,
                              gap_cost = NA_integer_,
                              gap_open_cost = NA_integer_) {
    distances <- dist_matrix(
      query,
      target,
      mode = mode,
      cost_matrix = cost_matrix,
      gap_cost = gap_cost,
      gap_open_cost = gap_open_cost,
      nthreads = NTHREADS,
      show_progress = FALSE
    )
    max_distance <- rep(max_distance, length.out = length(query))
    include_sizes <- mode %in% c("anchored", "an", "extension", "en")
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
    if(is.null(result)) empty_search_result(include_sizes) else sort_result(result)
  }

  expected_self_search <- function(query,
                                   max_distance,
                                   mode,
                                   cost_matrix = NULL,
                                   gap_cost = NA_integer_,
                                   gap_open_cost = NA_integer_) {
    distances <- dist_matrix(
      query,
      query,
      mode = mode,
      cost_matrix = cost_matrix,
      gap_cost = gap_cost,
      gap_open_cost = gap_open_cost,
      nthreads = NTHREADS,
      show_progress = FALSE
    )
    max_distance <- rep(max_distance, length.out = length(query))
    include_sizes <- mode %in% c("anchored", "an", "extension", "en")
    query_size <- attr(distances, "query_size")
    target_size <- attr(distances, "target_size")
    rows <- vector("list", length(query))

    for(i in seq_along(query)) {
      idx <- which(seq_along(query) < i & !is.na(distances[i, ]) & distances[i, ] <= max_distance[i])
      if(length(idx) == 0L) next
      row <- data.frame(
        query = query[i],
        target = query[idx],
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
    if(is.null(result)) empty_search_result(include_sizes) else sort_result(result)
  }

  check_identical_result <- function(label, actual, expected) {
    actual <- sort_result(actual)
    expected <- sort_result(expected)
    if(!identical(actual, expected)) {
      print(label)
      print(actual)
      print(expected)
    }
    stopifnot(identical(actual, expected))
  }

  print("Checking direct S7 interface")
  s7_tree <- radix_tree(c("ACGT", "ACGA"))
  stopifnot(is_valid(s7_tree))
  stopifnot(size(s7_tree) == 2L)
  stopifnot(identical(has_sequence(s7_tree, c("ACGT", "AAAA", NA_character_)), c(TRUE, FALSE, FALSE)))
  stopifnot(identical(insert(s7_tree, c("AAAA", "ACGT", NA_character_)), c(TRUE, FALSE, FALSE)))
  stopifnot(identical(erase(s7_tree, c("ACGA", "CCCC", NA_character_)), c(TRUE, FALSE, FALSE)))
  stopifnot(size(s7_tree) == 2L)
  stopifnot(isTRUE(insert(s7_tree, "ACGA")))
  stopifnot(all(c("ACGT", "ACGA") %in% prefix_search(s7_tree, "AC")$target))
  stopifnot(all(c("parent", "child") %in% names(plot_tree(s7_tree, plot = FALSE))))

  r6_tree <- RadixTree$new(to_vector(s7_tree))
  check_identical_result(
    "S7 tree matches R6 tree",
    align_search(s7_tree, c("ACGT", "ACGG"), max_distance = 1L, mode = "levenshtein", nthreads = NTHREADS),
    r6_tree$align_search(c("ACGT", "ACGG"), max_distance = 1L, mode = "levenshtein", nthreads = NTHREADS)
  )

  s7_forest <- radix_forest(c("ACGT", "ACGA", "AAAA"))
  stopifnot(is_valid(s7_forest))
  stopifnot(size(s7_forest) == 3L)
  stopifnot(identical(has_sequence(s7_forest, c("ACGT", "CCCC", NA_character_)), c(TRUE, FALSE, FALSE)))
  stopifnot(all(c("ACGT", "ACGA") %in% prefix_search(s7_forest, "AC")$target))

  r6_forest <- RadixForest$new(to_vector(s7_forest))
  check_identical_result(
    "S7 forest matches R6 forest",
    align_search(s7_forest, c("ACGT", "ACGG"), max_distance = 1L, mode = "levenshtein", nthreads = NTHREADS),
    r6_forest$align_search(c("ACGT", "ACGG"), max_distance = 1L, mode = "levenshtein", nthreads = NTHREADS)
  )

  print("Checking dist_search helpers")
  query <- c("ACGT", "ACGA", "TTTA")
  target <- c("ACGT", "ACGG", "ACGA", "TTTT", "GGGG")

  tree <- RadixTree$new(target)
  check_identical_result(
    "dist_search RadixTree matches direct search",
    dist_search(query, target, max_distance = 1L, mode = "levenshtein", tree_class = "RadixTree", nthreads = NTHREADS),
    tree$search(query, max_distance = 1L, mode = "levenshtein", nthreads = NTHREADS)
  )

  forest <- RadixForest$new(target)
  check_identical_result(
    "dist_search RadixForest matches direct search",
    dist_search(query, target, max_distance = 1L, mode = "levenshtein", tree_class = "RadixForest", nthreads = NTHREADS),
    forest$search(query, max_distance = 1L, mode = "levenshtein", nthreads = NTHREADS)
  )

  cost_matrix <- generate_cost_matrix("ACGT", match = 0L, mismatch = 2L)
  check_identical_result(
    "dist_search RadixForest linear custom cost",
    dist_search(query, target, max_distance = 2L, mode = "levenshtein", cost_matrix = cost_matrix, gap_cost = 1L, tree_class = "RadixForest", nthreads = NTHREADS),
    forest$search(query, max_distance = 2L, mode = "levenshtein", cost_matrix = cost_matrix, gap_cost = 1L, nthreads = NTHREADS)
  )

  self_query <- c("AAAA", "AAAT", "AATT", "ATTT", "TTTT")
  self_cost_matrix <- generate_cost_matrix("AT", match = 0L, mismatch = 2L)
  self_cases <- list(
    list(label = "hamming", mode = "hamming", max_distance = 1L, cost_matrix = NULL, gap_cost = NA_integer_, gap_open_cost = NA_integer_),
    list(label = "global unit", mode = "levenshtein", max_distance = c(0L, 1L, 2L, 1L, 0L), cost_matrix = NULL, gap_cost = NA_integer_, gap_open_cost = NA_integer_),
    list(label = "global linear", mode = "levenshtein", max_distance = 3L, cost_matrix = self_cost_matrix, gap_cost = 2L, gap_open_cost = NA_integer_),
    list(label = "global affine", mode = "levenshtein", max_distance = 4L, cost_matrix = self_cost_matrix, gap_cost = 1L, gap_open_cost = 2L)
  )

  for(case in self_cases) {
    for(tree_class in c("RadixTree", "RadixForest")) {
      check_identical_result(
        paste("dist_search self", case$label, tree_class),
        dist_search(
          self_query,
          max_distance = case$max_distance,
          mode = case$mode,
          cost_matrix = case$cost_matrix,
          gap_cost = case$gap_cost,
          gap_open_cost = case$gap_open_cost,
          tree_class = tree_class,
          nthreads = NTHREADS,
          show_progress = FALSE
        ),
        expected_self_search(
          self_query,
          max_distance = case$max_distance,
          mode = case$mode,
          cost_matrix = case$cost_matrix,
          gap_cost = case$gap_cost,
          gap_open_cost = case$gap_open_cost
        )
      )
    }
  }

  anchored_cases <- list(
    list(label = "anchored unit", mode = "anchored", max_distance = c(0L, 1L, 2L, 1L, 0L), cost_matrix = NULL, gap_cost = NA_integer_, gap_open_cost = NA_integer_),
    list(label = "anchored linear", mode = "anchored", max_distance = 3L, cost_matrix = self_cost_matrix, gap_cost = 2L, gap_open_cost = NA_integer_),
    list(label = "anchored affine", mode = "anchored", max_distance = 4L, cost_matrix = self_cost_matrix, gap_cost = 1L, gap_open_cost = 2L)
  )

  for(case in anchored_cases) {
    check_identical_result(
      paste("dist_search self", case$label),
      dist_search(
        self_query,
        max_distance = case$max_distance,
        mode = case$mode,
        cost_matrix = case$cost_matrix,
        gap_cost = case$gap_cost,
        gap_open_cost = case$gap_open_cost,
        tree_class = "RadixTree",
        nthreads = NTHREADS,
        show_progress = FALSE
      ),
      expected_self_search(
        self_query,
        max_distance = case$max_distance,
        mode = case$mode,
        cost_matrix = case$cost_matrix,
        gap_cost = case$gap_cost,
        gap_open_cost = case$gap_open_cost
      )
    )
  }

  check_identical_result(
    "dist_search self duplicate behavior",
    dist_search(c("AAAA", "AAAA", "AAAT"), max_distance = 1L, mode = "hamming", tree_class = "RadixTree"),
    data.frame(
      query = c("AAAA", "AAAT"),
      target = c("AAAA", "AAAA"),
      distance = c(0L, 1L),
      stringsAsFactors = FALSE
    )
  )

  print("Checking split_search")
  reverse_strings <- function(x) {
    vapply(strsplit(x, "", fixed = TRUE), function(chars) paste0(rev(chars), collapse = ""), character(1L))
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

    keep_target <- !duplicated(paste(target, target_split, sep = "\r"))
    target <- target[keep_target]
    target_split <- target_split[keep_target]

    query_parts <- split_parts(query, query_split, edge_trim)
    target_parts <- split_parts(target, target_split, edge_trim)

    left_dist <- dist_matrix(
      query_parts$left,
      target_parts$left,
      mode = "anchored",
      cost_matrix = cost_matrix,
      gap_cost = gap_cost,
      gap_open_cost = gap_open_cost,
      nthreads = NTHREADS
    )
    right_dist <- dist_matrix(
      query_parts$right,
      target_parts$right,
      mode = "anchored",
      cost_matrix = cost_matrix,
      gap_cost = gap_cost,
      gap_open_cost = gap_open_cost,
      nthreads = NTHREADS
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
    if(is.null(out)) empty_search_result(FALSE) else sort_result(out)
  }

  check_split_search <- function(label, query, target, query_split, target_split, ...) {
    actual <- sort_result(split_search(query, target, query_split, target_split, nthreads = NTHREADS, ...))
    expected <- split_search_expected(query, target, query_split, target_split, ...)
    if(!identical(actual, expected)) {
      print(label)
      print(actual)
      print(expected)
    }
    stopifnot(identical(actual, expected))
  }

  random_strings <- function(n, charset = "ACGT", min_len = 4L, max_len = 9L) {
    chars <- strsplit(charset, "", fixed = TRUE)[[1L]]
    len <- sample(min_len:max_len, n, replace = TRUE)
    vapply(len, function(x) paste0(sample(chars, x, replace = TRUE), collapse = ""), character(1L))
  }

  random_split <- function(sequences, edge_trim) {
    vapply(sequences, function(x) {
      if(runif(1) < 0.2) return(-1L)
      sample(seq.int(edge_trim + 1L, nchar(x) - edge_trim), 1L)
    }, integer(1L))
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
    "same target sequence with distinct split positions",
    query = "ABCD",
    target = c("ABCD", "ABCD"),
    query_split = 3L,
    target_split = c(2L, 3L),
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

  cost_matrix <- generate_cost_matrix("ACGT", match = 0L, mismatch = 2L)
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

  nonzero_diag <- generate_cost_matrix("ACGTN", match = 0L, mismatch = 2L)
  nonzero_diag["N", "N"] <- 2L
  for(tree_class in c("RadixTree", "RadixForest")) {
    check_identical_result(
      paste("nonzero diagonal cost matrix", tree_class),
      dist_search(
        query = c("NN", "NA", "AA"),
        target = c("NN", "NA", "AA"),
        max_distance = 2L,
        mode = "levenshtein",
        cost_matrix = nonzero_diag,
        gap_cost = 1L,
        tree_class = tree_class,
        nthreads = NTHREADS,
        show_progress = FALSE
      ),
      expected_search(
        query = c("NN", "NA", "AA"),
        target = c("NN", "NA", "AA"),
        max_distance = 2L,
        mode = "levenshtein",
        cost_matrix = nonzero_diag,
        gap_cost = 1L
      )
    )
  }

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

  print("Checking utility validation")
  amb <- generate_cost_matrix("ACGT", ambiguity_base = "N", match = 0L, mismatch = 2L)
  stopifnot(all(amb["N", c("A", "C", "G", "T")] == 0L))
  stopifnot(all(amb[c("A", "C", "G", "T"), "N"] == 0L))
  stopifnot(inherits(try(generate_cost_matrix("ACGT", ambiguity_base = "NN"), silent = TRUE), "try-error"))
  stopifnot(inherits(try(generate_cost_matrix("ACGT", match = 0.5, mismatch = 1L), silent = TRUE), "try-error"))
  stopifnot(inherits(try(dist_matrix("A", "A", mode = "bad"), silent = TRUE), "try-error"))
  stopifnot(inherits(try(dist_search("A", "A", max_distance = 1L, nthreads = 0L), silent = TRUE), "try-error"))

  bad_cost_matrix <- matrix(0L, nrow = 1L, ncol = 1L, dimnames = list("A", "A"))
  stopifnot(inherits(
    try(dist_matrix("C", "A", mode = "levenshtein", cost_matrix = bad_cost_matrix, gap_cost = 1L), silent = TRUE),
    "try-error"
  ))

  warned <- FALSE
  withCallingHandlers(
    dist_matrix("A", "A", mode = "hamming", cost_matrix = generate_cost_matrix("A"), gap_cost = 1L),
    warning = function(w) {
      warned <<- TRUE
      invokeRestart("muffleWarning")
    }
  )
  stopifnot(warned)
}

print(Sys.time() - runtime)
