print("Running test_RadixTree.R")
runtime <- Sys.time()

if(requireNamespace("seqtrie", quietly = TRUE)) {
  library(seqtrie)

  IS_LOCAL <- Sys.getenv("IS_LOCAL") != ""
  NTHREADS <- ifelse(IS_LOCAL, 4L, 2L)
  NITER <- ifelse(IS_LOCAL, 3L, 1L)
  NTARGET <- ifelse(IS_LOCAL, 10000L, 2500L)
  NQUERY <- ifelse(IS_LOCAL, 100L, 25L)
  MAXSEQLEN <- 80L
  CHARSET <- "ACGT"

  test_seed <- Sys.getenv("SEQTRIE_TEST_SEED")
  if(nzchar(test_seed)) {
    test_seed <- as.integer(test_seed)
  } else {
    test_seed <- as.integer(as.numeric(Sys.time())) %% .Machine$integer.max
  }
  cat("Test seed:", test_seed, "\n")
  set.seed(test_seed)

  charset_chars <- function(charset = CHARSET) {
    strsplit(charset, "", fixed = TRUE)[[1L]]
  }

  random_strings <- function(n, charset = CHARSET, min_len = 0L, max_len = MAXSEQLEN) {
    chars <- charset_chars(charset)
    len <- sample(min_len:max_len, n, replace = TRUE)
    vapply(len, function(x) paste0(sample(chars, x, replace = TRUE), collapse = ""), character(1L))
  }

  mutate_strings <- function(x, prob = 0.025, indel_prob = 0.025, charset = CHARSET) {
    chars <- charset_chars(charset)
    xsplit <- strsplit(x, "", fixed = TRUE)
    vapply(xsplit, function(a) {
      if(length(a) == 0L) return("")
      r <- runif(length(a)) < prob
      a[r] <- sample(chars, sum(r), replace = TRUE)
      ins <- runif(length(a)) < indel_prob
      a[ins] <- paste0(sample(chars, sum(ins), replace = TRUE), sample(chars, sum(ins), replace = TRUE))
      del <- runif(length(a)) < indel_prob
      a[del] <- ""
      paste0(a, collapse = "")
    }, character(1L))
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

  expected_from_matrix <- function(query,
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
    if(is.null(result)) empty_search_result(include_sizes) else sort_search_result(result)
  }

  check_identical_result <- function(label, actual, expected) {
    actual <- sort_search_result(actual)
    expected <- sort_search_result(expected)
    if(!identical(actual, expected)) {
      print(label)
      print(actual)
      print(expected)
    }
    stopifnot(identical(actual, expected))
  }

  random_cost_matrix <- function(charset = CHARSET, max_cost = 4L) {
    chars <- charset_chars(charset)
    cost_matrix <- matrix(
      sample.int(max_cost, length(chars) * length(chars), replace = TRUE),
      nrow = length(chars),
      dimnames = list(chars, chars)
    )
    diag(cost_matrix) <- 0L
    cost_matrix
  }

  make_active_target <- function() {
    inserted <- unique(c(random_strings(NTARGET, CHARSET), ""))
    erase_existing <- sample(inserted, min(length(inserted), max(1L, length(inserted) %/% 8L)))
    erase_missing <- random_strings(max(4L, length(inserted) %/% 20L), CHARSET)
    erase <- unique(c(erase_existing, erase_missing))
    active <- setdiff(inserted, erase)
    list(inserted = inserted, erase = erase, active = active)
  }

  make_query <- function(active, hamming = FALSE) {
    seed <- sample(active, min(length(active), NQUERY), replace = length(active) < NQUERY)
    indel_prob <- if(hamming) 0 else 0.04
    unique(c(seed, mutate_strings(seed, prob = 0.06, indel_prob = indel_prob, charset = CHARSET), ""))
  }

  check_tree_state_after_erase <- function(tree, inserted, erase, active) {
    stopifnot(tree$validate())
    stopifnot(tree$size() == length(unique(active)))
    stopifnot(setequal(tree$to_vector(), unique(active)))

    find_query <- unique(c(sample(active, min(30L, length(active))), random_strings(30L, CHARSET), ""))
    stopifnot(identical(tree$find(find_query), find_query %in% active))

    prefix_query <- unique(substr(c(sample(active, min(20L, length(active))), ""), 1L, 5L))
    expected <- lapply(prefix_query, function(prefix) {
      target <- active[startsWith(active, prefix)]
      if(length(target) == 0L) return(NULL)
      data.frame(query = prefix, target = target, stringsAsFactors = FALSE)
    })
    expected <- do.call(rbind, expected)
    if(is.null(expected)) {
      expected <- data.frame(query = character(), target = character(), stringsAsFactors = FALSE)
    }
    expected <- sort_search_result(expected)
    actual <- sort_search_result(tree$prefix_search(prefix_query))
    stopifnot(identical(actual, expected))

    if(length(erase) > 0L) {
      stopifnot(!any(tree$find(setdiff(erase, active))))
    }
    invisible(TRUE)
  }

  check_search_case <- function(label,
                                tree,
                                query,
                                target,
                                mode,
                                max_distance,
                                cost_matrix = NULL,
                                gap_cost = NA_integer_,
                                gap_open_cost = NA_integer_) {
    actual <- tree$search(
      query,
      max_distance = max_distance,
      mode = mode,
      cost_matrix = cost_matrix,
      gap_cost = gap_cost,
      gap_open_cost = gap_open_cost,
      nthreads = NTHREADS,
      show_progress = FALSE
    )
    expected <- expected_from_matrix(
      query,
      target,
      max_distance = max_distance,
      mode = mode,
      cost_matrix = cost_matrix,
      gap_cost = gap_cost,
      gap_open_cost = gap_open_cost
    )
    check_identical_result(label, actual, expected)
  }

  check_search_fraction <- function(label,
                                    tree,
                                    query,
                                    target,
                                    mode,
                                    max_fraction,
                                    cost_matrix = NULL,
                                    gap_cost = NA_integer_,
                                    gap_open_cost = NA_integer_) {
    actual <- tree$search(
      query,
      max_fraction = max_fraction,
      mode = mode,
      cost_matrix = cost_matrix,
      gap_cost = gap_cost,
      gap_open_cost = gap_open_cost,
      nthreads = NTHREADS,
      show_progress = FALSE
    )
    expected <- expected_from_matrix(
      query,
      target,
      max_distance = as.integer(nchar(query) * max_fraction),
      mode = mode,
      cost_matrix = cost_matrix,
      gap_cost = gap_cost,
      gap_open_cost = gap_open_cost
    )
    check_identical_result(label, actual, expected)
  }

  check_child_map_transitions <- function(prefix) {
    suffixes <- c("A", "C", "G", "T", "N", "X", "Y", "Z")
    sequences <- paste0(prefix, suffixes)
    tree <- RadixTree$new()
    expected <- character()

    for(seq in sequences) {
      stopifnot(identical(tree$insert(seq), TRUE))
      expected <- c(expected, seq)
      stopifnot(tree$validate())
      stopifnot(tree$size() == length(expected))
      stopifnot(setequal(tree$to_vector(), expected))
      stopifnot(all(tree$find(expected)))
    }

    for(seq in sequences) {
      stopifnot(identical(tree$erase(seq), TRUE))
      expected <- setdiff(expected, seq)
      stopifnot(tree$validate())
      stopifnot(tree$size() == length(expected))
      stopifnot(setequal(tree$to_vector(), expected))
      if(length(expected) > 0L) stopifnot(all(tree$find(expected)))
      stopifnot(!tree$find(seq))
    }
  }

  print("Checking RadixTree insert/erase/search fuzz cases")
  for(iter in seq_len(NITER)) {
    data <- make_active_target()
    tree <- RadixTree$new()
    inserted <- data$inserted
    active <- data$active
    tree$insert(inserted)
    stopifnot(tree$validate())
    tree$erase(data$erase)
    check_tree_state_after_erase(tree, inserted, data$erase, active)

    unit_query <- make_query(active)
    hamming_query <- make_query(active, hamming = TRUE)
    cost_matrix <- random_cost_matrix(CHARSET)
    gap_cost <- sample.int(3L, 1L)
    affine_gap_cost <- sample.int(2L, 1L)
    gap_open_cost <- affine_gap_cost + sample.int(3L, 1L)

    cases <- list(
      list(label = "hamming", query = hamming_query, mode = "hamming", max_distance = sample(0:4, length(hamming_query), replace = TRUE)),
      list(label = "global unit", query = unit_query, mode = "levenshtein", max_distance = sample(0:8, length(unit_query), replace = TRUE)),
      list(label = "anchored unit", query = unit_query, mode = "anchored", max_distance = sample(0:8, length(unit_query), replace = TRUE)),
      list(label = "global linear", query = unit_query, mode = "levenshtein", max_distance = 8L, cost_matrix = cost_matrix, gap_cost = gap_cost),
      list(label = "anchored linear", query = unit_query, mode = "anchored", max_distance = 8L, cost_matrix = cost_matrix, gap_cost = gap_cost),
      list(label = "global affine", query = unit_query, mode = "levenshtein", max_distance = 10L, cost_matrix = cost_matrix, gap_cost = affine_gap_cost, gap_open_cost = gap_open_cost),
      list(label = "anchored affine", query = unit_query, mode = "anchored", max_distance = 10L, cost_matrix = cost_matrix, gap_cost = affine_gap_cost, gap_open_cost = gap_open_cost)
    )

    for(case in cases) {
      check_search_case(
        paste("fuzz", iter, case$label),
        tree = tree,
        query = case$query,
        target = active,
        mode = case$mode,
        max_distance = case$max_distance,
        cost_matrix = case$cost_matrix,
        gap_cost = if(is.null(case$gap_cost)) NA_integer_ else case$gap_cost,
        gap_open_cost = if(is.null(case$gap_open_cost)) NA_integer_ else case$gap_open_cost
      )
    }

    check_search_fraction(
      paste("fuzz", iter, "max_fraction global"),
      tree = tree,
      query = unit_query,
      target = active,
      mode = "levenshtein",
      max_fraction = 0.08
    )
    check_search_fraction(
      paste("fuzz", iter, "max_fraction anchored"),
      tree = tree,
      query = unit_query,
      target = active,
      mode = "anchored",
      max_fraction = 0.08
    )
  }

  print("Checking RadixTree explicit edge cases")
  subtree_target <- c("", "A", "AA", "AAA", "AC", "CAA", "T")
  subtree <- RadixTree$new(subtree_target)
  check_identical_result(
    "anchored empty query should add the whole subtree",
    subtree$search("", max_distance = 0L, mode = "anchored", nthreads = NTHREADS, show_progress = FALSE),
    data.frame(
      query = rep("", length(subtree_target)),
      target = subtree_target,
      distance = 0L,
      query_size = 0L,
      target_size = 0L,
      stringsAsFactors = FALSE
    )
  )
  check_identical_result(
    "anchored prefix match should add terminal child descendants",
    subtree$search("A", max_distance = 0L, mode = "anchored", nthreads = NTHREADS, show_progress = FALSE),
    data.frame(
      query = rep("A", 5L),
      target = c("", "A", "AA", "AAA", "AC"),
      distance = 0L,
      query_size = c(0L, 1L, 1L, 1L, 1L),
      target_size = c(0L, 1L, 1L, 1L, 1L),
      stringsAsFactors = FALSE
    )
  )

  check_identical_result(
    "hamming should accumulate distance across compressed branches",
    RadixTree$new(c("AAAAC", "AAAAT", "AAAGT", "CCCC", "CCCT"))$search(
      c("AAAAG", "CCCC", "CCCA", "AAAGT"),
      max_distance = c(1L, 0L, 1L, 0L),
      mode = "hamming",
      nthreads = NTHREADS,
      show_progress = FALSE
    ),
    data.frame(
      query = c("AAAAG", "AAAAG", "AAAGT", "CCCC", "CCCA", "CCCA"),
      target = c("AAAAC", "AAAAT", "AAAGT", "CCCC", "CCCC", "CCCT"),
      distance = c(1L, 1L, 0L, 0L, 1L, 1L),
      stringsAsFactors = FALSE
    )
  )

  filter_tree <- RadixTree$new(c("AAAA", "AAAT", "AATT", "TTTT"))
  check_identical_result(
    "match_mode best should keep all nearest ties",
    filter_tree$search(
      c("AAAG", "TTTT"),
      max_distance = 3L,
      match_mode = "best",
      nthreads = NTHREADS,
      show_progress = FALSE
    ),
    data.frame(
      query = c("AAAG", "AAAG", "TTTT"),
      target = c("AAAA", "AAAT", "TTTT"),
      distance = c(1L, 1L, 0L),
      stringsAsFactors = FALSE
    )
  )
  check_identical_result(
    "lower_triangle should compare original query index to target insertion index",
    filter_tree$search(
      c(NA_character_, "AAAT", "AATT"),
      max_distance = 2L,
      lower_triangle = TRUE,
      nthreads = NTHREADS,
      show_progress = FALSE
    ),
    data.frame(
      query = c("AAAT", "AATT", "AATT"),
      target = c("AAAA", "AAAA", "AAAT"),
      distance = c(1L, 2L, 1L),
      stringsAsFactors = FALSE
    )
  )

  target <- c("", "A", "AA", "AAA", "AC", "CAA", "T", "TTTT", "ACCC")
  query <- c("", "A", "AA", "AAC", "TT", "CCCC")
  max_distance <- c(0L, 0L, 1L, 2L, 2L, 3L)
  edge_cost_matrix <- matrix(
    c(0L, 3L, 4L,
      2L, 0L, 3L,
      4L, 2L, 0L),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(c("A", "C", "T"), c("A", "C", "T"))
  )
  edge_tree <- RadixTree$new(target)
  check_search_case("global unit empty strings", edge_tree, query, target, "levenshtein", max_distance)
  check_search_case("anchored unit empty strings", edge_tree, query, target, "anchored", max_distance)
  check_search_case("global linear empty strings", edge_tree, query, target, "levenshtein", max_distance + 2L, edge_cost_matrix, 2L)
  check_search_case("anchored linear empty strings", edge_tree, query, target, "anchored", max_distance + 2L, edge_cost_matrix, 2L)
  check_search_case("global affine empty strings", edge_tree, query, target, "levenshtein", max_distance + 3L, edge_cost_matrix, 1L, 2L)
  check_search_case("anchored affine empty strings", edge_tree, query, target, "anchored", max_distance + 3L, edge_cost_matrix, 1L, 2L)

  check_child_map_transitions("")
  check_child_map_transitions("A")

  print("Checking RadixTree affine boundary cases")
  enumerate_strings <- function(chars, max_len) {
    out <- ""
    for(len in seq_len(max_len)) {
      grid <- expand.grid(rep(list(chars), len), stringsAsFactors = FALSE)
      out <- c(out, apply(grid, 1L, paste0, collapse = ""))
    }
    unique(out)
  }
  high_mismatch_cost <- function(chars, mismatch) {
    cost_matrix <- matrix(
      as.integer(mismatch),
      nrow = length(chars),
      ncol = length(chars),
      dimnames = list(chars, chars)
    )
    diag(cost_matrix) <- 0L
    cost_matrix
  }
  asymmetric_cost <- function(chars) {
    n <- length(chars)
    values <- ((seq_len(n * n) * 3L) %% 7L) + 1L
    cost_matrix <- matrix(as.integer(values), nrow = n, ncol = n)
    diag(cost_matrix) <- 0L
    dimnames(cost_matrix) <- list(chars, chars)
    cost_matrix
  }
  run_affine_boundary_case <- function(label, query, target, cost_matrix, gap_cost, gap_open_cost, thresholds) {
    tree <- RadixTree$new(unique(target))
    for(mode in c("levenshtein", "anchored")) {
      for(name in names(thresholds)) {
        check_search_case(
          paste(label, mode, name, sep = " / "),
          tree,
          unique(query),
          unique(target),
          mode,
          thresholds[[name]],
          cost_matrix,
          gap_cost,
          gap_open_cost
        )
      }
    }
  }

  small_strings <- enumerate_strings(c("A", "C"), 5L)
  run_affine_boundary_case(
    "small exhaustive high mismatch",
    small_strings,
    small_strings,
    high_mismatch_cost(c("A", "C"), 5L),
    gap_cost = 1L,
    gap_open_cost = 3L,
    thresholds = list(scalar_tight = 2L, scalar_wide = 5L, per_query = as.integer((seq_along(small_strings) - 1L) %% 6L))
  )
  run_affine_boundary_case(
    "small exhaustive asymmetric",
    small_strings,
    small_strings,
    asymmetric_cost(c("A", "C")),
    gap_cost = 2L,
    gap_open_cost = 1L,
    thresholds = list(scalar_tight = 2L, scalar_wide = 5L, per_query = as.integer((seq_along(small_strings) - 1L) %% 6L))
  )

  boundary_target <- c(
    "", "A", "AA", "AAA", "AAAA", "AAAAAA", "AAAAAAAA",
    "CAAAAAAAA", "AAAAAAAAC", "AAAACAAAA", "C", "CC", "CCCC",
    "CCCCCCCC", "CCCCAAAACCCC", "ACACACAC", "ACACACACAC",
    "TAAAAAAAAT", "GCCCCCCCG", "AAAAGGGGCCCC"
  )
  boundary_query <- c(
    "", "A", "AAAA", "AAAAAAA", "AAAAAAAAAA", "CAAAAAAAAC",
    "AAAAACAAAA", "AAAACAAA", "C", "CCCCC", "CCCCCCCCCC",
    "CCCCAACCCC", "ACACAC", "ACACACACACAC", "TAAAAT",
    "TAAAAAAAAAAT", "GCCCCG", "AAAAGGCCCC"
  )
  boundary_thresholds <- list(
    max_1 = 1L,
    max_3 = 3L,
    max_6 = 6L,
    per_query_boundary = as.integer(pmin(8L, 1L + (nchar(boundary_query) %% 8L)))
  )
  run_affine_boundary_case(
    "constructed boundary high mismatch",
    boundary_query,
    boundary_target,
    high_mismatch_cost(c("A", "C", "G", "T"), 6L),
    gap_cost = 1L,
    gap_open_cost = 2L,
    thresholds = boundary_thresholds
  )
  run_affine_boundary_case(
    "constructed boundary asymmetric",
    boundary_query,
    boundary_target,
    asymmetric_cost(c("A", "C", "G", "T")),
    gap_cost = 2L,
    gap_open_cost = 2L,
    thresholds = boundary_thresholds
  )

  print("Checking RadixTree single_gap_search")
  cost_mat <- generate_cost_matrix(charset = "ACGT", match = 0L, mismatch = 1L)
  drop_size_cols <- function(results) {
    results[, setdiff(names(results), c("query_size", "target_size")), drop = FALSE]
  }

  target <- c("ACGT", "ACGGT", "AACGT", "ACGTAA", "GGGG")
  query <- c("ACGT", "ACGGT", "ACGTAA", "GGGG")
  max_distance <- c(2L, 3L, 1L, 0L)
  tree <- RadixTree$new(target)
  single_gap <- sort_search_result(tree$single_gap_search(query, max_distance = max_distance, gap_cost = 2L, nthreads = NTHREADS, show_progress = FALSE))
  expected_rows <- lapply(seq_along(query), function(i) {
    result <- expected_from_matrix(query[i], target, max_distance[i], "anchored", cost_mat, gap_cost = 2L)
    drop_size_cols(result)
  })
  expected <- sort_search_result(do.call(rbind, expected_rows))
  stopifnot(identical(single_gap, expected))

  tree <- RadixTree$new("ABEFGHIJKLMNO")
  result <- tree$single_gap_search("ABCDEFGHIJKLMNO", max_distance = 5L, gap_cost = 2L, nthreads = NTHREADS, show_progress = FALSE)
  stopifnot(nrow(result) == 0L)
  result <- tree$single_gap_search("ABCDEFGHIJKLMNO", max_distance = 100L, gap_cost = 2L, nthreads = NTHREADS, show_progress = FALSE)
  stopifnot(identical(result$distance, nchar("ABEFGHIJKLMNO") - 2L))

  seqs <- gsub("G|T", "A", sample(covid_cdr3, 1000L))
  tree <- RadixTree$new(seqs)
  result <- sort_search_result(tree$single_gap_search(seqs, max_distance = 8L, gap_cost = 5L, nthreads = NTHREADS, show_progress = FALSE))
  expected <- sort_search_result(drop_size_cols(tree$search(seqs, cost_matrix = cost_mat, max_distance = 8L, gap_cost = 5L, mode = "anchored", nthreads = NTHREADS, show_progress = FALSE)))
  stopifnot(identical(result, expected))

  tree <- RadixTree$new("GACCC")
  stopifnot(identical(
    sort_search_result(tree$single_gap_search("GTT", max_distance = 3L, gap_cost = 1L, nthreads = NTHREADS, show_progress = FALSE)),
    sort_search_result(drop_size_cols(tree$search("GTT", max_distance = 3L, cost_matrix = cost_mat, gap_cost = 1L, mode = "anchored", nthreads = NTHREADS, show_progress = FALSE)))
  ))
}

print(Sys.time() - runtime)
