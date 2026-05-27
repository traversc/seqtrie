print("Running test_RadixForest.R")
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
    columns <- intersect(c("query", "target", "distance"), names(result))
    result <- result[, columns, drop = FALSE]
    if(nrow(result) > 0L) {
      result <- result[do.call(order, result[, columns, drop = FALSE]), , drop = FALSE]
    }
    rownames(result) <- NULL
    result
  }

  empty_search_result <- function() {
    data.frame(query = character(), target = character(), distance = integer(), stringsAsFactors = FALSE)
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

    result <- do.call(rbind, rows)
    if(is.null(result)) empty_search_result() else sort_search_result(result)
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

  check_forest_state_after_erase <- function(forest, erase, active) {
    stopifnot(forest$validate())
    stopifnot(forest$size() == length(unique(active)))
    stopifnot(setequal(forest$to_vector(), unique(active)))

    find_query <- unique(c(sample(active, min(30L, length(active))), random_strings(30L, CHARSET), ""))
    stopifnot(identical(forest$find(find_query), find_query %in% active))

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
    actual <- sort_search_result(forest$prefix_search(prefix_query))
    stopifnot(identical(actual, expected))

    if(length(erase) > 0L) {
      stopifnot(!any(forest$find(setdiff(erase, active))))
    }
    invisible(TRUE)
  }

  check_search_case <- function(label,
                                forest,
                                query,
                                target,
                                mode,
                                max_distance,
                                cost_matrix = NULL,
                                gap_cost = NA_integer_,
                                gap_open_cost = NA_integer_) {
    actual <- forest$search(
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
                                    forest,
                                    query,
                                    target,
                                    mode,
                                    max_fraction,
                                    cost_matrix = NULL,
                                    gap_cost = NA_integer_,
                                    gap_open_cost = NA_integer_) {
    actual <- forest$search(
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

  print("Checking RadixForest insert/erase/search fuzz cases")
  for(iter in seq_len(NITER)) {
    data <- make_active_target()
    forest <- RadixForest$new()
    active <- data$active
    forest$insert(data$inserted)
    stopifnot(forest$validate())
    forest$erase(data$erase)
    check_forest_state_after_erase(forest, data$erase, active)

    unit_query <- make_query(active)
    hamming_query <- make_query(active, hamming = TRUE)
    cost_matrix <- random_cost_matrix(CHARSET)
    gap_cost <- sample.int(3L, 1L)
    affine_gap_cost <- sample.int(2L, 1L)
    gap_open_cost <- affine_gap_cost + sample.int(3L, 1L)

    cases <- list(
      list(label = "hamming", query = hamming_query, mode = "hamming", max_distance = sample(0:4, length(hamming_query), replace = TRUE)),
      list(label = "global unit", query = unit_query, mode = "levenshtein", max_distance = sample(0:8, length(unit_query), replace = TRUE)),
      list(label = "global linear", query = unit_query, mode = "levenshtein", max_distance = 8L, cost_matrix = cost_matrix, gap_cost = gap_cost),
      list(label = "global affine", query = unit_query, mode = "levenshtein", max_distance = 10L, cost_matrix = cost_matrix, gap_cost = affine_gap_cost, gap_open_cost = gap_open_cost)
    )

    for(case in cases) {
      check_search_case(
        paste("fuzz", iter, case$label),
        forest = forest,
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
      forest = forest,
      query = unit_query,
      target = active,
      mode = "levenshtein",
      max_fraction = 0.08
    )
  }

  print("Checking RadixForest custom-cost and edge cases")
  target <- c("ACGT", "ACGGT", "AGT", "TTTT", "ACCCCCGT", "")
  query <- c("ACGT", "ACT", "", NA_character_)
  query_keep <- query[!is.na(query)]
  max_distance <- c(4L, 5L, 2L, 100L)
  max_distance_keep <- max_distance[!is.na(query)]
  cost_matrix <- generate_cost_matrix("ACGT", match = 0L, mismatch = 2L)
  forest <- RadixForest$new(target)
  check_search_case("linear custom cost with empty and NA query", forest, query_keep, target, "global", max_distance_keep, cost_matrix, 2L)
  check_search_case("affine custom cost with empty and NA query", forest, query_keep, target, "global", max_distance_keep, cost_matrix, 1L, 3L)

  edge_target <- c("", "A", "AA", "AAA", "AC", "CAA", "T", "TTTT", "ACCC")
  edge_query <- c("", "A", "AA", "AAC", "TT", "CCCC")
  edge_max_distance <- c(0L, 0L, 1L, 2L, 2L, 3L)
  edge_forest <- RadixForest$new(edge_target)
  check_search_case("hamming edge cases", edge_forest, edge_query, edge_target, "hamming", edge_max_distance)
  check_search_case("global edge cases", edge_forest, edge_query, edge_target, "levenshtein", edge_max_distance)

  filter_forest <- RadixForest$new(c("AAAA", "AAAT", "AATT", "TTTT"))
  check_identical_result(
    "match_mode best should keep all nearest ties",
    filter_forest$search(
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
    filter_forest$search(
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

  anchored_error <- tryCatch({
    forest$search("ACGT", max_distance = 1L, mode = "anchored")
    FALSE
  }, error = function(e) TRUE)
  stopifnot(anchored_error)
}

print(Sys.time() - runtime)
