print("Running test_StarTree.R")
runtime <- Sys.time()

if(requireNamespace("seqtrie", quietly = TRUE)) {
  library(seqtrie)

  IS_LOCAL <- Sys.getenv("IS_LOCAL") != ""
  NTHREADS <- ifelse(IS_LOCAL, 4L, 2L)

  sort_result <- function(x) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    x <- x[, c("query", "target", "distance"), drop = FALSE]
    if(nrow(x) > 0L) {
      x <- x[do.call(order, x), , drop = FALSE]
    }
    rownames(x) <- NULL
    x
  }

  canon_pairs <- function(x) {
    x <- sort_result(x)
    if(nrow(x) == 0L) {
      return(character())
    }
    sort(paste(pmin(x$query, x$target), pmax(x$query, x$target), x$distance, sep = "\t"))
  }

  expect_error_message <- function(expr, pattern) {
    err <- tryCatch(
      {
        force(expr)
        NULL
      },
      error = identity
    )
    if(is.null(err)) {
      stop("expected an error matching: ", pattern)
    }
    if(!grepl(pattern, conditionMessage(err))) {
      stop("unexpected error message: ", conditionMessage(err))
    }
    invisible(TRUE)
  }

  normalize_unique <- function(x) {
    x <- toupper(x)
    x <- unique(x)
    x[order(nchar(x), x)]
  }

  weighted_distance <- function(target, query, mismatch_cost = 1L, gap_cost = 1L) {
    target <- strsplit(target, "", fixed = TRUE)[[1L]]
    query <- strsplit(query, "", fixed = TRUE)[[1L]]
    n <- length(target)
    m <- length(query)
    dp <- matrix(0L, nrow = n + 1L, ncol = m + 1L)
    dp[, 1L] <- seq.int(0L, n) * gap_cost
    dp[1L, ] <- seq.int(0L, m) * gap_cost
    for(i in seq_len(n)) {
      for(j in seq_len(m)) {
        subst <- if(target[i] == query[j] && target[i] != "N") 0L else mismatch_cost
        dp[i + 1L, j + 1L] <- min(
          dp[i, j] + subst,
          dp[i, j + 1L] + gap_cost,
          dp[i + 1L, j] + gap_cost
        )
      }
    }
    dp[n + 1L, m + 1L]
  }

  expected_self <- function(sequences, max_distance, mismatch_cost = 1L, gap_cost = 1L) {
    target <- normalize_unique(sequences)
    rows <- list()
    k <- 0L
    for(i in seq_along(target)) {
      if(i == 1L) next
      for(j in seq_len(i - 1L)) {
        d <- weighted_distance(target[j], target[i], mismatch_cost, gap_cost)
        if(d > 0L && d <= max_distance) {
          k <- k + 1L
          rows[[k]] <- data.frame(query = target[i], target = target[j], distance = d, stringsAsFactors = FALSE)
        }
      }
    }
    out <- do.call(rbind, rows)
    if(is.null(out)) {
      out <- data.frame(query = character(), target = character(), distance = integer(), stringsAsFactors = FALSE)
    }
    sort_result(out)
  }

  expected_search <- function(query, target, max_distance, mismatch_cost = 1L, gap_cost = 1L) {
    query <- toupper(query[!is.na(query)])
    query <- query[order(nchar(query), query)]
    target <- normalize_unique(target)
    rows <- list()
    k <- 0L
    for(q in query) {
      for(t in target) {
        d <- weighted_distance(t, q, mismatch_cost, gap_cost)
        if(d <= max_distance) {
          k <- k + 1L
          rows[[k]] <- data.frame(query = q, target = t, distance = d, stringsAsFactors = FALSE)
        }
      }
    }
    out <- do.call(rbind, rows)
    if(is.null(out)) {
      out <- data.frame(query = character(), target = character(), distance = integer(), stringsAsFactors = FALSE)
    }
    sort_result(out)
  }

  random_dna <- function(n, min_len = 1L, max_len = 14L) {
    lens <- sample(seq.int(min_len, max_len), n, replace = TRUE)
    out <- character(n)
    for(len in unique(lens)) {
      idx <- which(lens == len)
      chars <- matrix(sample(c("A", "C", "G", "T", "N"), length(idx) * len, replace = TRUE),
                      nrow = len, ncol = length(idx))
      out[idx] <- apply(chars, 2L, paste0, collapse = "")
    }
    out
  }

  unique_random_dna <- function(n, min_len = 1L, max_len = 14L) {
    out <- character()
    while(length(out) < n) {
      out <- unique(c(out, random_dna(max(2L * n, 8L), min_len, max_len)))
    }
    out[seq_len(n)]
  }

  mutate_dna <- function(sequence, edits = 1L) {
    chars <- strsplit(sequence, "", fixed = TRUE)[[1L]]
    pos <- sample(seq_along(chars), min(edits, length(chars)))
    for(i in pos) {
      chars[i] <- sample(setdiff(c("A", "C", "G", "T", "N"), chars[i]), 1L)
    }
    paste0(chars, collapse = "")
  }

  stress_dna <- function(n, len = 32L, cluster_count = 1000L) {
    seeds <- unique_random_dna(cluster_count, len, len)
    cluster <- unlist(lapply(seeds, function(x) {
      c(x, mutate_dna(x, 1L), mutate_dna(x, 2L))
    }), use.names = FALSE)
    out <- unique(c(cluster, unique_random_dna(n - length(unique(cluster)), len, len)))
    while(length(out) < n) {
      out <- unique(c(out, unique_random_dna(n - length(out), len, len)))
    }
    out[seq_len(n)]
  }

  startree_cost_matrix <- function(sequences, query, mismatch_cost) {
    chars <- sort(unique(unlist(strsplit(c(sequences, query), "", fixed = TRUE))))
    cost_matrix <- generate_cost_matrix(paste0(chars, collapse = ""), match = 0L, mismatch = mismatch_cost)
    if("N" %in% chars) {
      cost_matrix["N", "N"] <- mismatch_cost
    }
    cost_matrix
  }

  check_fuzz_case <- function(sequences,
                              query,
                              max_distance,
                              compare_class,
                              mismatch_cost = 1L,
                              gap_cost = 1L) {
    sequences <- normalize_unique(sequences)
    query <- toupper(query)
    cost_matrix <- startree_cost_matrix(sequences, query, mismatch_cost)

    st <- star_tree(
      sequences,
      max_distance = max_distance,
      mismatch_cost = mismatch_cost,
      gap_cost = gap_cost,
      nthreads = NTHREADS
    )

    expected_self <- dist_search(
      sequences,
      max_distance = max_distance,
      mode = "levenshtein",
      cost_matrix = cost_matrix,
      gap_cost = gap_cost,
      tree_class = compare_class,
      nthreads = NTHREADS,
      show_progress = FALSE
    )
    actual_self <- canon_pairs(result(st))
    expected_self <- canon_pairs(expected_self)
    if(!identical(actual_self, expected_self)) {
      print(list(
        comparison = "self",
        compare_class = compare_class,
        max_distance = max_distance,
        mismatch_cost = mismatch_cost,
        gap_cost = gap_cost,
        n_sequences = length(sequences)
      ))
      print(head(setdiff(actual_self, expected_self), 20L))
      print(head(setdiff(expected_self, actual_self), 20L))
      dput(sequences)
      stop("StarTree self-search fuzz mismatch")
    }

    expected_query <- dist_search(
      query,
      sequences,
      max_distance = max_distance,
      mode = "levenshtein",
      cost_matrix = cost_matrix,
      gap_cost = gap_cost,
      tree_class = compare_class,
      nthreads = NTHREADS,
      show_progress = FALSE
    )
    actual_query <- sort_result(align_search(st, query))
    expected_query <- sort_result(expected_query)
    if(!identical(actual_query, expected_query)) {
      print(list(
        comparison = "query",
        compare_class = compare_class,
        max_distance = max_distance,
        mismatch_cost = mismatch_cost,
        gap_cost = gap_cost,
        n_sequences = length(sequences),
        n_query = length(query)
      ))
      actual_key <- do.call(paste, c(actual_query, sep = "\t"))
      expected_key <- do.call(paste, c(expected_query, sep = "\t"))
      print(head(setdiff(actual_key, expected_key), 20L))
      print(head(setdiff(expected_key, actual_key), 20L))
      dput(sequences)
      dput(query)
      stop("StarTree query-search fuzz mismatch")
    }
  }

  print("Checking StarTree unit-cost self and query searches")
  target <- c("ACGT", "ACGA", "ACGG", "AAAA", "AAAT", "ACGTA", "acgt", "NNNN", "NNNA")
  tree <- star_tree(target, max_distance = 2L, nthreads = NTHREADS)
  stopifnot(size(tree) == length(normalize_unique(target)))
  stopifnot(identical(to_vector(tree), normalize_unique(target)))
  stopifnot(identical(canon_pairs(result(tree)), canon_pairs(expected_self(target, 2L))))

  query <- c("ACGT", "ACG", "AAAC", "NNNN", NA_character_)
  stopifnot(identical(
    sort_result(align_search(tree, query)),
    expected_search(query, target, 2L)
  ))

  print("Checking StarTree weighted mismatch and gap costs")
  weighted_target <- c("ACGT", "ACG", "AGGT", "TTTT", "ACGTT")
  weighted_tree <- star_tree(
    weighted_target,
    max_distance = 2L,
    mismatch_cost = 2L,
    gap_cost = 1L,
    nthreads = NTHREADS
  )
  stopifnot(identical(
    canon_pairs(result(weighted_tree)),
    canon_pairs(expected_self(weighted_target, 2L, mismatch_cost = 2L, gap_cost = 1L))
  ))
  stopifnot(identical(
    sort_result(align_search(weighted_tree, c("ACGT", "ACGG", "AC"))),
    expected_search(c("ACGT", "ACGG", "AC"), weighted_target, 2L, mismatch_cost = 2L, gap_cost = 1L)
  ))

  print("Checking StarTree zero-distance query search")
  zero_target <- c("TCC", "AAAAAA")
  zero_tree <- star_tree(
    zero_target,
    max_distance = 0L,
    mismatch_cost = 3L,
    gap_cost = 2L,
    nthreads = NTHREADS
  )
  stopifnot(identical(canon_pairs(result(zero_tree)), character()))
  stopifnot(identical(
    sort_result(align_search(zero_tree, "TCC")),
    expected_search("TCC", zero_target, 0L, mismatch_cost = 3L, gap_cost = 2L)
  ))

  print("Checking StarTree leading-pad edit distances")
  pad_target <- c("NAGTNGA", "AGCTGGCTA")
  pad_tree <- star_tree(pad_target, max_distance = 5L, nthreads = NTHREADS)
  stopifnot(identical(
    canon_pairs(result(pad_tree)),
    "AGCTGGCTA\tNAGTNGA\t5"
  ))
  stopifnot(identical(
    sort_result(align_search(pad_tree, "AGCTGGCTA")),
    expected_search("AGCTGGCTA", pad_target, 5L)
  ))
  pad_tree_weighted <- star_tree(
    pad_target,
    max_distance = 11L,
    mismatch_cost = 3L,
    gap_cost = 2L,
    nthreads = NTHREADS
  )
  stopifnot(identical(
    canon_pairs(result(pad_tree_weighted)),
    "AGCTGGCTA\tNAGTNGA\t11"
  ))

  print("Checking StarTree R6 wrapper")
  tree_r6 <- StarTree$new(target, max_distance = 2L, nthreads = NTHREADS)
  stopifnot(tree_r6$size() == size(tree))
  stopifnot(identical(tree_r6$to_vector(), to_vector(tree)))
  stopifnot(identical(canon_pairs(tree_r6$result()), canon_pairs(result(tree))))
  stopifnot(identical(sort_result(tree_r6$search(query)), expected_search(query, target, 2L)))
  stopifnot(is.null(tree_r6$insert))
  stopifnot(is.null(tree_r6$erase))

  print("Checking dist_search StarTree path")
  stopifnot(identical(
    canon_pairs(dist_search(target, max_distance = 2L, tree_class = "StarTree", nthreads = NTHREADS)),
    canon_pairs(expected_self(target, 2L))
  ))
  stopifnot(identical(
    sort_result(dist_search(query, target, max_distance = 2L, tree_class = "star_tree", nthreads = NTHREADS)),
    expected_search(query, target, 2L)
  ))

  print("Checking StarTree against alternating RadixTree and RadixForest")
  set.seed(72219)
  for(i in seq_len(80L)) {
    compare_class <- if(i %% 2L == 1L) "RadixTree" else "RadixForest"
    mismatch_cost <- sample(1:4, 1L)
    gap_cost <- sample(1:4, 1L)
    sequences <- unique_random_dna(sample(20:120, 1L), max_len = sample(8:18, 1L))
    query_fuzz <- random_dna(sample(10:40, 1L), max_len = sample(8:18, 1L))
    check_fuzz_case(
      sequences = sequences,
      query = query_fuzz,
      max_distance = sample(0:(8L * min(mismatch_cost, gap_cost)), 1L),
      compare_class = compare_class,
      mismatch_cost = mismatch_cost,
      gap_cost = gap_cost
    )
  }

  for(i in seq_len(15L)) {
    compare_class <- if(i %% 2L == 1L) "RadixTree" else "RadixForest"
    mismatch_cost <- sample(1:4, 1L)
    gap_cost <- sample(1:4, 1L)
    n_sequences <- sample(10000:20000, 1L)
    sequence_len <- sample(24:48, 1L)
    cluster_count <- sample(500:2000, 1L)
    sequences <- stress_dna(n_sequences, len = sequence_len, cluster_count = cluster_count)
    query_fuzz <- c(
      sample(sequences, 100L),
      vapply(sample(sequences, 100L), mutate_dna, character(1L), edits = 1L),
      random_dna(300L, min_len = sequence_len - 2L, max_len = sequence_len + 2L)
    )
    check_fuzz_case(
      sequences = sequences,
      query = query_fuzz,
      max_distance = sample(1:(4L * min(mismatch_cost, gap_cost)), 1L),
      compare_class = compare_class,
      mismatch_cost = mismatch_cost,
      gap_cost = gap_cost
    )
  }

  print("Checking StarTree restrictions")
  expect_error_message(star_tree(c("ACGT", "AXGT"), max_distance = 1L), "DNA")
  expect_error_message(star_tree(c("ACGT", ""), max_distance = 1L), "empty")
  expect_error_message(star_tree(c("ACGT", NA_character_), max_distance = 1L), "missing")
  expect_error_message(star_tree(c("ACGT", "ACGA"), max_distance = 9L), "max_distance")
  expect_error_message(
    star_tree(c("ACGT", "ACGA"), max_distance = .Machine$integer.max,
              mismatch_cost = .Machine$integer.max, gap_cost = .Machine$integer.max),
    "max_distance"
  )
  expect_error_message(align_search(tree, paste(rep("A", 1024L), collapse = "")), "sequence length")
  expect_error_message(align_search(tree, "ACGT", mode = "hamming"), "global")
  expect_error_message(align_search(tree, "ACGT", max_fraction = 0.1), "max_fraction")
  expect_error_message(align_search(tree, "ACGT", gap_open_cost = 1L), "affine")
  expect_error_message(align_search(tree, "ACGT", lower_triangle = TRUE), "lower_triangle")
  expect_error_message(align_search(tree, "ACGT", match_mode = "best"), "best")
  expect_error_message(dist_search(target, max_fraction = 0.1, tree_class = "StarTree"), "max_fraction")

  print(Sys.time() - runtime)
}
