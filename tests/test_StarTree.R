print("Running test_StarTree.R")
runtime <- Sys.time()

if(requireNamespace("seqtrie", quietly = TRUE)) {
  library(seqtrie)

  IS_LOCAL <- Sys.getenv("IS_LOCAL") != ""
  NTHREADS <- 2L

  # CRAN runs at 2 threads. Locally, vary the thread count per fuzz iteration to
  # surface any thread-count-dependent behaviour in the parallel self-join.
  fuzz_nthreads <- function() if(IS_LOCAL) sample(2:6, 1L) else 2L

  sort_result <- function(x) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    x <- x[, c("query", "target", "distance"), drop = FALSE]
    if(nrow(x) > 0L) {
      x <- x[do.call(order, x), , drop = FALSE]
    }
    rownames(x) <- NULL
    x
  }

  sort_anchored_result <- function(x) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    x <- x[, c("query", "target", "distance", "query_size", "target_size"), drop = FALSE]
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
                              gap_cost = 1L,
                              nthreads = NTHREADS) {
    sequences <- normalize_unique(sequences)
    query <- toupper(query)
    cost_matrix <- startree_cost_matrix(sequences, query, mismatch_cost)

    st <- star_tree(
      sequences,
      max_distance = max_distance,
      mismatch_cost = mismatch_cost,
      gap_cost = gap_cost,
      nthreads = nthreads
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
    actual_query <- sort_result(align_search(st, query, nthreads = nthreads))
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

  check_anchored_fuzz_case <- function(sequences,
                                       query,
                                       max_distance,
                                       compare_class,
                                       mismatch_cost = 1L,
                                       gap_cost = 1L,
                                       nthreads = NTHREADS) {
    sequences <- normalize_unique(sequences)
    query <- toupper(query)
    cost_matrix <- startree_cost_matrix(sequences, query, mismatch_cost)

    st <- star_tree(
      sequences,
      max_distance = max_distance,
      mode = "anchored",
      mismatch_cost = mismatch_cost,
      gap_cost = gap_cost,
      nthreads = nthreads
    )

    expected_self <- dist_search(
      sequences,
      max_distance = max_distance,
      mode = "anchored",
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
        comparison = "anchored self",
        compare_class = compare_class,
        max_distance = max_distance,
        mismatch_cost = mismatch_cost,
        gap_cost = gap_cost,
        n_sequences = length(sequences)
      ))
      print(head(setdiff(actual_self, expected_self), 20L))
      print(head(setdiff(expected_self, actual_self), 20L))
      dput(sequences)
      stop("anchored StarTree self-search fuzz mismatch")
    }

    expected_query <- dist_search(
      query,
      sequences,
      max_distance = max_distance,
      mode = "anchored",
      cost_matrix = cost_matrix,
      gap_cost = gap_cost,
      tree_class = compare_class,
      nthreads = NTHREADS,
      show_progress = FALSE
    )
    actual_query <- sort_anchored_result(align_search(st, query, nthreads = nthreads))
    expected_query <- sort_anchored_result(expected_query)
    if(!identical(actual_query, expected_query)) {
      print(list(
        comparison = "anchored query",
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
      stop("anchored StarTree query-search fuzz mismatch")
    }
  }

  # Hamming oracle: Levenshtein with an N-always-mismatch cost matrix and a gap
  # so costly that no indel can fit within max_distance. For equal-length pairs
  # this is exactly substitution-only Hamming with N mismatching even itself;
  # unequal-length pairs need >= 1 gap and are therefore excluded.
  check_hamming_fuzz_case <- function(sequences,
                                      query,
                                      max_distance,
                                      compare_class,
                                      nthreads = NTHREADS) {
    sequences <- normalize_unique(sequences)
    query <- toupper(query)
    cost_matrix <- startree_cost_matrix(sequences, query, 1L)
    big_gap <- max_distance + 1L

    st <- star_tree(
      sequences,
      max_distance = max_distance,
      mode = "hamming",
      nthreads = nthreads
    )

    expected_self <- dist_search(
      sequences,
      max_distance = max_distance,
      mode = "levenshtein",
      cost_matrix = cost_matrix,
      gap_cost = big_gap,
      tree_class = compare_class,
      nthreads = NTHREADS,
      show_progress = FALSE
    )
    actual_self <- canon_pairs(result(st))
    expected_self <- canon_pairs(expected_self)
    if(!identical(actual_self, expected_self)) {
      print(list(
        comparison = "hamming self",
        compare_class = compare_class,
        max_distance = max_distance,
        n_sequences = length(sequences)
      ))
      print(head(setdiff(actual_self, expected_self), 20L))
      print(head(setdiff(expected_self, actual_self), 20L))
      dput(sequences)
      stop("hamming StarTree self-search fuzz mismatch")
    }

    expected_query <- dist_search(
      query,
      sequences,
      max_distance = max_distance,
      mode = "levenshtein",
      cost_matrix = cost_matrix,
      gap_cost = big_gap,
      tree_class = compare_class,
      nthreads = NTHREADS,
      show_progress = FALSE
    )
    actual_query <- sort_result(align_search(st, query, nthreads = nthreads))
    expected_query <- sort_result(expected_query)
    if(!identical(actual_query, expected_query)) {
      print(list(
        comparison = "hamming query",
        compare_class = compare_class,
        max_distance = max_distance,
        n_sequences = length(sequences),
        n_query = length(query)
      ))
      actual_key <- do.call(paste, c(actual_query, sep = "\t"))
      expected_key <- do.call(paste, c(expected_query, sep = "\t"))
      print(head(setdiff(actual_key, expected_key), 20L))
      print(head(setdiff(expected_key, actual_key), 20L))
      dput(sequences)
      dput(query)
      stop("hamming StarTree query-search fuzz mismatch")
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
  stopifnot(identical(
    sort_result(tree_r6$search(query, nthreads = NTHREADS, show_progress = FALSE)),
    expected_search(query, target, 2L)
  ))
  stopifnot(identical(
    sort_result(tree_r6$align_search(query, nthreads = NTHREADS, show_progress = FALSE)),
    expected_search(query, target, 2L)
  ))
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

  print("Checking StarTree hamming-mode self and query searches")
  hamming_target <- c("ACGT", "ACGA", "ACGG", "TTTT", "ACG", "ACGTT", "NNNN", "ACGN")
  hamming_cost_matrix <- startree_cost_matrix(hamming_target, character(), 1L)
  hamming_tree <- star_tree(
    hamming_target,
    max_distance = 2L,
    mode = "hamming",
    nthreads = NTHREADS
  )
  hamming_expected_self <- dist_search(
    normalize_unique(hamming_target),
    max_distance = 2L,
    mode = "levenshtein",
    cost_matrix = hamming_cost_matrix,
    gap_cost = 3L,
    tree_class = "RadixTree",
    nthreads = NTHREADS
  )
  stopifnot(identical(S7::prop(hamming_tree, "mode"), "hamming"))
  stopifnot(size(hamming_tree) == length(normalize_unique(hamming_target)))
  stopifnot(identical(to_vector(hamming_tree), normalize_unique(hamming_target)))
  stopifnot(identical(canon_pairs(result(hamming_tree)), canon_pairs(hamming_expected_self)))

  hamming_query <- c("ACGT", "ACG", "ACGTT", "NNNN", "TCGA", NA_character_)
  hamming_expected_query <- dist_search(
    hamming_query,
    normalize_unique(hamming_target),
    max_distance = 2L,
    mode = "levenshtein",
    cost_matrix = hamming_cost_matrix,
    gap_cost = 3L,
    tree_class = "RadixTree",
    nthreads = NTHREADS
  )
  stopifnot(identical(
    sort_result(align_search(hamming_tree, hamming_query)),
    sort_result(hamming_expected_query)
  ))

  print("Checking StarTree hamming-mode R6 wrapper and dist_search path")
  hamming_tree_r6 <- StarTree$new(
    hamming_target,
    max_distance = 2L,
    mode = "hamming",
    nthreads = NTHREADS
  )
  stopifnot(identical(hamming_tree_r6$mode, "hamming"))
  stopifnot(identical(canon_pairs(hamming_tree_r6$result()), canon_pairs(result(hamming_tree))))
  stopifnot(identical(
    sort_result(hamming_tree_r6$align_search(hamming_query, nthreads = NTHREADS, show_progress = FALSE)),
    sort_result(hamming_expected_query)
  ))
  stopifnot(identical(
    canon_pairs(dist_search(hamming_target, max_distance = 2L, mode = "hamming",
                            tree_class = "StarTree", nthreads = NTHREADS)),
    canon_pairs(hamming_expected_self)
  ))

  print("Checking StarTree hamming N-vs-N mismatch")
  # N mismatches every base including another N. "NACG" and "NTCG" therefore
  # differ at the shared N position and the A/T position: distance 2, not 1.
  nvn_target <- c("NACG", "NTCG")
  stopifnot(identical(
    canon_pairs(result(star_tree(nvn_target, max_distance = 1L, mode = "hamming", nthreads = NTHREADS))),
    character()
  ))
  stopifnot(identical(
    canon_pairs(result(star_tree(nvn_target, max_distance = 2L, mode = "hamming", nthreads = NTHREADS))),
    "NACG\tNTCG\t2"
  ))

  print("Checking StarTree hamming mode against RadixTree and RadixForest")
  set.seed(40517)
  for(i in seq_len(40L)) {
    compare_class <- if(i %% 2L == 1L) "RadixTree" else "RadixForest"
    sequences <- unique_random_dna(sample(20:120, 1L), max_len = sample(8:18, 1L))
    query_fuzz <- random_dna(sample(10:40, 1L), max_len = sample(8:18, 1L))
    check_hamming_fuzz_case(
      sequences = sequences,
      query = query_fuzz,
      max_distance = sample(0:8, 1L),
      compare_class = compare_class,
      nthreads = fuzz_nthreads()
    )
  }

  print("Checking StarTree anchored-mode unit-cost self and query searches")
  anchored_target <- c("ACGT", "ACG", "ACGG", "AAAA", "AA", "acgt", "NNNN", "NNNA")
  anchored_cost_matrix <- startree_cost_matrix(anchored_target, character(), 1L)
  anchored_tree <- star_tree(
    anchored_target,
    max_distance = 2L,
    mode = "anchored",
    nthreads = NTHREADS
  )
  anchored_expected_self <- dist_search(
    normalize_unique(anchored_target),
    max_distance = 2L,
    mode = "anchored",
    cost_matrix = anchored_cost_matrix,
    gap_cost = 1L,
    tree_class = "RadixTree",
    nthreads = NTHREADS
  )
  stopifnot(size(anchored_tree) == length(normalize_unique(anchored_target)))
  stopifnot(identical(S7::prop(anchored_tree, "mode"), "anchored"))
  stopifnot(identical(sort(to_vector(anchored_tree)), sort(normalize_unique(anchored_target))))
  stopifnot(identical(canon_pairs(result(anchored_tree)), canon_pairs(anchored_expected_self)))

  anchored_query <- c("ACGT", "AC", "AAAAA", "NNNN", NA_character_)
  anchored_expected_query <- dist_search(
    anchored_query,
    normalize_unique(anchored_target),
    max_distance = 2L,
    mode = "anchored",
    cost_matrix = anchored_cost_matrix,
    gap_cost = 1L,
    tree_class = "RadixTree",
    nthreads = NTHREADS
  )
  stopifnot(identical(
    sort_anchored_result(align_search(anchored_tree, anchored_query)),
    sort_anchored_result(anchored_expected_query)
  ))

  print("Checking StarTree anchored-mode weighted mismatch and gap costs")
  anchored_weighted_target <- c("ACGT", "ACG", "AGGT", "TTTT", "ACGTT", "ANNN")
  anchored_weighted_cost_matrix <- startree_cost_matrix(anchored_weighted_target, character(), 2L)
  anchored_weighted_tree <- star_tree(
    anchored_weighted_target,
    max_distance = 4L,
    mode = "anchored",
    mismatch_cost = 2L,
    gap_cost = 2L,
    nthreads = NTHREADS
  )
  anchored_weighted_expected_self <- dist_search(
    normalize_unique(anchored_weighted_target),
    max_distance = 4L,
    mode = "anchored",
    cost_matrix = anchored_weighted_cost_matrix,
    gap_cost = 2L,
    tree_class = "RadixTree",
    nthreads = NTHREADS
  )
  stopifnot(identical(
    canon_pairs(result(anchored_weighted_tree)),
    canon_pairs(anchored_weighted_expected_self)
  ))
  anchored_weighted_expected_query <- dist_search(
    c("ACGT", "ACGG", "AC", "AN"),
    normalize_unique(anchored_weighted_target),
    max_distance = 4L,
    mode = "anchored",
    cost_matrix = anchored_weighted_cost_matrix,
    gap_cost = 2L,
    tree_class = "RadixTree",
    nthreads = NTHREADS
  )
  stopifnot(identical(
    sort_anchored_result(align_search(anchored_weighted_tree, c("ACGT", "ACGG", "AC", "AN"))),
    sort_anchored_result(anchored_weighted_expected_query)
  ))

  print("Checking StarTree anchored-mode R6 wrapper and dist_search path")
  anchored_tree_r6 <- StarTree$new(
    anchored_target,
    max_distance = 2L,
    mode = "anchored",
    nthreads = NTHREADS
  )
  stopifnot(anchored_tree_r6$size() == size(anchored_tree))
  stopifnot(identical(anchored_tree_r6$mode, "anchored"))
  stopifnot(identical(sort(anchored_tree_r6$to_vector()), sort(to_vector(anchored_tree))))
  stopifnot(identical(canon_pairs(anchored_tree_r6$result()), canon_pairs(result(anchored_tree))))
  stopifnot(identical(
    sort_anchored_result(anchored_tree_r6$search(
      anchored_query,
      nthreads = NTHREADS,
      show_progress = FALSE
    )),
    sort_anchored_result(anchored_expected_query)
  ))
  stopifnot(identical(
    sort_anchored_result(anchored_tree_r6$align_search(
      anchored_query,
      nthreads = NTHREADS,
      show_progress = FALSE
    )),
    sort_anchored_result(anchored_expected_query)
  ))
  stopifnot(is.null(anchored_tree_r6$insert))
  stopifnot(is.null(anchored_tree_r6$erase))
  stopifnot(identical(
    canon_pairs(dist_search(anchored_target, max_distance = 2L, mode = "anchored",
                            tree_class = "StarTree", nthreads = NTHREADS)),
    canon_pairs(anchored_expected_self)
  ))
  stopifnot(identical(
    sort_anchored_result(dist_search(anchored_query, anchored_target, max_distance = 2L,
                                     mode = "anchored", tree_class = "star_tree",
                                     nthreads = NTHREADS)),
    sort_anchored_result(anchored_expected_query)
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
      gap_cost = gap_cost,
      nthreads = fuzz_nthreads()
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
      gap_cost = gap_cost,
      nthreads = fuzz_nthreads()
    )
  }

  # RadixForest does not implement anchored mode, so anchored uses RadixTree as
  # the only differential oracle.
  print("Checking StarTree anchored mode against RadixTree")
  set.seed(81227)
  for(i in seq_len(40L)) {
    mismatch_cost <- sample(1:4, 1L)
    gap_cost <- sample(1:4, 1L)
    sequences <- unique_random_dna(sample(20:100, 1L), max_len = sample(8:18, 1L))
    query_fuzz <- random_dna(sample(10:40, 1L), max_len = sample(8:18, 1L))
    check_anchored_fuzz_case(
      sequences = sequences,
      query = query_fuzz,
      max_distance = sample(0:(8L * min(mismatch_cost, gap_cost)), 1L),
      compare_class = "RadixTree",
      mismatch_cost = mismatch_cost,
      gap_cost = gap_cost,
      nthreads = fuzz_nthreads()
    )
  }

  print("Checking StarTree anchored mode at scale")
  # Anchored has no prefilter, so segment length is irrelevant; what needs scale
  # is the LCP-reuse and lower-triangle prune logic against clustered data with
  # length-varied queries (the case anchored mode targets).
  set.seed(60413)
  for(i in seq_len(4L)) {
    mismatch_cost <- sample(1:3, 1L)
    gap_cost <- sample(1:3, 1L)
    sequence_len <- sample(20:30, 1L)
    sequences <- stress_dna(sample(2000:4000, 1L), len = sequence_len,
                            cluster_count = sample(300:800, 1L))
    query_fuzz <- c(
      sample(sequences, 100L),
      vapply(sample(sequences, 100L), mutate_dna, character(1L), edits = 1L),
      random_dna(100L, min_len = sequence_len - 2L, max_len = sequence_len + 2L)
    )
    check_anchored_fuzz_case(
      sequences = sequences,
      query = query_fuzz,
      max_distance = sample(1:(3L * min(mismatch_cost, gap_cost)), 1L),
      compare_class = "RadixTree",
      mismatch_cost = mismatch_cost,
      gap_cost = gap_cost,
      nthreads = fuzz_nthreads()
    )
  }

  print("Checking StarTree prefilter with large segments (K > 14)")
  # Long sequences with small tau give segments far longer than the old K=14
  # bitmap cap (and many beyond 32), exercising the uncapped wyhash prefilter on
  # both the unit-cost and weighted paths, plus the hamming trie.
  set.seed(99001)
  for(i in seq_len(10L)) {
    compare_class <- if(i %% 2L == 1L) "RadixTree" else "RadixForest"
    seq_length <- sample(64:90, 1L)
    n_sequences <- sample(1500:2500, 1L)
    cluster_count <- sample(150:400, 1L)
    sequences <- stress_dna(n_sequences, len = seq_length, cluster_count = cluster_count)
    query_fuzz <- c(
      sample(sequences, 50L),
      vapply(sample(sequences, 50L), mutate_dna, character(1L), edits = 2L),
      random_dna(80L, min_len = seq_length, max_len = seq_length)
    )
    max_distance <- sample(1:3, 1L)
    nthreads <- fuzz_nthreads()
    check_fuzz_case(sequences, query_fuzz, max_distance, compare_class,
                    nthreads = nthreads)
    check_fuzz_case(sequences, query_fuzz, max_distance, compare_class,
                    mismatch_cost = 2L, gap_cost = 3L, nthreads = nthreads)
    check_hamming_fuzz_case(sequences, query_fuzz, max_distance, compare_class,
                            nthreads = nthreads)
  }

  print("Checking StarTree at the tau = 8 boundary")
  # tau = max_distance / min(mismatch_cost, gap_cost); 8 is the maximum allowed.
  # Pin it deterministically for the global and anchored paths.
  set.seed(8088)
  boundary_seqs <- unique_random_dna(80L, min_len = 16L, max_len = 20L)
  boundary_query <- random_dna(40L, min_len = 16L, max_len = 20L)
  check_fuzz_case(boundary_seqs, boundary_query, max_distance = 8L, "RadixTree")
  check_anchored_fuzz_case(boundary_seqs, boundary_query, max_distance = 8L, "RadixTree")

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
  expect_error_message(align_search(tree, "ACGT", mode = "hamming"), "mode")
  expect_error_message(align_search(tree, "ACGT", max_distance = 2L), "max_distance")
  expect_error_message(align_search(tree, "ACGT", max_fraction = 0.1), "max_fraction")
  expect_error_message(align_search(tree, "ACGT", gap_cost = 1L), "gap_cost")
  expect_error_message(align_search(tree, "ACGT", gap_open_cost = 1L), "gap_open_cost")
  expect_error_message(align_search(tree, "ACGT", lower_triangle = TRUE), "lower_triangle")
  expect_error_message(align_search(tree, "ACGT", match_mode = "best"), "match_mode")
  expect_error_message(dist_search(target, max_fraction = 0.1, tree_class = "StarTree"), "max_fraction")

  print("Checking StarTree anchored-mode restrictions")
  expect_error_message(star_tree(c("ACGT", "AXGT"), max_distance = 1L, mode = "anchored"), "DNA")
  expect_error_message(star_tree(c("ACGT", ""), max_distance = 1L, mode = "anchored"), "empty")
  expect_error_message(star_tree(c("ACGT", NA_character_), max_distance = 1L, mode = "anchored"), "missing")
  expect_error_message(star_tree(c("ACGT", "ACGA"), max_distance = 9L, mode = "anchored"), "max_distance")
  expect_error_message(
    star_tree(c("ACGT", "ACGA"), max_distance = .Machine$integer.max,
              mode = "anchored", mismatch_cost = .Machine$integer.max,
              gap_cost = .Machine$integer.max),
    "max_distance"
  )
  expect_error_message(align_search(anchored_tree, paste(rep("A", 1024L), collapse = "")), "sequence length")
  expect_error_message(align_search(anchored_tree, "ACGT", mode = "global"), "mode")
  expect_error_message(align_search(anchored_tree, "ACGT", max_distance = 2L), "max_distance")
  expect_error_message(align_search(anchored_tree, "ACGT", max_fraction = 0.1), "max_fraction")
  expect_error_message(align_search(anchored_tree, "ACGT", gap_cost = 1L), "gap_cost")
  expect_error_message(align_search(anchored_tree, "ACGT", gap_open_cost = 1L), "gap_open_cost")
  expect_error_message(align_search(anchored_tree, "ACGT", lower_triangle = TRUE), "lower_triangle")
  expect_error_message(align_search(anchored_tree, "ACGT", match_mode = "best"), "match_mode")
  expect_error_message(dist_search(anchored_target, max_fraction = 0.1,
                                   mode = "anchored", tree_class = "StarTree"), "max_fraction")

  print(Sys.time() - runtime)
}
