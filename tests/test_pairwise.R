# This test file tests the `dist_matrix` and `dist_pairwise` functions
# These two functions are simple dynamic programming algorithms for computing pairwise distances and are themselves used to validate
# the RadixTree imeplementation (see test_radix_tree.R)

library(seqtrie)
library(stringi)
library(stringdist)
library(Biostrings)
library(dplyr)

# Use 2 threads on github actions and CRAN, 4 threads locally
IS_LOCAL  <- Sys.getenv("IS_LOCAL") != ""
NTHREADS  <- ifelse(IS_LOCAL, 4, 2)
NITER     <- ifelse(IS_LOCAL, 4, 1)
NSEQS     <- 10000
MAXSEQLEN <- 200
CHARSET   <- "ACGT"

random_strings <- function(N, charset = "abcdefghijklmnopqrstuvwxyz") {
  charset_stri <- paste0("[", charset, "]")
  len <- sample(0:MAXSEQLEN, N, replace=TRUE)
  result <- lapply(0:MAXSEQLEN, function(x) {
    nx <- sum(len == x)
    if(nx == 0) return(character())
    stringi::stri_rand_strings(nx, x, pattern = charset_stri)
  })
  sample(unlist(result))
}

mutate_strings <- function(x, prob = 0.025, indel_prob = 0.025, charset = "abcdefghijklmnopqrstuvwxyz") {
  charset <- unlist(strsplit(charset, ""))
  xsplit <- strsplit(x, "")
  sapply(xsplit, function(a) {
    r <- runif(length(a)) < prob
    a[r] <- sample(charset, sum(r), replace=TRUE)
    ins <- runif(length(a)) < indel_prob
    a[ins] <- paste0(sample(charset, sum(ins), replace=TRUE), sample(charset, sum(ins), replace=TRUE))
    del <- runif(length(a)) < indel_prob
    a[del] <- ""
    paste0(a, collapse = "")
  })
}

# Biostrings notes:
# subject (target) must be of length 1 or equal to pattern (query)
# To get a distance matrix, iterate over target and perform a column bind
# special_zero_case -- if both query and target are empty, Biostrings fails with an error
pairwiseAlignmentFix <- function(pattern, subject, ...) {
    results <- rep(0, length(subject))
    special_zero_case <- nchar(pattern) == 0 & nchar(subject) == 0
    if(all(special_zero_case)) {
      results
    } else {
      results[!special_zero_case] <- Biostrings::pairwiseAlignment(pattern=pattern[!special_zero_case], subject=subject[!special_zero_case], ...)
      results
    }
}

biostrings_matrix_global <- function(query, target, cost_matrix, gap_cost, gap_open_cost = 0) {
  substitutionMatrix <- -cost_matrix
  lapply(query, function(x) {
    query2 <- rep(x, length(target))
    -pairwiseAlignmentFix(pattern=query2, subject=target, substitutionMatrix = substitutionMatrix, gapOpening=gap_open_cost, gapExtension=gap_cost, scoreOnly=TRUE, type="global")
  }) %>% do.call(rbind, .)
}

biostrings_pairwise_global <- function(query, target, cost_matrix, gap_cost, gap_open_cost = 0) {
  substitutionMatrix <- -cost_matrix
  -pairwiseAlignment(pattern=query, subject=target, substitutionMatrix = substitutionMatrix,gapOpening=gap_open_cost, gapExtension=gap_cost, scoreOnly=TRUE, type="global")
}

biostrings_matrix_anchored <- function(query, target, query_size, target_size, cost_matrix, gap_cost, gap_open_cost = 0) {
  substitutionMatrix <- -cost_matrix
  lapply(seq_along(query), function(i) {
    query2 <- substring(query[i], 1, query_size[i,,drop=TRUE])
    target2 <- substring(target, 1, target_size[i,,drop=TRUE])
    -pairwiseAlignmentFix(pattern=query2, subject=target2, substitutionMatrix = substitutionMatrix, gapOpening=gap_open_cost, gapExtension=gap_cost, scoreOnly=TRUE, type="global")
  }) %>% do.call(rbind, .)
}

biostrings_pairwise_anchored <- function(query, target, query_size, target_size, cost_matrix, gap_cost, gap_open_cost = 0) {
  substitutionMatrix <- -cost_matrix
  query2 <- substring(query, 1, query_size)
  target2 <- substring(target, 1, target_size)
  -pairwiseAlignmentFix(pattern=query2, subject=target2, substitutionMatrix = substitutionMatrix, gapOpening=gap_open_cost, gapExtension=gap_cost, scoreOnly=TRUE, type="global")
}

for(. in 1:NITER) {

    print("Checking hamming search correctness")
    local({
      # Note: seqtrie returns `NA_integer_` for hamming distance when the lengths are different
      # whereas stringdist returns `Inf`
      # This is why we need to replace `NA_integer_` with `Inf` when comparing results

      target <- c(random_strings(NSEQS, CHARSET),"") %>% unique
      query <- sample(c(sample(target, NSEQS/1000), random_strings(NSEQS/1000, CHARSET)))
      query <- c(mutate_strings(query, indel_prob=0, charset = CHARSET), "") %>% unique

      # Check matrix results
      results_seqtrie <- dist_matrix(query, target, mode = "hamming", nthreads=NTHREADS)
      results_seqtrie[is.na(results_seqtrie)] <- Inf
      results_stringdist <- stringdist::stringdistmatrix(query, target, method = "hamming", nthread=NTHREADS)
      stopifnot(all(results_seqtrie == results_stringdist))

      # Check pairwise results
      query_pairwise <- mutate_strings(target, prob=0.025, indel_prob=0.05, charset = CHARSET)
      results_seqtrie <- dist_pairwise(query_pairwise, target, mode = "hamming", nthreads=NTHREADS)
      results_seqtrie[is.na(results_seqtrie)] <- Inf
      results_stringdist <- stringdist::stringdist(query_pairwise, target, method = "hamming", nthread=NTHREADS)
      stopifnot(all(results_seqtrie == results_stringdist))
    })

    print("Checking levenshtein search correctness")
    local({
      target <- c(random_strings(NSEQS, CHARSET),"") %>% unique
      query <- sample(c(sample(target, NSEQS/1000), random_strings(NSEQS/1000, CHARSET)))
      query <- c(mutate_strings(query, indel_prob=0, charset = CHARSET), "") %>% unique

      # Check matrix results
      results_seqtrie <- dist_matrix(query, target, mode = "levenshtein", nthreads=NTHREADS)
      results_stringdist <- stringdist::stringdistmatrix(query, target, method = "lv", nthread=NTHREADS)
      stopifnot(all(results_seqtrie == results_stringdist))

      # Check pairwise results
      query_pairwise <- mutate_strings(target, prob=0.025, indel_prob=0.05, charset = CHARSET)
      results_seqtrie <- dist_pairwise(query_pairwise, target, mode = "levenshtein", nthreads=NTHREADS)
      results_stringdist <- stringdist::stringdist(query_pairwise, target, method = "lv", nthread=NTHREADS)
      stopifnot(all(results_seqtrie == results_stringdist))
    })

    print("Checking anchored search correctness")
    local({
      # There is no anchored search in stringdist (or elsewhere). To get the same results, we substring the query and target sequences
      # By the results of the seqtrie anchored search and then compare the results

      target <- c(random_strings(NSEQS, CHARSET),"") %>% unique
      query <- sample(c(sample(target, NSEQS/1000), random_strings(NSEQS/1000, CHARSET)))
      query <- c(mutate_strings(query, indel_prob=0, charset = CHARSET), "") %>% unique

      # Check matrix results
      results_seqtrie <- dist_matrix(query, target, mode = "anchored", nthreads=NTHREADS)
      query_size <- attr(results_seqtrie, "query_size")
      target_size <- attr(results_seqtrie, "target_size")
      results_stringdist <- lapply(seq_along(query), function(i) {
        query_size2 <- query_size[i,,drop=TRUE]
        target_size2 <- target_size[i,,drop=TRUE]
        query2 <- substring(query[i], 1, query_size2) # query[i] is recycled
        target2 <- substring(target, 1, target_size2)
        stringdist::stringdist(query2, target2, method = "lv", nthread=NTHREADS)
      }) %>% do.call(rbind, .)
      stopifnot(all(results_seqtrie == results_stringdist))

      # Check pairwise results
      query_pairwise <- mutate_strings(target, prob=0.025, indel_prob=0.05, charset = CHARSET)
      results_seqtrie <- dist_pairwise(query_pairwise, target, mode = "anchored", nthreads=NTHREADS)
      query_size <- attr(results_seqtrie, "query_size")
      target_size <- attr(results_seqtrie, "target_size")
      query2 <- substring(query_pairwise, 1, query_size)
      target2 <- substring(target, 1, target_size)
      results_stringdist <- stringdist::stringdist(query2, target2, method = "lv", nthread=NTHREADS)
      stopifnot(all(results_seqtrie == results_stringdist))
    })

    print("Checking global search with linear gap for correctness")
    local({
      target <- c(random_strings(NSEQS, CHARSET),"") %>% unique
      query <- sample(c(sample(target, NSEQS/1000), random_strings(NSEQS/1000, CHARSET)))
      query <- c(mutate_strings(query, indel_prob=0, charset = CHARSET), "") %>% unique

      # Check matrix results
      cost_matrix <- matrix(sample(1:3, size = nchar(CHARSET)^2, replace=TRUE), nrow=nchar(CHARSET))
      diag(cost_matrix) <- 0
      colnames(cost_matrix) <- rownames(cost_matrix) <- strsplit(CHARSET, "")[[1]]
      gap_cost <- sample(1:3, size = 1)
      results_seqtrie <- dist_matrix(query, target, mode = "levenshtein", cost_matrix = cost_matrix, gap_cost = gap_cost, nthreads=NTHREADS)
      results_biostrings <- biostrings_matrix_global(query, target, cost_matrix = cost_matrix, gap_cost = gap_cost)
      stopifnot(all(results_seqtrie == results_biostrings))

      # Check pairwise results
      query_pairwise <- mutate_strings(target, prob=0.025, indel_prob=0.05, charset = CHARSET)
      results_seqtrie <- dist_pairwise(query_pairwise, target, mode = "levenshtein", cost_matrix = cost_matrix, gap_cost = gap_cost, nthreads=NTHREADS)
      results_biostrings <- biostrings_pairwise_global(query_pairwise, target, cost_matrix = cost_matrix, gap_cost = gap_cost)
      stopifnot(all(results_seqtrie == results_biostrings))
    })

    print("Checking anchored search with linear gap for correctness")
    local({
      target <- c(random_strings(NSEQS, CHARSET),"") %>% unique
      query <- sample(c(sample(target, NSEQS/1000), random_strings(NSEQS/1000, CHARSET)))
      query <- c(mutate_strings(query, indel_prob=0, charset = CHARSET), "") %>% unique

      # Check matrix results
      cost_matrix <- matrix(sample(1:3, size = nchar(CHARSET)^2, replace=TRUE), nrow=nchar(CHARSET))
      diag(cost_matrix) <- 0
      colnames(cost_matrix) <- rownames(cost_matrix) <- strsplit(CHARSET, "")[[1]]
      gap_cost <- sample(1:3, size = 1)
      results_seqtrie <- dist_matrix(query, target, mode = "anchored", cost_matrix = cost_matrix, gap_cost = gap_cost, nthreads=NTHREADS)
      query_size <- attr(results_seqtrie, "query_size")
      target_size <- attr(results_seqtrie, "target_size")
      results_biostrings <- biostrings_matrix_anchored(query, target, query_size, target_size, cost_matrix = cost_matrix, gap_cost = gap_cost)
      stopifnot(all(results_seqtrie == results_biostrings))

      # Check pairwise results
      query_pairwise <- mutate_strings(target, prob=0.025, indel_prob=0.05, charset = CHARSET)
      results_seqtrie <- dist_pairwise(query_pairwise, target, mode = "anchored", cost_matrix = cost_matrix, gap_cost = gap_cost, nthreads=NTHREADS)
      query_size <- attr(results_seqtrie, "query_size")
      target_size <- attr(results_seqtrie, "target_size")
      results_biostrings <- biostrings_pairwise_anchored(query_pairwise, target, query_size, target_size, cost_matrix = cost_matrix, gap_cost = gap_cost)
      stopifnot(all(results_seqtrie == results_biostrings))
    })


    print("Checking global search with affine gap for correctness")
    local({
      target <- c(random_strings(NSEQS, CHARSET),"") %>% unique
      query <- sample(c(sample(target, NSEQS/1000), random_strings(NSEQS/1000, CHARSET)))
      query <- c(mutate_strings(query, indel_prob=0, charset = CHARSET), "") %>% unique

      # Check matrix results
      cost_matrix <- matrix(sample(1:3, size = nchar(CHARSET)^2, replace=TRUE), nrow=nchar(CHARSET))
      diag(cost_matrix) <- 0
      colnames(cost_matrix) <- rownames(cost_matrix) <- strsplit(CHARSET, "")[[1]]
      gap_cost <- sample(1:3, size = 1)
      gap_open_cost <- sample(1:3, size = 1)
      results_seqtrie <- dist_matrix(query, target, mode = "levenshtein", cost_matrix = cost_matrix, gap_cost = gap_cost, gap_open_cos=gap_open_cost, nthreads=NTHREADS)
      results_biostrings <- biostrings_matrix_global(query, target, cost_matrix = cost_matrix, gap_cost = gap_cost, gap_open_cost=gap_open_cost)
      stopifnot(all(results_seqtrie == results_biostrings))

      # Check pairwise results
      query_pairwise <- mutate_strings(target, prob=0.025, indel_prob=0.05, charset = CHARSET)
      results_seqtrie <- dist_pairwise(query_pairwise, target, mode = "levenshtein", cost_matrix = cost_matrix, gap_cost = gap_cost, gap_open_cost=gap_open_cost, nthreads=NTHREADS)
      results_biostrings <- biostrings_pairwise_global(query_pairwise, target, cost_matrix = cost_matrix, gap_cost = gap_cost, gap_open_cost=gap_open_cost)
      stopifnot(all(results_seqtrie == results_biostrings))
    })

    print("Checking anchored search with affine gap for correctness")
    local({
      target <- c(random_strings(NSEQS, CHARSET),"") %>% unique
      query <- sample(c(sample(target, NSEQS/1000), random_strings(NSEQS/1000, CHARSET)))
      query <- c(mutate_strings(query, indel_prob=0, charset = CHARSET), "") %>% unique

      # Check matrix results
      cost_matrix <- matrix(sample(1:3, size = nchar(CHARSET)^2, replace=TRUE), nrow=nchar(CHARSET))
      diag(cost_matrix) <- 0
      colnames(cost_matrix) <- rownames(cost_matrix) <- strsplit(CHARSET, "")[[1]]
      gap_cost <- sample(1:3, size = 1)
      gap_open_cost <- sample(1:3, size = 1)
      results_seqtrie <- dist_matrix(query, target, mode = "anchored", cost_matrix = cost_matrix, gap_cost = gap_cost, gap_open_cost=gap_open_cost, nthreads=NTHREADS)
      query_size <- attr(results_seqtrie, "query_size")
      target_size <- attr(results_seqtrie, "target_size")
      results_biostrings <- biostrings_matrix_anchored(query, target, query_size, target_size, cost_matrix = cost_matrix, gap_cost = gap_cost, gap_open_cost=gap_open_cost)
      stopifnot(all(results_seqtrie == results_biostrings))

      # Check pairwise results
      query_pairwise <- mutate_strings(target, prob=0.025, indel_prob=0.05, charset = CHARSET)
      results_seqtrie <- dist_pairwise(query_pairwise, target, mode = "anchored", cost_matrix = cost_matrix, gap_cost = gap_cost, gap_open_cost=gap_open_cost, nthreads=NTHREADS)
      query_size <- attr(results_seqtrie, "query_size")
      target_size <- attr(results_seqtrie, "target_size")
      results_biostrings <- biostrings_pairwise_anchored(query_pairwise, target, query_size, target_size, cost_matrix = cost_matrix, gap_cost = gap_cost, gap_open_cost=gap_open_cost)
      stopifnot(all(results_seqtrie == results_biostrings))
    })
}
