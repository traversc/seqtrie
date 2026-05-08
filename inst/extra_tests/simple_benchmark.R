suppressPackageStartupMessages({
  library(seqtrie)
  library(dplyr)
})

data(covid_cdr3, package = "seqtrie")

set.seed(314156)

NITER <- 5L
NSEQS <- 30000L
NTHREADS <- 4L
MAX_DISTANCE <- 3L
LINEAR_GAP_COST <- 2L
AFFINE_GAP_COST <- 1L
AFFINE_GAP_OPEN_COST <- 2L

make_cost_matrix <- function(seqs) {
  charset <- sort(unique(unlist(strsplit(seqs, "", fixed = TRUE), use.names = FALSE)))
  seqtrie::generate_cost_matrix(paste0(charset, collapse = ""), match = 0L, mismatch = 1L)
}

run_radixtree_levenshtein_search <- function(query, target, cost_matrix) {
  tree <- seqtrie::RadixTree$new(target)
  tree$search(
    query,
    max_distance = MAX_DISTANCE,
    mode = "levenshtein",
    nthreads = NTHREADS,
    show_progress = FALSE
  )
}

run_radixforest_levenshtein_search <- function(query, target, cost_matrix) {
  forest <- seqtrie::RadixForest$new(target)
  forest$search(
    query,
    max_distance = MAX_DISTANCE,
    mode = "levenshtein",
    nthreads = NTHREADS,
    show_progress = FALSE
  )
}

run_radixtree_hamming_search <- function(query, target, cost_matrix) {
  tree <- seqtrie::RadixTree$new(target)
  tree$search(
    query,
    max_distance = MAX_DISTANCE,
    mode = "hamming",
    nthreads = NTHREADS,
    show_progress = FALSE
  )
}

run_radixforest_hamming_search <- function(query, target, cost_matrix) {
  forest <- seqtrie::RadixForest$new(target)
  forest$search(
    query,
    max_distance = MAX_DISTANCE,
    mode = "hamming",
    nthreads = NTHREADS,
    show_progress = FALSE
  )
}

run_radixtree_anchored_search <- function(query, target, cost_matrix) {
  tree <- seqtrie::RadixTree$new(target)
  tree$search(
    query,
    max_distance = MAX_DISTANCE,
    mode = "anchored",
    nthreads = NTHREADS,
    show_progress = FALSE
  )
}

run_radixtree_global_linear_search <- function(query, target, cost_matrix) {
  tree <- seqtrie::RadixTree$new(target)
  tree$search(
    query,
    max_distance = MAX_DISTANCE,
    mode = "levenshtein",
    cost_matrix = cost_matrix,
    gap_cost = LINEAR_GAP_COST,
    nthreads = NTHREADS,
    show_progress = FALSE
  )
}

run_radixtree_anchored_linear_search <- function(query, target, cost_matrix) {
  tree <- seqtrie::RadixTree$new(target)
  tree$search(
    query,
    max_distance = MAX_DISTANCE,
    mode = "anchored",
    cost_matrix = cost_matrix,
    gap_cost = LINEAR_GAP_COST,
    nthreads = NTHREADS,
    show_progress = FALSE
  )
}

run_radixtree_global_affine_search <- function(query, target, cost_matrix) {
  tree <- seqtrie::RadixTree$new(target)
  tree$search(
    query,
    max_distance = MAX_DISTANCE,
    mode = "levenshtein",
    cost_matrix = cost_matrix,
    gap_cost = AFFINE_GAP_COST,
    gap_open_cost = AFFINE_GAP_OPEN_COST,
    nthreads = NTHREADS,
    show_progress = FALSE
  )
}

run_radixtree_anchored_affine_search <- function(query, target, cost_matrix) {
  tree <- seqtrie::RadixTree$new(target)
  tree$search(
    query,
    max_distance = MAX_DISTANCE,
    mode = "anchored",
    cost_matrix = cost_matrix,
    gap_cost = AFFINE_GAP_COST,
    gap_open_cost = AFFINE_GAP_OPEN_COST,
    nthreads = NTHREADS,
    show_progress = FALSE
  )
}

run_single_gap_search <- function(query, target, cost_matrix) {
  tree <- seqtrie::RadixTree$new(target)
  tree$single_gap_search(
    query,
    max_distance = MAX_DISTANCE,
    gap_cost = 1L,
    nthreads = NTHREADS,
    show_progress = FALSE
  )
}

methods <- list(
  "RadixTree$search levenshtein" = run_radixtree_levenshtein_search,
  "RadixForest$search levenshtein" = run_radixforest_levenshtein_search,
  "RadixTree$search hamming" = run_radixtree_hamming_search,
  "RadixForest$search hamming" = run_radixforest_hamming_search,
  "RadixTree$search anchored" = run_radixtree_anchored_search,
  "RadixTree$search global linear" = run_radixtree_global_linear_search,
  "RadixTree$search anchored linear" = run_radixtree_anchored_linear_search,
  "RadixTree$search global affine" = run_radixtree_global_affine_search,
  "RadixTree$search anchored affine" = run_radixtree_anchored_affine_search,
  "RadixTree$single_gap_search" = run_single_gap_search
)

grid <- expand.grid(
  iter = seq_len(NITER),
  method = names(methods),
  stringsAsFactors = FALSE
)
grid <- grid[sample.int(nrow(grid)), ]
grid$time <- NA_real_
grid$matches <- NA_integer_

for (i in seq_len(nrow(grid))) {
  row <- grid[i, ]

  set.seed(row$iter)
  x <- sample(covid_cdr3, size = NSEQS)
  cost_matrix <- make_cost_matrix(x)

  elapsed <- system.time({
    result <- methods[[row$method]](x, x, cost_matrix)
  })[["elapsed"]]

  grid$time[i] <- elapsed
  grid$matches[i] <- nrow(result)

  rm(x, result, cost_matrix)
  gc(full = TRUE)
}

summary <- grid %>%
  group_by(method) %>%
  summarize(
    mean_time = mean(time),
    median_time = median(time),
    mean_matches = mean(matches),
    median_matches = median(matches),
    .groups = "drop"
  ) %>%
  arrange(method)

summary$mean_time <- round(summary$mean_time, 3L)
summary$median_time <- round(summary$median_time, 3L)
summary$mean_matches <- round(summary$mean_matches, 1L)
summary$median_matches <- round(summary$median_matches, 1L)

old_width <- getOption("width")
options(width = max(old_width, 200L))
print(as.data.frame(summary), row.names = FALSE, right = FALSE)
options(width = old_width)
