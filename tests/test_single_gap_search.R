print("Running test_single_gap_search.R")

if(requireNamespace("seqtrie", quietly=TRUE) &&
   requireNamespace("dplyr", quietly=TRUE)
) {

NTHREADS <- 2
cost_mat <- seqtrie::generate_cost_matrix(charset = "ACGT", match = 0L, mismatch = 1L)

dist_matrix_search <- function(query, target, cost_matrix = NULL, gap_cost = NA_integer_, gap_open_cost = NA_integer_, mode = "anchored") {
  results <- seqtrie::dist_matrix(query, target, mode = mode, cost_matrix, gap_cost, gap_open_cost, nthreads=NTHREADS)
  if(mode == "anchored") {
    results <- data.frame(query = rep(query, times=length(target)),
                          target = rep(target, each=length(query)), 
                          distance = as.vector(results),
                          query_size = as.vector(attr(results, "query_size")),
                          target_size = as.vector(attr(results, "target_size")),
                          stringsAsFactors = F)
  } else {
    results <- data.frame(query = rep(query, times=length(target)),
                      target = rep(target, each=length(query)), 
                      distance = as.vector(results),
                      stringsAsFactors = F)
  }
  results <- dplyr::filter(results, is.finite(distance))
  results$distance <- as.integer(results$distance)
  results %>% 
    dplyr::arrange(query, target) %>%
    dplyr::select(-target_size, -query_size)
}

target <- c("ACGT", "ACGGT", "AACGT", "ACGTAA", "GGGG")
query <- c("ACGT", "ACGGT", "ACGTAA", "GGGG")
max_distance <- c(2L, 3L, 1L, 0L)

tree <- RadixTree$new()
tree$insert(target)
result <- tree$single_gap_search(query, max_distance = max_distance, gap_cost = 2L, nthreads = NTHREADS, show_progress = FALSE)
result <- dplyr::arrange(result, query, target)

expected <- lapply(1:length(query), function(i) {
  dist_matrix_search(query[i], target, cost_matrix = cost_mat, gap_cost = 2L, mode = "anchored") %>%
    filter(distance <= max_distance[i])
}) %>% do.call(rbind, .) %>% dplyr::arrange(query, target) %>% as.data.frame

stopifnot(identical(result, expected))

# Check that the single gap condition is strictly enforced
query <- "ABCDEFGHIJKLMNO"
target <- "ABEFGHIJKLMNO" # 2 gaps, missing "CD"
tree <- RadixTree$new()
tree$insert(target)
result <- tree$single_gap_search(query, max_distance = 5L, gap_cost = 2L, nthreads = NTHREADS, show_progress = FALSE)
stopifnot(nrow(result) == 0)

result <- tree$single_gap_search(query, max_distance = 100L, gap_cost = 2L, nthreads = NTHREADS, show_progress = FALSE)
expected_dist <- nchar(target) - 2L # No gap, match AB mismatch EFGHIJKLMNO
stopifnot(identical(result$distance, expected_dist))

tree <- RadixTree$new()
seqs <- sample(covid_cdr3, 1000) %>% gsub("G|T", "A", .) # reduce charset to A and C to get more matches
tree$insert(seqs)
result <- tree$single_gap_search(seqs, max_distance = 8L, gap_cost = 5L, nthreads = NTHREADS, show_progress = FALSE)
expected <- tree$search(seqs, cost_matrix = cost_mat, max_distance = 8L, gap_cost = 5L, mode = "anchored", nthreads = NTHREADS, show_progress = FALSE) %>%
  dplyr::select(-query_size, -target_size)
stopifnot(identical(result, expected))


query <- "GTT"
target <- "GACCC"
tree <- RadixTree$new()
tree$insert(target)
stopifnot(identical(
  tree$single_gap_search(query, max_distance = 3L, gap_cost = 1L, nthreads = NTHREADS, show_progress = FALSE),
  tree$search(query, max_distance = 3L, cost_matrix = cost_mat, gap_cost = 1L, mode = "anchored", nthreads = NTHREADS, show_progress = FALSE) %>% dplyr::select(-query_size, -target_size)
))

}