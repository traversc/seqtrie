print("Running test_search_hook.R")

if(requireNamespace("seqtrie", quietly = TRUE)) {
  library(seqtrie)

  expect_head_equal <- function(res) {
    full <- res$full
    early <- res$early
    stopifnot(nrow(early) == 1L)
    stopifnot(identical(early$target, full$target[1]))
    stopifnot(identical(early$distance, full$distance[1]))
    stopifnot(identical(early$query, full$query[1]))
  }

  res <- seqtrie:::test_search_hook()
  lapply(res, expect_head_equal)
}
