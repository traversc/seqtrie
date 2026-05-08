suppressPackageStartupMessages({
  library(seqtrie)
})

data(covid_cdr3, package = "seqtrie")

bench_param <- function(name, default) {
  as.integer(Sys.getenv(paste0("SEQTRIE_BENCH_", name), default))
}

NITER <- bench_param("NITER", 5L)
NSEQS <- bench_param("NSEQS", 10000L)
NTHREADS <- bench_param("NTHREADS", 4L)
MAX_DISTANCE <- bench_param("MAX_DISTANCE", 3L)

run_radixtree_levenshtein_search <- function(query, target) {
  tree <- seqtrie::RadixTree$new(target)
  tree$search(
    query,
    max_distance = MAX_DISTANCE,
    mode = "levenshtein",
    nthreads = NTHREADS,
    show_progress = FALSE
  )
}

run_radixforest_levenshtein_search <- function(query, target) {
  forest <- seqtrie::RadixForest$new(target)
  forest$search(
    query,
    max_distance = MAX_DISTANCE,
    mode = "levenshtein",
    nthreads = NTHREADS,
    show_progress = FALSE
  )
}

methods <- list(
  "RadixTree$search levenshtein" = run_radixtree_levenshtein_search,
  "RadixForest$search levenshtein" = run_radixforest_levenshtein_search
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

  elapsed <- system.time({
    result <- methods[[row$method]](x, x)
  })[["elapsed"]]

  grid$time[i] <- elapsed
  grid$matches[i] <- nrow(result)

  rm(x, result)
  gc(full = TRUE)
}

summary <- do.call(rbind, lapply(unique(grid$method), function(method) {
  rows <- grid[grid$method == method, ]
  data.frame(
    method = method,
    mean_time = round(mean(rows$time), 3L),
    median_time = round(median(rows$time), 3L),
    mean_matches = round(mean(rows$matches), 1L),
    median_matches = round(median(rows$matches), 1L),
    stringsAsFactors = FALSE
  )
}))
summary <- summary[order(summary$method), ]

old_width <- getOption("width")
options(width = max(old_width, 200L))
print(summary, row.names = FALSE, right = FALSE)
options(width = old_width)
