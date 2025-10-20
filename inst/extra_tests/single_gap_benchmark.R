library(seqtrie)
library(stringdist)
library(Rcpp)
library(dplyr)
library(ggplot2)

NITER <- 3

tic <- function() { .time <<- Sys.time() }
toc <- function() { as.numeric(Sys.time() - .time, units = "secs") }

run_search <- function(query, target, max_distance, nthreads = 8) {
  x <- seqtrie::RadixTree$new()
  x$insert(target)
  max_distance <- seqtrie:::recycle_arg(max_distance, query)
  seqtrie:::RadixTree_search(x$root_pointer, query, max_distance, mode = "anchored", cost_matrix = NULL, gap_cost = 1L, gap_open_cost = NA_integer_, nthreads = nthreads, show_progress = FALSE)
}

run_single_gap_search <- function(query, target, max_distance, nthreads = 8) {
  x <- seqtrie::RadixTree$new()
  x$insert(target)
  max_distance <- seqtrie:::recycle_arg(max_distance, query)
  seqtrie:::RadixTree_single_gap_search(x$root_pointer, query, max_distance, gap_cost = 1L, nthreads, show_progress = FALSE)
}

methods <- list(run_search, run_single_gap_search)
names(methods) <- c("run_search", "run_single_gap_search")

################################################################################

grid <- expand.grid(nseqs = c(100,300,1000,3000,10000), gap_cost = 1, maxdist = 1, iter = 1:NITER, method = names(methods)) %>% sample_n(nrow(.))
grid$time <- rep(0, nrow(grid))
for(i in 1:nrow(grid)) {
  print(grid[i,])
  set.seed(grid$iter[i])
  x <- sample(covid_cdr3, size = grid$nseqs[i])
  tic()
  methods[[grid$method[i]]](x, x, max_distance = grid$maxdist[i])
  grid$time[i] <- toc()
}
grid %>% group_by(nseqs, method, maxdist) %>% summarize(time = mean(time)) %>% as.data.frame %>% print