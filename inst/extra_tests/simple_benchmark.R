suppressPackageStartupMessages({
  library(seqtrie)
  library(dplyr)
})
data(covid_cdr3)
set.seed(314156)
NITER = 5
NT = 4

grid <- expand.grid(nseqs = c(133034), maxfrac = c(0.035), iter = 1:NITER, method = c("RadixForest")) %>% sample_n(nrow(.))
grid$time <- rep(0, nrow(grid))
for(i in 1:nrow(grid)) {
  x <- sample(covid_cdr3, size = grid$nseqs[i])
  time <- Sys.time()
  r <- seqtrie::dist_search(x, x, max_fraction = grid$maxfrac[i], show_progres = FALSE, tree_class = grid$method[i], nthread=NT)
  grid$time[i] <- as.numeric(Sys.time() - time, units = "secs")
  rm(x, r)
  gc(full=TRUE)
}

cat(mean(grid$time), "\n")
