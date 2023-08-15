# This test file tests the RadixForest class in R/r6_class.r
# 1) That insertion and deletion produce the correct results with random strings
# 2) That search (hamming, levenshtein and anchored) produce the same results as the internal `dist_matrix` and `dist_pairwise` functions

library(seqtrie)
library(stringdist)
library(stringfish)
library(dplyr)

# Use 2 threads on github actions and CRAN, 4 threads locally
IS_LOCAL  <- Sys.getenv("IS_LOCAL") != ""
NTHREADS  <- ifelse(IS_LOCAL, 4, 2)
NITER     <- ifelse(IS_LOCAL, 4, 2)
NSEQS     <- 10000 # must be larger than 1000
MAXSEQLEN <- 200
MAXDIST   <- MAXSEQLEN * 0.05
MAXFRAC   <- 0.05
CHARSET   <- "ACGT"

forest_equal <- function(x, y) {
  xs <- gsub("[0-9]+", "#", x$to_string())
  ys <- gsub("[0-9]+", "#", y$to_string())
  if(!setequal(xs, ys)) return(FALSE)
  if(x$size() != y$size()) return(FALSE)
  return(TRUE)
}

random_strings <- function(N, charset = "abcdefghijklmnopqrstuvwxyz") {
  len <- sample(0:MAXSEQLEN, N, replace=T)
  result <- lapply(0:MAXSEQLEN, function(x) {
    nx <- sum(len == x)
    stringfish::random_strings(nx, x, charset = charset, vector_mode = "normal")
  })
  sample(unlist(result))
}

mutate_strings <- function(x, prob = 0.025, indel_prob = 0.025, charset = "abcdefghijklmnopqrstuvwxyz") {
  charset <- unlist(strsplit(charset, ""))
  xsplit <- strsplit(x, "")
  sapply(xsplit, function(a) {
    r <- runif(length(a)) < prob
    a[r] <- sample(charset, sum(r), replace=T)
    ins <- runif(length(a)) < indel_prob
    a[ins] <- paste0(sample(charset, sum(ins), replace=T), sample(charset, sum(ins), replace=T))
    del <- runif(length(a)) < indel_prob
    a[del] <- ""
    paste0(a, collapse = "")
  })
}

sd_search <- function(query, target, method = "lv") {
  results <- stringdist::stringdistmatrix(query, target, method = method, nthread=4)
  results <- data.frame(query = rep(query, times=length(target)), 
                        target = rep(target, each=length(query)), 
                        distance = as.vector(results), stringsAsFactors = F)
  results <- dplyr::filter(results, is.finite(distance))
  results$distance <- as.integer(results$distance)
  dplyr::arrange(results, query, target)
}

tt <- "RadixForest"
for(. in 1:NITER) {
  print(paste0("Checking correct insert/erase methods for ", tt))
  local({
    x <- RadixForest$new()
    y <- RadixForest$new()
    ins <- c(random_strings(NSEQS, CHARSET),"")
    era <- c(sample(c(sample(ins, NSEQS/10), random_strings(NSEQS/10, CHARSET))),"")
    x$insert(ins)
    stopifnot(x$validate())
    stopifnot(x$size() == n_distinct(ins))
    x$erase(era)
    stopifnot(x$validate())
    stopifnot(x$size() == n_distinct(ins[!ins %in% era]))
    y$insert(ins[!ins %in% era])
    stopifnot(y$validate())
    stopifnot(forest_equal(x, y))
  })

  print(paste0('Checking find for ', tt))
  local({
    x <- RadixForest$new()
    ins <- c(random_strings(NSEQS, CHARSET),"")
    era <- c(sample(c(sample(ins, NSEQS/10), random_strings(NSEQS/10, CHARSET))))
    fin <- c(sample(c(sample(ins, NSEQS/10), random_strings(NSEQS/10, CHARSET))),"")
    expected <- fin %in% setdiff(ins, era)
    x$insert(ins)
    stopifnot(x$validate())
    x$erase(era)
    stopifnot(x$validate())
    results <- x$find(fin)
    stopifnot(identical(results, expected))
  })

  print(paste0('Checking prefix_search for ', tt))
  local({
    x <- RadixForest$new()
    ins <- c(random_strings(NSEQS, CHARSET),"")
    era <- c(sample(c(sample(ins, NSEQS/10), random_strings(NSEQS/10, CHARSET))))
    fin <- c(sample(c(sample(ins, NSEQS/1000), random_strings(NSEQS/1000, CHARSET))),"") %>% substr(1,5)
    fin <- c(fin, paste0(fin, substr(fin,1,1)))
    ins2 <- setdiff(ins, era)
    expected <- lapply(fin, function(f) {
      ex <- grep(paste0("^", f), ins2, value=T)
      if(length(ex) == 0) return(NULL)
      data.frame(query = f, target = ex, stringsAsFactors = F)
    })
    expected <- do.call(rbind, expected)
    expected <- dplyr::arrange(expected, query, target)
    x$insert(ins)
    stopifnot(x$validate())
    x$erase(era)
    stopifnot(x$validate())
    results <- x$prefix_search(fin) %>% dplyr::arrange(query, target)
    stopifnot(identical(results, expected))
  })

  print(paste0("Checking multithreaded hamming search correctness for ", tt))
  local({
    x <- RadixForest$new()
    target <- c(random_strings(NSEQS, CHARSET),"") %>% unique
    query <- sample(c(sample(target, NSEQS/1000), random_strings(NSEQS/1000, CHARSET)))
    query <- c(mutate_strings(query, indel_prob=0, charset = CHARSET), "") %>% unique
    x$insert(target)
    stopifnot(x$validate())
    results_dist <- x$search(query, max_distance = MAXDIST, mode = "hamming", nthreads=NTHREADS, show_progress=TRUE) %>% dplyr::arrange(query, target)
    results_frac <- x$search(query, max_fraction = MAXFRAC, mode = "hamming", nthreads=NTHREADS, show_progress=TRUE) %>% dplyr::arrange(query, target)
    sd_results <- sd_search(query, target, method = "hamming")
    sd_dist <- dplyr::filter(sd_results, distance <= MAXDIST)
    sd_frac <- dplyr::filter(sd_results, distance <= nchar(query) * MAXFRAC)
    stopifnot(identical(results_dist, sd_dist))
    stopifnot(identical(results_frac, sd_frac))
  })

    print(paste0("Checking multithreaded levenshtein search correctness for ", tt))
    local({
      x <- RadixForest$new()
      target <- c(random_strings(NSEQS, CHARSET),"") %>% unique
      query <- sample(c(sample(target, NSEQS/1000), random_strings(NSEQS/1000, CHARSET)))
      query <- c(mutate_strings(query, charset = CHARSET), "") %>% unique
      x$insert(target)
      stopifnot(x$validate())
      results_dist <- x$search(query, max_distance = MAXDIST, mode = "levenshtein", nthreads=NTHREADS, show_progress=TRUE) %>% dplyr::arrange(query, target)
      results_frac <- x$search(query, max_fraction = MAXFRAC, mode = "levenshtein", nthreads=NTHREADS, show_progress=TRUE) %>% dplyr::arrange(query, target)
      sd_results <- sd_search(query, target, method = "lv")
      sd_dist <- dplyr::filter(sd_results, distance <= MAXDIST)
      sd_frac <- dplyr::filter(sd_results, distance <= nchar(query) * MAXFRAC)
      stopifnot(identical(results_dist, sd_dist))
      stopifnot(identical(results_frac, sd_frac))
    })  
}
