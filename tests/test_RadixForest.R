# This test file tests the RadixForest class in R/r6_class.r
# 1) That insertion and deletion produce the correct results with random strings
# 2) That search (hamming, levenshtein and anchored) produce the same results as the internal `dist_matrix` and `dist_pairwise` functions

if(requireNamespace("seqtrie", quietly=TRUE)) {


library(seqtrie)

# Use 2 threads on github actions and CRAN, 4 threads locally
IS_LOCAL  <- Sys.getenv("IS_LOCAL") != ""
NTHREADS  <- ifelse(IS_LOCAL, 4, 2)
NITER     <- ifelse(IS_LOCAL, 4, 1)
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
  charset <- unlist(strsplit(charset, "", fixed = TRUE))
  len <- sample(0:MAXSEQLEN, N, replace=TRUE)
  vapply(len, function(n) {
    paste0(sample(charset, n, replace = TRUE), collapse = "")
  }, character(1))
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

arrange_result <- function(results) {
  if(is.null(results)) {
    return(data.frame())
  }
  results <- as.data.frame(results, stringsAsFactors = FALSE)
  if(nrow(results) > 0L) {
    results <- results[order(results$query, results$target), , drop = FALSE]
  }
  rownames(results) <- NULL
  results
}

dist_matrix_search <- function(query, target, mode = "levenshtein") {
  results <- seqtrie::dist_matrix(query, target, mode = mode, nthreads=NTHREADS)
  results <- data.frame(query = rep(query, times=length(target)), 
                        target = rep(target, each=length(query)), 
                        distance = as.vector(results), stringsAsFactors = FALSE)
  results <- results[is.finite(results$distance), , drop = FALSE]
  results$distance <- as.integer(results$distance)
  arrange_result(results)
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
    stopifnot(x$size() == length(unique(ins)))
    x$erase(era)
    stopifnot(x$validate())
    stopifnot(x$size() == length(unique(ins[!ins %in% era])))
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
    fin <- substr(c(sample(c(sample(ins, NSEQS/1000), random_strings(NSEQS/1000, CHARSET))), ""), 1, 5)
    fin <- c(fin, paste0(fin, substr(fin,1,1)))
    ins2 <- setdiff(ins, era)
    expected <- lapply(fin, function(f) {
      ex <- grep(paste0("^", f), ins2, value=TRUE)
      if(length(ex) == 0) return(NULL)
      data.frame(query = f, target = ex, stringsAsFactors = F)
    })
    expected <- do.call(rbind, expected)
    expected <- arrange_result(expected)
    x$insert(ins)
    stopifnot(x$validate())
    x$erase(era)
    stopifnot(x$validate())
    results <- arrange_result(x$prefix_search(fin))
    stopifnot(identical(results, expected))
  })

  print(paste0("Checking multithreaded hamming search correctness for ", tt))
  local({
    x <- RadixForest$new()
    target <- unique(c(random_strings(NSEQS, CHARSET),""))
    query <- sample(c(sample(target, NSEQS/1000), random_strings(NSEQS/1000, CHARSET)))
    query <- unique(c(mutate_strings(query, indel_prob=0, charset = CHARSET), ""))
    x$insert(target)
    stopifnot(x$validate())
    results_dist <- arrange_result(x$search(query, max_distance = MAXDIST, mode = "hamming", nthreads=NTHREADS, show_progress=TRUE))
    results_frac <- arrange_result(x$search(query, max_fraction = MAXFRAC, mode = "hamming", nthreads=NTHREADS, show_progress=TRUE))
    sd_results <- dist_matrix_search(query, target, mode = "hamming")
    sd_dist <- arrange_result(sd_results[sd_results$distance <= MAXDIST, , drop = FALSE])
    sd_frac <- arrange_result(sd_results[sd_results$distance <= nchar(sd_results$query) * MAXFRAC, , drop = FALSE])
    stopifnot(identical(results_dist, sd_dist))
    stopifnot(identical(results_frac, sd_frac))
  })

    print(paste0("Checking multithreaded levenshtein search correctness for ", tt))
    local({
      x <- RadixForest$new()
      target <- unique(c(random_strings(NSEQS, CHARSET),""))
      query <- sample(c(sample(target, NSEQS/1000), random_strings(NSEQS/1000, CHARSET)))
      query <- unique(c(mutate_strings(query, charset = CHARSET), ""))
      x$insert(target)
      stopifnot(x$validate())
      results_dist <- arrange_result(x$search(query, max_distance = MAXDIST, mode = "levenshtein", nthreads=NTHREADS, show_progress=TRUE))
      results_frac <- arrange_result(x$search(query, max_fraction = MAXFRAC, mode = "levenshtein", nthreads=NTHREADS, show_progress=TRUE))
      sd_results <- dist_matrix_search(query, target, mode = "levenshtein")
      sd_dist <- arrange_result(sd_results[sd_results$distance <= MAXDIST, , drop = FALSE])
      sd_frac <- arrange_result(sd_results[sd_results$distance <= nchar(sd_results$query) * MAXFRAC, , drop = FALSE])
      stopifnot(identical(results_dist, sd_dist))
      stopifnot(identical(results_frac, sd_frac))
    })  
}

}
