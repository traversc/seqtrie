# export("DNATree_create", "RadixTree_create", "PrefixTree_create", 
#        "DNATree_size", "RadixTree_size", "PrefixTree_size",
#        "DNATree_insert", "RadixTree_insert", "PrefixTree_insert",
#        "DNATree_erase", "RadixTree_erase", "PrefixTree_erase",
#        "DNATree_levenshtein", "RadixTree_levenshtein", "PrefixTree_levenshtein",
#        "DNATree_hamming", "RadixTree_hamming", "PrefixTree_hamming")


library(seqtrie)
library(stringdist)
library(stringfish)
library(Rcpp)
library(dplyr)

NSEQS     <- 10000 # must be larger than 1000
NITER     <- 3
MAXSEQLEN <- 200
MAXDIST   <- MAXSEQLEN * 0.05
MAXFRAC   <- 0.05

tic <- function() { .time <<- Sys.time() }
toc <- function() { as.numeric(Sys.time() - .time, units = "secs") }


n_distinct <- function(x) {
  length(unique(x))
}

tree_equal <- function(x, y) {
  xs <- gsub("[0-9]+", "#", x$to_string())
  ys <- gsub("[0-9]+", "#", y$to_string())
  if(xs != ys) return(FALSE)
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
  x <- strsplit(x, "")
  sapply(x, function(a) {
    r <- runif(length(a)) < prob
    a[r] <- sample(charset, sum(r), replace=T)
    ins <- runif(length(a) < indel_prob)
    a[ins] <- paste0(sample(charset, sum(ins), replace=T), sample(charset, sum(ins), replace=T))
    del <- runif(length(a) < indel_prob)
    a[del] <- ""
    paste0(a, collapse = "")
  })
}

sd_search <- function(query, target, method = "lv") {
  results <- stringdist::stringdistmatrix(query, target, method = method, nthread=1)
  results <- data.frame(query = rep(query, times=length(target)), 
                        target = rep(target, each=length(query)), 
                        distance = as.vector(results), stringsAsFactors = F)
  results <- dplyr::filter(results, is.finite(distance))
  results$distance <- as.integer(results$distance)
  dplyr::arrange(results, query, target)
}

tt <- "RadixTree"

for(. in 1:NITER) {

  print(paste0("Checking correct insert/erase methods for ", tt))

    x <- RadixTree$new()
    y <- RadixTree$new()

    ins <- c(random_strings(NSEQS, "ACGT"),"")
    era <- c(sample(c(sample(ins, NSEQS/10), random_strings(NSEQS/10, "ACGT"))),"")
    x$insert(ins)
    stopifnot(x$validate())
    stopifnot(x$size() == n_distinct(ins))
    x$erase(era)
    stopifnot(x$validate())
    stopifnot(x$size() == n_distinct(ins[!ins %in% era]))
    
    y$insert(ins[!ins %in% era])
    stopifnot(y$validate())
    stopifnot(tree_equal(x, y))

  print(paste0('Checking find for ', tt))

    x <- RadixTree$new()

    ins <- c(random_strings(NSEQS, "ACGT"),"")
    era <- c(sample(c(sample(ins, NSEQS/10), random_strings(NSEQS/10, "ACGT"))))
    fin <- c(sample(c(sample(ins, NSEQS/10), random_strings(NSEQS/10, "ACGT"))),"")
    expected <- fin %in% setdiff(ins, era)
    x$insert(ins)
    stopifnot(x$validate())
    x$erase(era)
    stopifnot(x$validate())
    results <- x$find(fin)
    stopifnot(identical(results, expected))

  print(paste0('Checking prefix_search for ', tt))

    x <- RadixTree$new()
    
    ins <- c(random_strings(NSEQS, "ACGT"),"")
    era <- c(sample(c(sample(ins, NSEQS/10), random_strings(NSEQS/10, "ACGT"))))
    fin <- c(sample(c(sample(ins, NSEQS/1000), random_strings(NSEQS/1000, "ACGT"))),"") %>% substr(1,5)
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

  
  print(paste0("Checking levenshtein search correctness for ", tt))

    x <- RadixTree$new()

    target <- c(random_strings(NSEQS, "ACGT"),"") %>% unique
    query <- sample(c(sample(target, NSEQS/1000), random_strings(NSEQS/1000, "ACGT")))
    query <- c(mutate_strings(query, charset = "ACGT"), "") %>% unique
    x$insert(target)
    stopifnot(x$validate())
    results_dist <- x$search(query, max_distance = MAXDIST, mode = "levenshtein", show_progress=TRUE) %>% dplyr::arrange(query, target)
    results_frac <- x$search(query, max_fraction = MAXFRAC, mode = "levenshtein", show_progress=TRUE) %>% dplyr::arrange(query, target)
    sd_results <- sd_search(query, target, method = "lv")
    sd_dist <- dplyr::filter(sd_results, distance <= MAXDIST)
    sd_frac <- dplyr::filter(sd_results, distance <= nchar(query) * MAXFRAC)
    stopifnot(identical(results_dist, sd_dist))
    stopifnot(identical(results_frac, sd_frac))
  
  print(paste0("Checking hamming search correctness for ", tt))

    x <- RadixTree$new()

    target <- c(random_strings(NSEQS, "ACGT"),"") %>% unique
    query <- sample(c(sample(target, NSEQS/1000), random_strings(NSEQS/1000, "ACGT")))
    query <- c(mutate_strings(query, indel_prob=0, charset = "ACGT"), "") %>% unique
    x$insert(target)
    stopifnot(x$validate())
    results_dist <- x$search(query, max_distance = MAXDIST, mode = "hamming", show_progress=TRUE) %>% dplyr::arrange(query, target)
    results_frac <- x$search(query, max_fraction = MAXFRAC, mode = "hamming", show_progress=TRUE) %>% dplyr::arrange(query, target)
    sd_results <- sd_search(query, target, method = "hamming")
    sd_dist <- dplyr::filter(sd_results, distance <= MAXDIST)
    sd_frac <- dplyr::filter(sd_results, distance <= nchar(query) * MAXFRAC)
    stopifnot(identical(results_dist, sd_dist))
    stopifnot(identical(results_frac, sd_frac))

  
  print(paste0("Checking multithreaded levenshtein search correctness for ", tt))

    x <- RadixTree$new()

    target <- c(random_strings(NSEQS, "ACGT"),"") %>% unique
    query <- sample(c(sample(target, NSEQS/1000), random_strings(NSEQS/1000, "ACGT")))
    query <- c(mutate_strings(query, charset = "ACGT"), "") %>% unique
    x$insert(target)
    stopifnot(x$validate())
    results_dist <- x$search(query, max_distance = MAXDIST, mode = "levenshtein", nthreads=4, show_progress=TRUE) %>% dplyr::arrange(query, target)
    results_frac <- x$search(query, max_fraction = MAXFRAC, mode = "levenshtein", nthreads=4, show_progress=TRUE) %>% dplyr::arrange(query, target)
    sd_results <- sd_search(query, target, method = "lv")
    sd_dist <- dplyr::filter(sd_results, distance <= MAXDIST)
    sd_frac <- dplyr::filter(sd_results, distance <= nchar(query) * MAXFRAC)
    stopifnot(identical(results_dist, sd_dist))
    stopifnot(identical(results_frac, sd_frac))
  
  print(paste0("Checking multithreaded hamming search correctness for ", tt))

    x <- RadixTree$new()

    target <- c(random_strings(NSEQS, "ACGT"),"") %>% unique
    query <- sample(c(sample(target, NSEQS/1000), random_strings(NSEQS/1000, "ACGT")))
    query <- c(mutate_strings(query, indel_prob=0, charset = "ACGT"), "") %>% unique
    x$insert(target)
    stopifnot(x$validate())
    results_dist <- x$search(query, max_distance = MAXDIST, mode = "hamming", nthreads=4, show_progress=TRUE) %>% dplyr::arrange(query, target)
    results_frac <- x$search(query, max_fraction = MAXFRAC, mode = "hamming", nthreads=4, show_progress=TRUE) %>% dplyr::arrange(query, target)
    sd_results <- sd_search(query, target, method = "hamming")
    sd_dist <- dplyr::filter(sd_results, distance <= MAXDIST)
    sd_frac <- dplyr::filter(sd_results, distance <= nchar(query) * MAXFRAC)
    stopifnot(identical(results_dist, sd_dist))
    stopifnot(identical(results_frac, sd_frac))
  
  
}
