library(dplyr)
library(tidyr)
library(stringdist)
library(stringr)
library(data.table)
library(treedist)
library(stringfish)
library(Biostrings)


x <- td_prefix_tree(c("ABCDEF", "ABDEF", "ABCZF"))
td_levenshtein(c("ABCDEF", "ABDEF", "ABCZF"), x, symmetric=F, max_distance = 10)

s2 <- c("", random_strings(1e3, 40, charset = "GCAT", vector_mode = "normal"))

mat2df <- function(res) {
  colnames(res) <- colnames(res) <- 1:nrow(res)
  res %>% 
    as.data.frame %>%
    tibble::rownames_to_column("target") %>% 
    pivot_longer(cols = -target, names_to = "query", values_to = "distance") %>%
    mutate(query = as.integer(query), target = as.integer(target)) %>%
    filter(query > target) %>% 
    dplyr::select(query, target, distance) %>%
    arrange(query, target) %>% as.data.frame
}

partial_hamming_pairwise <- function(query, target) {
  smat <- nucleotideSubstitutionMatrix(match = 0, mismatch = -1, baseOnly = TRUE)
  results <- list()
  for(i in 1:length(target)) {
    dist <- pairwiseAlignment(query, target[i], scoreOnly=T, substitutionMatrix=smat, gapOpening=1000, gapExtension=1000, type = "global-local")
    results[[i]] <- data.frame(query = 1:length(query), target = i, distance = -as.integer(dist))
  }
  rbindlist(results) %>% arrange(query, target) %>% as.data.frame
}

# levenshtein
res <- stringdistmatrix(s2, method = "lv", nthread = 8)
res <- mat2df(as.matrix(res)) %>% arrange(query, target)
x <- td_prefix_tree(s2)
res2 <- td_levenshtein(s2, x, nthreads=8, symmetric = T, max_distance = .Machine$integer.max)
res2 <- res2 %>% arrange(query, target)
stopifnot(all(res == res2))

# hamming
res <- stringdistmatrix(s2, method = "hamming", nthread = 8)
res <- mat2df(as.matrix(res)) %>% arrange(query, target) %>% filter(is.finite(distance))
x <- td_prefix_tree(s2)
res2 <- td_hamming(s2, x, nthreads=8, symmetric = T, max_distance = .Machine$integer.max)
res2 <- res2 %>% arrange(query, target)
stopifnot(all(res == res2))

# partial hamming
s2 <- random_strings(1e3, 40, charset = "GCAT", vector_mode = "normal")
x <- td_prefix_tree(s2)
s3 <- random_strings(100, 20, charset = "GCAT", vector_mode = "normal")
res2 <- td_partial_hamming(s3, x, nthreads=8, max_distance = 100)
res2 <- res2 %>% arrange(query, target, distance) %>% distinct(query, target, .keep_all=T)
res <- partial_hamming_pairwise(s3, s2)
stopifnot(all(res == res2))



# partial levenshtein
s2 <- random_strings(10, 20, charset = "GCAT", vector_mode = "normal")
x <- td_prefix_tree(s2)
s3 <- str_sub(s2, 2, -4)
td_partial_levenshtein(s3, x, anchor = "right", max_distance=3, nthreads=1)
td_partial_levenshtein(s3, x, anchor = "left", max_distance=3, nthreads=1)
td_partial_levenshtein(s3, x, anchor = "none", max_distance=3, nthreads=1)




#################################################

# grid <- expand.grid(test_size = c(100,300,1000,3000,10000), 
#                     method = c("levenshtein pairwise", "levenshtein trie max_dist=Inf", "levenshtein trie max_dist=3", "levenshtein trie max_dist=5", "levenshtein trie max_dist=1"), 
#                     nthread = 8, rep = 1:3)
# grid$time <- NA_real_
# for(i in 1:nrow(grid)) {
#   print(i)
#   z <- sample(seqs, grid$test_size[i])
#   if(grid$method[i] == "levenshtein pairwise") {
#     time <- Sys.time()
#     res <- stringdistmatrix(seqs[1:grid$test_size[i]], method = "lv", nthread = grid$nthread[i])
#     res <- as.matrix(res)
#     grid$time[i] <- as.numeric(Sys.time() - time, units = "secs")
#   } else if(grid$method[i] == "levenshtein trie max_dist=Inf") {
#     time <- Sys.time()
#     x <- create_trie(seqs[1:grid$test_size[i]])
#     res <- trie_levenshtein_matrix(x, seqs[1:grid$test_size[i]], nthreads = grid$nthread[i])
#     grid$time[i] <- as.numeric(Sys.time() - time, units = "secs")
#   } else if(grid$method[i] == "levenshtein trie max_dist=3") {
#     time <- Sys.time()
#     x <- create_trie(seqs[1:grid$test_size[i]])
#     res <- trie_levenshtein_sparse(x, seqs[1:grid$test_size[i]], max_distance = 3, nthreads = grid$nthread[i])
#     grid$time[i] <- as.numeric(Sys.time() - time, units = "secs")
#   } else if(grid$method[i] == "levenshtein trie max_dist=5") {
#     time <- Sys.time()
#     x <- create_trie(seqs[1:grid$test_size[i]])
#     res <- trie_levenshtein_sparse(x, seqs[1:grid$test_size[i]], max_distance = 5, nthreads = grid$nthread[i])
#     grid$time[i] <- as.numeric(Sys.time() - time, units = "secs")
#   } else if(grid$method[i] == "levenshtein trie max_dist=1") {
#     time <- Sys.time()
#     x <- create_trie(seqs[1:grid$test_size[i]])
#     res <- trie_levenshtein_sparse(x, seqs[1:grid$test_size[i]], max_distance = 1, nthreads = grid$nthread[i])
#     grid$time[i] <- as.numeric(Sys.time() - time, units = "secs")
#   }
#   rm(res)
#   gc()
# }
# 
# 
# ggplot(grid, aes(x = test_size, y = time, color = method, lty = ifelse(method == "levenshtein pairwise", 1, 3))) + 
#   geom_point() + geom_smooth(method = "lm", fill = NA, lwd=0.6) + theme_bw() +
#   scale_linetype_identity() + 
#   scale_x_log10() + scale_y_log10() + 
#   labs(x = "# of CDR3 sequences", y = "Time (s)", subtitle = "levenshtein distance, nthreads = 8")
