suppressPackageStartupMessages({
  library(seqtrie)
})

data(covid_cdr3, package = "seqtrie")

bench_param <- function(name, default) {
  value <- Sys.getenv(paste0("SEQTRIE_BENCH_", name), unset = "")
  if (!nzchar(value)) {
    return(default)
  }
  out <- suppressWarnings(as.integer(value))
  if (is.na(out)) {
    stop("SEQTRIE_BENCH_", name, " must be an integer")
  }
  out
}

translate_cdr3_to_aa <- function(seqs) {
  codon_table <- c(
    TTT = "F", TTC = "F", TTA = "L", TTG = "L",
    TCT = "S", TCC = "S", TCA = "S", TCG = "S",
    TAT = "Y", TAC = "Y", TAA = "*", TAG = "*",
    TGT = "C", TGC = "C", TGA = "*", TGG = "W",
    CTT = "L", CTC = "L", CTA = "L", CTG = "L",
    CCT = "P", CCC = "P", CCA = "P", CCG = "P",
    CAT = "H", CAC = "H", CAA = "Q", CAG = "Q",
    CGT = "R", CGC = "R", CGA = "R", CGG = "R",
    ATT = "I", ATC = "I", ATA = "I", ATG = "M",
    ACT = "T", ACC = "T", ACA = "T", ACG = "T",
    AAT = "N", AAC = "N", AAA = "K", AAG = "K",
    AGT = "S", AGC = "S", AGA = "R", AGG = "R",
    GTT = "V", GTC = "V", GTA = "V", GTG = "V",
    GCT = "A", GCC = "A", GCA = "A", GCG = "A",
    GAT = "D", GAC = "D", GAA = "E", GAG = "E",
    GGT = "G", GGC = "G", GGA = "G", GGG = "G"
  )

  seqs <- toupper(seqs)
  keep <- !is.na(seqs) & nchar(seqs) > 0L & nchar(seqs) %% 3L == 0L & grepl("^[ACGT]*$", seqs)
  seqs <- seqs[keep]

  translate_one <- function(seq) {
    starts <- seq.int(1L, nchar(seq), by = 3L)
    codons <- substring(seq, starts, starts + 2L)
    aa <- unname(codon_table[codons])
    if (any(is.na(aa)) || any(aa == "*")) {
      return(NA_character_)
    }
    paste0(aa, collapse = "")
  }

  aa <- vapply(seqs, translate_one, character(1))
  unique(aa[!is.na(aa)])
}

sample_vec <- function(x, n) {
  sample(x, n, replace = length(x) < n)
}

make_cost_matrix <- function(query, target) {
  chars <- unlist(strsplit(c(query, target), "", fixed = TRUE), use.names = FALSE)
  chars <- sort(unique(chars[nzchar(chars)]))
  seqtrie::generate_cost_matrix(paste0(chars, collapse = ""), match = 0L, mismatch = 1L)
}

summarize_result <- function(result) {
  list(
    matches = nrow(result),
    distance_sum = if (nrow(result)) sum(result$distance) else 0L
  )
}

NITER <- bench_param("NITER", 3L)
NTARGET <- bench_param("NTARGET", bench_param("NSEQS", 10000L))
NQUERY <- bench_param("NQUERY", bench_param("NSEQS", 10000L))
NTHREADS <- bench_param("NTHREADS", 4L)
MAX_DISTANCE <- bench_param("MAX_DISTANCE", 3L)
LINEAR_GAP_COST <- bench_param("LINEAR_GAP_COST", 2L)
AFFINE_GAP_COST <- bench_param("AFFINE_GAP_COST", 1L)
AFFINE_GAP_OPEN_COST <- bench_param("AFFINE_GAP_OPEN_COST", 2L)
SEED <- bench_param("SEED", 314156L)
VARIANT <- Sys.getenv("SEQTRIE_VARIANT", unset = "local")

datasets <- list(
  cdr3_nt = unique(covid_cdr3[!is.na(covid_cdr3) & nzchar(covid_cdr3)]),
  cdr3_aa = translate_cdr3_to_aa(covid_cdr3)
)

methods <- list(
  RadixTree_search_hamming = function(query, tree, forest, cost_matrix) {
    tree$search(query, max_distance = MAX_DISTANCE, mode = "hamming",
                nthreads = NTHREADS, show_progress = FALSE)
  },
  RadixForest_search_hamming = function(query, tree, forest, cost_matrix) {
    forest$search(query, max_distance = MAX_DISTANCE, mode = "hamming",
                  nthreads = NTHREADS, show_progress = FALSE)
  },
  RadixTree_search_global_unit = function(query, tree, forest, cost_matrix) {
    tree$search(query, max_distance = MAX_DISTANCE, mode = "levenshtein",
                nthreads = NTHREADS, show_progress = FALSE)
  },
  RadixForest_search_global_unit = function(query, tree, forest, cost_matrix) {
    forest$search(query, max_distance = MAX_DISTANCE, mode = "levenshtein",
                  nthreads = NTHREADS, show_progress = FALSE)
  },
  RadixTree_search_anchored_unit = function(query, tree, forest, cost_matrix) {
    tree$search(query, max_distance = MAX_DISTANCE, mode = "anchored",
                nthreads = NTHREADS, show_progress = FALSE)
  },
  RadixTree_search_global_linear = function(query, tree, forest, cost_matrix) {
    tree$search(query, max_distance = MAX_DISTANCE, mode = "global",
                cost_matrix = cost_matrix, gap_cost = LINEAR_GAP_COST,
                nthreads = NTHREADS, show_progress = FALSE)
  },
  RadixForest_search_global_linear = function(query, tree, forest, cost_matrix) {
    forest$search(query, max_distance = MAX_DISTANCE, mode = "global",
                  cost_matrix = cost_matrix, gap_cost = LINEAR_GAP_COST,
                  nthreads = NTHREADS, show_progress = FALSE)
  },
  RadixTree_search_anchored_linear = function(query, tree, forest, cost_matrix) {
    tree$search(query, max_distance = MAX_DISTANCE, mode = "anchored",
                cost_matrix = cost_matrix, gap_cost = LINEAR_GAP_COST,
                nthreads = NTHREADS, show_progress = FALSE)
  },
  RadixTree_search_global_affine = function(query, tree, forest, cost_matrix) {
    tree$search(query, max_distance = MAX_DISTANCE, mode = "global",
                cost_matrix = cost_matrix, gap_cost = AFFINE_GAP_COST,
                gap_open_cost = AFFINE_GAP_OPEN_COST,
                nthreads = NTHREADS, show_progress = FALSE)
  },
  RadixForest_search_global_affine = function(query, tree, forest, cost_matrix) {
    forest$search(query, max_distance = MAX_DISTANCE, mode = "global",
                  cost_matrix = cost_matrix, gap_cost = AFFINE_GAP_COST,
                  gap_open_cost = AFFINE_GAP_OPEN_COST,
                  nthreads = NTHREADS, show_progress = FALSE)
  },
  RadixTree_search_anchored_affine = function(query, tree, forest, cost_matrix) {
    tree$search(query, max_distance = MAX_DISTANCE, mode = "anchored",
                cost_matrix = cost_matrix, gap_cost = AFFINE_GAP_COST,
                gap_open_cost = AFFINE_GAP_OPEN_COST,
                nthreads = NTHREADS, show_progress = FALSE)
  },
  RadixTree_single_gap_search = function(query, tree, forest, cost_matrix) {
    tree$single_gap_search(query, max_distance = MAX_DISTANCE, gap_cost = 1L,
                           nthreads = NTHREADS, show_progress = FALSE)
  }
)

grid <- expand.grid(
  iter = seq_len(NITER),
  dataset = names(datasets),
  method = names(methods),
  stringsAsFactors = FALSE
)
set.seed(SEED)
grid <- grid[sample.int(nrow(grid)), ]
grid$variant <- VARIANT
grid$n_target <- NTARGET
grid$n_query <- NQUERY
grid$nthreads <- NTHREADS
grid$max_distance <- MAX_DISTANCE
grid$build_time <- NA_real_
grid$elapsed <- NA_real_
grid$matches <- NA_integer_
grid$distance_sum <- NA_real_

for (i in seq_len(nrow(grid))) {
  row <- grid[i, ]
  set.seed(SEED + row$iter * 1000L + match(row$dataset, names(datasets)) * 100L)

  pool <- datasets[[row$dataset]]
  target <- sample_vec(pool, NTARGET)
  query <- sample_vec(pool, NQUERY)
  cost_matrix <- make_cost_matrix(query, target)

  build_time <- system.time({
    tree <- seqtrie::RadixTree$new(target)
    forest <- seqtrie::RadixForest$new(target)
  })[["elapsed"]]

  elapsed <- system.time({
    result <- methods[[row$method]](query, tree, forest, cost_matrix)
  })[["elapsed"]]
  result_summary <- summarize_result(result)

  grid$build_time[i] <- build_time
  grid$elapsed[i] <- elapsed
  grid$matches[i] <- result_summary$matches
  grid$distance_sum[i] <- result_summary$distance_sum

  cat(sprintf("%s %s %s iter=%d elapsed=%.3f matches=%d\n",
              VARIANT, row$dataset, row$method, row$iter, elapsed, result_summary$matches))
  flush.console()

  rm(target, query, cost_matrix, tree, forest, result, result_summary)
  gc(full = TRUE)
}

summary <- aggregate(
  cbind(elapsed, matches) ~ dataset + method,
  data = grid,
  FUN = function(x) c(mean = mean(x), median = median(x))
)
summary <- data.frame(
  dataset = summary$dataset,
  method = summary$method,
  mean_time = round(summary$elapsed[, "mean"], 3L),
  median_time = round(summary$elapsed[, "median"], 3L),
  mean_matches = round(summary$matches[, "mean"], 1L),
  median_matches = round(summary$matches[, "median"], 1L),
  stringsAsFactors = FALSE
)
summary <- summary[order(summary$dataset, summary$method), , drop = FALSE]

outfile <- Sys.getenv("SEQTRIE_BENCH_OUT", unset = "")
if (nzchar(outfile)) {
  write.table(grid, file = outfile, sep = ",", row.names = FALSE,
              col.names = !file.exists(outfile), append = file.exists(outfile))
}

old_width <- getOption("width")
options(width = max(old_width, 200L))
print(summary, row.names = FALSE, right = FALSE)
options(width = old_width)
