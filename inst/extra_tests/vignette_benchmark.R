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

sample_dna <- function(x, n, seed) {
  set.seed(seed)
  sample(x, min(n, length(x)), replace = FALSE)
}

prepare_dna_pool <- function(x) {
  x <- unique(toupper(x[!is.na(x) & nzchar(x)]))
  x <- x[grepl("^[ACGTN]+$", x)]
  if (!length(x)) {
    stop("No DNA sequences available for benchmarking")
  }
  x
}

bench_sizes <- function(max_nseqs) {
  value <- Sys.getenv("SEQTRIE_BENCH_SIZES", unset = "")
  if (nzchar(value)) {
    out <- suppressWarnings(as.integer(strsplit(value, "[,[:space:]]+")[[1L]]))
    if (any(is.na(out)) || any(out <= 0L)) {
      stop("SEQTRIE_BENCH_SIZES must be a comma- or space-separated list of positive integers")
    }
    if (any(out > max_nseqs)) {
      stop("SEQTRIE_BENCH_SIZES must not contain values larger than SEQTRIE_BENCH_MAX_NSEQS")
    }
    return(sort(unique(out)))
  }

  anchors <- c(1000L, 2000L, 5000L, 10000L, 20000L, 50000L, 100000L, max_nseqs)
  sort(unique(anchors[anchors <= max_nseqs]))
}

ensure_parent_dir <- function(path) {
  parent <- dirname(path)
  if (nzchar(parent) && parent != ".") {
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  }
}

format_count <- function(x) {
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

log_breaks <- function(x, multipliers = c(1, 3)) {
  x <- x[x > 0]
  powers <- seq(floor(log10(min(x))), ceiling(log10(max(x))))
  breaks <- as.vector(outer(multipliers, 10^powers))
  sort(unique(breaks[breaks >= min(x) & breaks <= max(x)]))
}

format_seconds <- function(x) {
  trim_number <- function(value, digits) {
    sub("\\.?0+$", "", sprintf(paste0("%.", digits, "f"), value))
  }
  ifelse(
    x < 1,
    trim_number(x, 3L),
    trim_number(x, 1L)
  )
}

plot_summary <- function(summary, outfile, sizes, niter, nthreads, max_distance, max_available) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required to generate the vignette benchmark plot")
  }

  classes <- c("radix_tree", "radix_forest", "star_tree")
  colors <- c("#2563EB", "#16A34A", "#F97316")
  names(colors) <- classes
  shapes <- c(radix_tree = 16, radix_forest = 17, star_tree = 15)

  positive_elapsed <- summary$median_elapsed[summary$median_elapsed > 0]
  if (!length(positive_elapsed)) {
    stop("All benchmark timings were zero; increase sequence sizes or iterations")
  }
  y_floor <- min(positive_elapsed) / 2
  y_values <- pmax(summary$median_elapsed, y_floor)
  y_lim <- range(y_values, na.rm = TRUE)
  y_lim <- c(y_lim[1] / 1.25, y_lim[2] * 1.4)
  y_ticks <- log_breaks(y_lim, multipliers = c(1, 3))
  x_ticks <- log_breaks(range(sizes), multipliers = c(1, 2, 5))

  summary$tree_class <- factor(summary$tree_class, levels = classes)

  ensure_parent_dir(outfile)

  plot <- ggplot2::ggplot(
    summary,
    ggplot2::aes(
      x = nseqs,
      y = pmax(median_elapsed, y_floor),
      color = tree_class,
      shape = tree_class,
      group = tree_class
    )
  ) +
    ggplot2::geom_line(linewidth = 0.9, alpha = 0.9) +
    ggplot2::geom_point(size = 2.6, stroke = 0.4) +
    ggplot2::scale_x_log10(
      breaks = x_ticks,
      labels = format_count,
      expand = ggplot2::expansion(mult = c(0.04, 0.06))
    ) +
    ggplot2::scale_y_log10(
      breaks = y_ticks,
      labels = format_seconds,
      limits = y_lim,
      expand = ggplot2::expansion(mult = c(0.02, 0.04))
    ) +
    ggplot2::scale_color_manual(values = colors, drop = FALSE) +
    ggplot2::scale_shape_manual(values = shapes, drop = FALSE) +
    ggplot2::labs(
      title = "Global edit-distance self-join",
      x = "Number of sequences",
      y = "Median elapsed time (seconds)",
      color = NULL,
      shape = NULL,
      caption = sprintf(
        "max_distance = %d; nthreads = %d; %d iteration%s; shown to %s of %s CDR3 sequences",
        max_distance,
        nthreads,
        niter,
        if (niter == 1L) "" else "s",
        format_count(max(sizes)),
        format_count(max_available)
      )
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 17,
        hjust = 0.5,
        color = "#111827",
        margin = ggplot2::margin(b = 14)
      ),
      axis.title = ggplot2::element_text(color = "#111827", size = 12),
      axis.text = ggplot2::element_text(color = "#374151", size = 10),
      panel.grid.major.y = ggplot2::element_line(color = "#E5E7EB", linewidth = 0.45),
      panel.grid.major.x = ggplot2::element_line(color = "#F1F5F9", linewidth = 0.35),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "#9CA3AF", linewidth = 0.35),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.text = ggplot2::element_text(color = "#111827", size = 10),
      legend.key.width = grid::unit(22, "pt"),
      plot.caption = ggplot2::element_text(hjust = 0.5, size = 9, color = "#4B5563"),
      plot.margin = ggplot2::margin(14, 18, 12, 14)
    )

  ggplot2::ggsave(
    outfile,
    plot = plot,
    width = 8.2,
    height = 5.2,
    dpi = 300
  )
}

NITER <- bench_param("NITER", 3L)
NTHREADS <- bench_param("NTHREADS", 8L)
MAX_DISTANCE <- bench_param("MAX_DISTANCE", 3L)
SEED <- bench_param("SEED", 314156L)

pool <- prepare_dna_pool(covid_cdr3)
MAX_NSEQS <- bench_param("MAX_NSEQS", length(pool))
MAX_NSEQS <- min(MAX_NSEQS, length(pool))
SIZES <- bench_sizes(MAX_NSEQS)

OUT_CSV <- Sys.getenv("SEQTRIE_BENCH_OUT", unset = "")
OUT_PNG <- Sys.getenv("SEQTRIE_BENCH_PNG", unset = "vignettes/vignette_benchmark.png")

methods <- c(
  radix_tree = "radix_tree",
  radix_forest = "radix_forest",
  star_tree = "star_tree"
)

grid <- expand.grid(
  iter = seq_len(NITER),
  nseqs = SIZES,
  tree_class = names(methods),
  stringsAsFactors = FALSE
)
set.seed(SEED)
grid <- grid[sample.int(nrow(grid)), ]
grid$nthreads <- NTHREADS
grid$max_distance <- MAX_DISTANCE
grid$elapsed <- NA_real_
grid$matches <- NA_integer_
grid$distance_sum <- NA_real_

samples <- lapply(seq_len(NITER), function(iter) {
  sample_dna(pool, MAX_NSEQS, SEED + iter)
})

for (i in seq_len(nrow(grid))) {
  row <- grid[i, ]
  seqs <- samples[[row$iter]][seq_len(row$nseqs)]

  elapsed <- system.time({
    result <- dist_search(
      seqs,
      max_distance = MAX_DISTANCE,
      mode = "global",
      tree_class = methods[[row$tree_class]],
      nthreads = NTHREADS,
      show_progress = FALSE
    )
  })[["elapsed"]]

  grid$elapsed[i] <- elapsed
  grid$matches[i] <- nrow(result)
  grid$distance_sum[i] <- if (nrow(result)) sum(result$distance) else 0L

  cat(sprintf(
    "%s iter=%d nseqs=%d elapsed=%.3f matches=%d\n",
    row$tree_class,
    row$iter,
    length(seqs),
    elapsed,
    nrow(result)
  ))
  flush.console()

  rm(seqs, result)
  gc(full = TRUE)
}

result_check <- aggregate(
  cbind(matches, distance_sum) ~ iter + nseqs,
  data = grid,
  FUN = function(x) length(unique(x))
)
if (any(result_check$matches != 1L | result_check$distance_sum != 1L)) {
  warning("Benchmark methods did not return identical result summaries")
}

summary <- aggregate(
  elapsed ~ tree_class + nseqs,
  data = grid,
  FUN = function(x) c(mean = mean(x), median = median(x))
)
summary <- data.frame(
  tree_class = summary$tree_class,
  nseqs = summary$nseqs,
  mean_elapsed = round(summary$elapsed[, "mean"], 3L),
  median_elapsed = round(summary$elapsed[, "median"], 3L),
  stringsAsFactors = FALSE
)
summary <- summary[order(match(summary$tree_class, names(methods)), summary$nseqs), ]

if (nzchar(OUT_CSV)) {
  ensure_parent_dir(OUT_CSV)
  write.table(grid, file = OUT_CSV, sep = ",", row.names = FALSE)
}
plot_summary(
  summary,
  outfile = OUT_PNG,
  sizes = SIZES,
  niter = NITER,
  nthreads = NTHREADS,
  max_distance = MAX_DISTANCE,
  max_available = length(pool)
)

old_width <- getOption("width")
options(width = max(old_width, 120L))
print(summary, row.names = FALSE, right = FALSE)
options(width = old_width)
