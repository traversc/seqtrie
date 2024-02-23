library(dplyr)
library(data.table)
library(ggplot2)
library(patchwork)
library(this.path)
setwd(dirname(this.path())) # ...seqtrie/inst/extra_tests

array_sizes <- seq(8,64,by = 8) %>% rep(each = 20) %>% sample
lapply(array_sizes, function(AS) {
  Sys.setenv("SEQTRIE_SMALL_ARRAY_SIZE"=AS)
  system("cd ../../ && make install-fast")
  res <- system2("/usr/bin/time",  args=c("-v", "Rscript simple_benchmark.R"), stdout=T, stderr=T)
  mem <- grep("Maximum resident set size", res, value=T) %>%
    gsub(".+:", "", .) %>%
    as.numeric
  data.frame(array_size=AS, time = as.numeric(res[1]), mem_usage = mem)
}) %>% rbindlist -> results

results <- filter(results, array_size != 1000)

results2 <- arrange(results, array_size) %>%
  group_by(array_size) %>%
  summarize(mem_usage = median(mem_usage), time = median(time))

g1 <- ggplot(results, aes(x = array_size, y = time)) + 
  geom_point(pch=21, color = "chartreuse") + 
  geom_line(data=results2, color = "chartreuse") + 
  scale_x_continuous(breaks = array_sizes) + 
  theme_bw(base_size=14)

g2 <- ggplot(results, aes(x = array_size, y = mem_usage)) + 
  geom_point(pch=21, color = "darkorange") + 
  geom_line(data=results2, color = "darkorange") + 
  scale_x_continuous(breaks = array_sizes) + 
  theme_bw(base_size=14)

g <- g1 + g2 + plot_layout(ncol=1)
ggsave(g, file = "small_array_size_bench.png", width=6, height=8)


