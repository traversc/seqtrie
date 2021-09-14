treedist
================

<!-- <img src="hex.png" width = "130" height = "150" align="right" style="border:0px;padding:15px"> -->

<!-- [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/treedist)](https://cran.r-project.org/package=treedist) [![CRAN\_Downloads\_Badge](https://cranlogs.r-pkg.org/badges/treedist)](https://cran.r-project.org/package=treedist) -->

<!-- [![CRAN\_Downloads\_Total\_Badge](https://cranlogs.r-pkg.org/badges/grand-total/treedist)](https://cran.r-project.org/package=treedist) -->

`treedist` is a collection of Radix Tree and Trie-based algorithms for
calculating string distances (Hamming and Levenshtein distances).

There are three `R6` classes in this package:

  - `RadixTree` a Radix Tree implementation. See:
    <https://en.wikipedia.org/wiki/Radix_tree>
  - `PrefixTree` a Prefix Tree implementation (aka Trie). See:
    <https://en.wikipedia.org/wiki/Trie>
  - `DNATree` a Radix Tree specialization for DNA sequences (ACGT
    characters). More memory efficient.

All classes have exactly the same interface.

The primary class functions for interacting with a tree are `$insert()`
for inserting sequences on the tree, `$erase()` for erasing sequences
from the tree and `$search()` for finding similar sequences stored on
the tree.

### Levenshtein “edit distance” search

Below is an example using COVID19 T-cell data from Adaptive
Biotechnologies.
(<https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7418738/>. This data is
licensed under CC 4.0.)

Here, we find highly similar sequences within a fixed edit distance.

``` r
library(treedist)
library(dplyr)

data(covid_cdr3) # ~130,000 "CDR3" sequences

# average CDR3 length
nchar(covid_cdr3) %>% mean
# [1] 43.56821

tree <- DNATree$new()

# insert sequences into the tree
tree$insert(covid_cdr3)

# search for sequences within an edit distance of 3 -- this takes a few minutes
results <- tree$search(covid_cdr3, max_distance = 3, mode = "levenshtein", nthreads=8)
```

### The output

The output is a data.frame mapping query (search sequences) and target
(sequences inserted into the tree).

``` r
head(results)
                                          # query                                        target distance
# 1          TGTGCCAGCGGGCGTCATAATTCACCCCTCCACTTT             TGTGCCAGCGGCTATAATTCACCCCTCCACTTT        3
# 2          TGTGCCAGCGGGCGTCATAATTCACCCCTCCACTTT          TGTGCCAGCGGGCGTCATAATTCACCCCTCCACTTT        0
# 3          TGTGCCAGCGGGCGTCATAATTCACCCCTCCACTTT          TGTGCCAGCGGGGACAGTAATTCACCCCTCCACTTT        3
# 4    TGTGCCAGCAGGCCAGCGGGAGTCAGGAATGAGCAGTTCTTC    TGTGCCAGCAGGCCAGCGGGAGTCAGGAATGAGCAGTTCTTC        0
# 5    TGCAGTGCTAGCCCTAGGACCTTATGGGGGGAGCAGTTCTTC    TGCAGTGCTAGCCCTAGGACCTTATGGGGGGAGCAGTTCTTC        0
# 6 TGCGCCAGCAGCTTGGTAACTAGCGTTAAAGATACGCAGTATTTT TGCGCCAGCAGCTTGGTAACTAGCGTTAAAGATACGCAGTATTTT        0
```

### Search parameters

The `$search()` function contains two mutually exclusive parameters:
`max_distance` and `max_fraction`.

The former parameter sets an absolute threshold for finding similar
sequences, and the latter parameter sets a threshold relative to the
query sequence length.

The search time is monotonically increasing with the distance threshold,
logarithmically increasing with tree size and linearly increasing with
the number of query sequences and the length of each sequence.

Overall, the algorithm is significantly faster than a pairwise edit
distance calculation. However, *care still needs to be taken when
setting parameters for searching a large number of sequences
(\~100,000+).*

#### Some additional examples using the `max_fraction` parameter.

**max\_fraction = 0.035 takes several seconds**

``` r
results <- tree$search(covid_cdr3, max_fraction = 0.035, mode = "levenshtein", nthreads=8)
```

**max\_fraction = 0.06 takes \~1 minute**

``` r
results <- tree$search(covid_cdr3, max_fraction = 0.06, mode = "levenshtein", nthreads=8)
```

**max\_fraction = 0.15 takes \~15 minutes**

``` r
results <- tree$search(covid_cdr3, max_fraction = 0.15, mode = "levenshtein", nthreads=8)
```

### Hamming distance search

Hamming distance is similar to Levenshtein distance, but does not allow
insertions or deletions. Sequences must be the same length.

Because of this restriction, Hamming distance is generally a lot faster.

**max\_fraction = 0.035 takes \~1 second**

``` r
results <- tree$search(covid_cdr3, max_fraction = 0.035, mode = "hamming", nthreads=8)
```

**max\_fraction = 0.06 takes a few seconds**

``` r
results <- tree$search(covid_cdr3, max_fraction = 0.06, mode = "hamming", nthreads=8)
```

**max\_fraction = 0.15 takes \~1.5 minutes**

``` r
results <- tree$search(covid_cdr3, max_fraction = 0.15, mode = "hamming", nthreads=8)
```
