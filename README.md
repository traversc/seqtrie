seqtrie
================

`seqtrie` is a collection of Radix Tree algorithms. These include some
classic algorithms like prefix lookup and more bioinformatics focused
algorithms such as alignment algorithms and sequence distance
calculations (Hamming and Levenshtein distances).

A trie (aka Prefix Tree) is a data structure used in many applications,
notably for bioinformatics in the BLAST algorithm for K-mer lookup. It
is a space efficient data structure where a collection of sequences are
stored in a tree structure, where each leaf represents one sequence, and
each node holds one character representing a shared prefix of all
sequences descending from it.

A Radix Tree (aka Compact Prefix Tree) is an improvement on a trie,
using less memory and being generally faster. In a Radix Tree, each node
is able to represent multiple characters instead of just one.

Tries and Radix Trees have similar complexity to a hashmap. Storing a
sequence within the tree or looking it up are O(k) where k is the length
of the sequence. The main advantage of a trie is that it can be used to
quickly find similar sequences by various algorithms and metrics,
whereas a hashmap does not contain any sequence similarity information.

See also: <https://en.wikipedia.org/wiki/Radix_tree>

In `seqtrie`, there are two `R6` classes:

`RadixTree` is the primary class in this package. There are three main
methods. The `$insert()` method is used to store sequences on the tree,
`$erase()` for erasing sequences from the tree and `$search()` for
finding similar sequences stored on the tree.

The second `R6` class is `RadixForest`, a derivative data structure
where separate trees are constructed for each sequence length. This data
structure has advantages and disadvantages, discussed later.

### Install

``` r
devtools::install_github("traversc/seqtrie")
```

### Simple example

To demonstrate the interface, below is a simple example where we insert
some sequences (strings), erase one and then plot out the tree.

``` r
tree <- RadixTree$new()
tree$insert(c("cargo", "cart", "carburetor", "carbuncle", "bar", "zebra"))
tree$erase("zebra")
tree$graph()
```

<img src="vignettes/simple_tree.png" title="simple_tree" width="576" />

### Levenshtein “edit distance” search

Below is an example using COVID19 T-cell data from Adaptive
Biotechnologies.
(<https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7418738/>. This data is
licensed under CC 4.0.)

Here, we find highly similar sequences within a fixed edit distance.

``` r
# create a new tree and insert 130,000 "CDR3" sequences
data(covid_cdr3)
tree <- RadixTree$new()
tree$insert(covid_cdr3)

# search for similar sequences within an edit distance of 2 ~ this may take a minute or two
results <- tree$search(covid_cdr3, max_distance=2, mode="levenshtein", nthreads=8)
```

### The output

The output is a data.frame mapping query (search sequences) and target
(sequences inserted into the tree).

``` r
tail(results)
#                                                query                                        target distance
# 262287          TGTGCCAGCAGCCCCGGGGACAATGAGCAGTTCTTC          TGTGCCAGCAGCCCGGGGACAATTGAGCAGTTCTTC        2
# 262288          TGTGCCAGCAGCCCCGGGGACAATGAGCAGTTCTTC          TGTGCCAGCAGCCCGGGGACTAATGAGCAGTTCTTC        2
# 262289          TGTGCCAGCAGCCCCGGGGACAATGAGCAGTTCTTC          TGTGCCAGCAGCCTCGGGTACAATGAGCAGTTCTTC        2
# 262290          TGTGCCAGCAGCCCCGGGGACAATGAGCAGTTCTTC          TGTGCCAGCAGCGCCCGGGACAATGAGCAGTTCTTC        2
# 262291 TGTGCCAGCAGCCCTCGACTAGTCTCCTACAATGAGCAGTTCTTC TGTGCCAGCAGCCCTCGACTAGTCTCCTACAATGAGCAGTTCTTC        0
# 262292 TGTGCCAGCAGCCTCGACAGGGGAAACGCCTACGAGCAGTACTTC TGTGCCAGCAGCCTCGACAGGGGAAACGCCTACGAGCAGTACTTC        0
```

### Search parameters

The `$search()` function contains two mutually exclusive parameters:
`max_distance` and `max_fraction`.

The former parameter sets an absolute threshold that limits the distance
between pairs of sequences in the output, and the latter parameter sets
a threshold relative to the query sequence length.

The search time is monotonically increasing with the distance threshold,
logarithmically increasing with tree size and linearly increasing with
the number of query sequences and the length of each sequence.

Overall, the algorithm is significantly faster than a pairwise/matrix
edit distance calculation for finding similar sequences.

However, *care still needs to be taken when setting parameters for
searching a large number of sequences (~100,000+).*

#### Some additional examples using the `max_fraction` parameter.

**Some additional examples using the max_fraction parameter.**

``` r
results <- tree$search(covid_cdr3, max_fraction=0.035, mode="levenshtein", nthreads=8) # ~several seconds
results <- tree$search(covid_cdr3, max_fraction=0.06, mode="levenshtein", nthreads=8) # ~1 minute
results <- tree$search(covid_cdr3, max_fraction=0.15, mode="levenshtein", nthreads=8) # ~15-20 minutes
```

### Hamming distance search

Hamming distance is similar to Levenshtein distance, but does not allow
insertions or deletions. Sequences must be the same length.

Because of this restriction, Hamming distance is generally a lot faster.

``` r
results <- tree$search(covid_cdr3, max_fraction=0.035, mode="hamming", nthreads=8) # ~1 second
results <- tree$search(covid_cdr3, max_fraction=0.06, mode="hamming", nthreads=8) # ~ several seconds
results <- tree$search(covid_cdr3, max_fraction=0.15, mode="hamming", nthreads=8) # ~ 1.5 minutes
```

### Anchored alignment searches

An anchored alignment is a form of semi-global alignment, where the
query sequence is “anchored” (global) to the beginning of both the query
and target sequences, but is semi-global in that the end of the either
the query sequence or target sequence (but not both) can be unaligned.
This type of alignment is sometimes called an “extension” alignment in
literature.

``` r
tree <- RadixTree$new()
tree$insert("CARTON")
tree$insert("CAR")
tree$insert("CARBON")
tree$search("CART", max_distance = 0, mode = "anchored")
print(results)
#   query target distance query_size target_size
# 1  CART    CAR        0          3           3
# 2  CART CARTON        0          4           4
```

Because the alignment is semi-global at the end of the alignment, the
query of “CART” finds “CAR” and “CARTON” but not “CARBON” given a max
distance of 0. Additionally, the output of an anchored search also
returns the position of the query and target at the ends. Either the
query or the target must fully align, so at least one of the end
positions will be the full length of the sequence.

This type of alignment is frequently useful in biology e.g. if you are
trying to align multiple reads that are variable in length but start at
the same genomic position or primer site.

### Custom distance searches and affine gap alignment

`seqtrie` supports custom alignment parameters, including an affine gap
parameter. The interface is similar to the
`Biostrings::pairwiseAlignment`. Note: we are calculating distance
(higher is worse) and not alignment score (higher is better).

``` r
tree <- RadixTree$new()
tree$insert(covid_cdr3)

# define a custom distance matrix - generate_cost_matrix is a convienence function
# gap and gap_open can be defined directly in the cost_matrix or as search method parameters
cost_mat <- generate_cost_matrix("ACGT", match=0, mismatch=1, gap=2, gap_open=1)
print(cost_mat)
#          A C G T gap gap_open
# A        0 1 1 1   2        1
# C        1 0 1 1   2        1
# G        1 1 0 1   2        1
# T        1 1 1 0   2        1
# gap      2 2 2 2  NA       NA
# gap_open 1 1 1 1  NA       NA

# Perform a search. "Mode" can be either global or anchored. 
results <- tree$search(covid_cdr3, max_distance=3, cost_matrix=cost_mat, mode="global", nthreads=8)
```

### Radix Forest for faster Levenshtein searches

The `RadixForest` class is a data structure holding a collection of
Radix Trees, where a separate tree is constructed for each sequence
length. The primary advantage of `RadixForest` is significantly faster
Levenshtein searches, because you can know sequence length up front. The
disadvantages are higher memory usage, and no support for custom
distance searches.

Below is a brief comparison:

``` r
# RadixTree ~ 45 seconds
tree <- RadixTree$new()
tree$insert(covid_cdr3)
results_tree <- tree$search(covid_cdr3, max_distance=2, mode="levenshtein", nthreads=8)

# RadixForest ~ 19 seconds
frst <- RadixForest$new()
frst$insert(covid_cdr3)
results_frst <- frst$search(covid_cdr3, max_distance=2, mode="levenshtein", nthreads=8)

# The results are the same, but order is not guaranteed
identical(
  results_tree %>% dplyr::arrange(query, target),
  results_frst %>% dplyr::arrange(query, target) )
```

### Finding strings that start with a pattern

The `$find_prefix()` function can be used to find similar sequences that
start with a pattern. This is one of the classic use cases of trie data
structures, for use as a database lookup and predictive text.

``` r
tree <- RadixTree$new()
tree$insert(c("cargo", "cart", "carburetor", "carbuncle", "bar"))
tree$prefix_search("car")

#   query     target
# 1   car  carbuncle
# 2   car carburetor
# 3   car       cart
# 4   car      cargo
```

### Why not just use Bowtie2, BWA or other fast alignment software?

There are no apples-to-apples comparisons thus no benchmark comparisons.

With NGS alignment software, you are looking for alignments of reads
(queries) *within* a genome reference (target). Here, we’re looking for
alignments from the query to the *full* target.

Compared to pairwise alignment software, calculating all alignment pairs
takes much longer, but on the other hand gives you more information.

Use the right tool for the right task. Trie-based methods are
competitive on the tasks they are designed for.
