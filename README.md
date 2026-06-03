seqtrie
================

### Background

`seqtrie` is a collection of radix tree algorithms for finding similar
strings, with a focus on biological sequence data. Here, “similar”
usually means close under Hamming distance, Levenshtein distance, or
semi-global alignment.

A radix tree, also called a compressed trie, is an efficient data
structure for searching large sets of sequences. It stores shared
prefixes once, so a search can skip whole branches instead of checking
every sequence one by one.

For bounded-distance searches, this can be sub-linear in the number of
stored sequences because the tree can prune whole branches.

### Basic usage

One common job is a similarity join: find every pair of sequences within
a fixed distance. The example below has 133,033 TCR beta CDR3 DNA
sequences (Nolan et al. 2020).

``` r
data(covid_cdr3)
results <- dist_search(covid_cdr3, max_distance = 3,
                       nthreads = 8, tree_class = "StarTree")
```

Trie-based algorithms are a good fit for this kind of search. The tree
structure lets `seqtrie` rule out large groups of targets without
checking every sequence against every other sequence.

### Tree classes

`seqtrie` has three tree classes.

- **RadixTree** (`radix_tree()`): the classic radix tree and the most
  general-purpose class. It supports insertion, deletion, prefix lookup,
  Hamming search, Levenshtein search, custom costs, single-gap search,
  and anchored search.
- **RadixForest** (`radix_forest()`): many radix trees, split by
  sequence length. This can be faster for Hamming and Levenshtein
  searches, especially when length is informative.
- **StarTree** (`star_tree()`): a specialized DNA similarity join class.
  It is a modified version of the Starcode all-pairs search algorithm
  from Zorita, Cuscó, and Filion (2015), adapted to operate over a radix
  trie rather than being a direct reimplementation. It supports fixed
  global/Levenshtein joins, fixed anchored joins (`mode = "anchored"`),
  and fixed Hamming joins (`mode = "hamming"`). For small edit-distance
  DNA self-joins it typically outperforms the other two classes (see the
  benchmark below).

`dist_search()` is a convenience wrapper around these classes. It can
search for similarity between `query` and `target`, or within `query`
when `target` is `NULL`.

Below is a benchmark applying `dist_search` to the `covid_cdr3` dataset.

![](vignettes/vignette_benchmark.png "vignette_benchmark")

### Understanding Radix Trees

Here is a small radix tree. We insert a few strings, erase one, and plot
the result.

The useful part is the shared prefix. `cargo`, `cart`, `carburetor`, and
`carbuncle` all start with `car`, so the tree stores `car` once. After
that, it only stores the pieces where the strings split.

``` r
tree <- radix_tree()
insert(tree, c("cargo", "cart", "carburetor", "carbuncle", "bar", "zebra"))
erase(tree, "zebra")
# plot_tree requires igraph and ggplot2
set.seed(1); plot_tree(tree)
```

![](vignettes/simple_tree.png "simple_tree")

### Hamming distance search

Hamming distance is similar to Levenshtein distance, but it does not
allow insertions or deletions. Sequences must be the same length.
Because of this restriction, Hamming distance is generally a lot faster.

Use `max_fraction` to set the distance threshold relative to query
length. For example, `max_fraction = 0.06` keeps matches within roughly
6% Hamming distance.

``` r
results <- align_search(tree, covid_cdr3, max_fraction = 0.035,
                        mode = "hamming", nthreads = 2)
results <- align_search(tree, covid_cdr3, max_fraction = 0.06,
                        mode = "hamming", nthreads = 2)
results <- align_search(tree, covid_cdr3, max_fraction = 0.15,
                        mode = "hamming", nthreads = 2)
```

### Anchored alignment searches

An anchored alignment is a form of semi-global alignment. The query is
anchored to the beginning of both the query and target sequences, but
the end of either sequence can be unaligned.

``` r
tree <- radix_tree()
insert(tree, "CARTON")
insert(tree, "CAR")
insert(tree, "CARBON")
align_search(tree, "CART", max_distance = 0, mode = "anchored")
```

    ##   query target distance query_size target_size
    ## 1  CART    CAR        0          3           3
    ## 2  CART CARTON        0          4           4

With `max_distance = 0`, the query `"CART"` finds `"CAR"` and `"CARTON"`
but not `"CARBON"`. Anchored search also returns where the alignment
ends in the query and target. This is useful when reads have variable
length but start from the same genomic position or primer site.

### Custom distance searches and affine gap alignment

`seqtrie` supports custom substitution costs and affine gap penalties.

``` r
tree <- radix_tree()
insert(tree, covid_cdr3)

# Define a custom substitution matrix. Use generate_cost_matrix for convenience.
cost_mat <- generate_cost_matrix("ACGT", match = 0, mismatch = 5)
print(cost_mat)
```

    ##   A C G T
    ## A 0 5 5 5
    ## C 5 0 5 5
    ## G 5 5 0 5
    ## T 5 5 5 0

``` r
# Set gap penalties via parameters (not in the matrix):
# - Linear gaps: set gap_cost only
# - Affine gaps: set both gap_cost and gap_open_cost

# Linear example
results_linear <- align_search(tree, covid_cdr3, max_distance = 8,
                               mode = "global",
                               cost_matrix = cost_mat,
                               gap_cost = 2,
                               nthreads = 2)

# Affine example
results_affine <- align_search(tree, covid_cdr3, max_distance = 8,
                               mode = "global",
                               cost_matrix = cost_mat,
                               gap_cost = 2,
                               gap_open_cost = 5,
                               nthreads = 2)

results_linear[results_linear$query != results_linear$target, , drop = FALSE]
```

    ##                                              query
    ## 130           TGTGCCAGCAGTGGGGGGAACACTGAAGCTTTCTTT
    ## 159  TGTGCCAGCAGCTCTAGCGGGGTGAGCACAGATACGCAGTATTTT
    ## 165  TGTGCCAGCAGCTCTACAGGGGTGAGCACAGATACGCAGTATTTT
    ## 203  TGTGCCAGTAGTCTCGGACAGGGCCTCTCCTACGAGCAGTACTTC
    ## 204  TGTGCCAGTAGTCTCGGACAGGGCCTCTCCTACGAGCAGTACTTC
    ## 258        TGTGCCAGCAGTTATTCGAACACCGGGGAGCTGTTTTTT
    ## 286  TGTGCCAGTAGTGTCGGACAGGGCCTCTCCTACGAGCAGTACTTC
    ## 287  TGTGCCAGTAGTGTCGGACAGGGCCTCTCCTACGAGCAGTACTTC
    ## 292        TGTGCCAGCAGCTCGGGACTGTCCTACGAGCAGTACTTC
    ## 293        TGTGCCAGCAGCTCGGGACTGTCCTACGAGCAGTACTTC
    ## 313  TGTGCCAGTAGTATGGGACAGGGGCTATCCTACGAGCAGTACTTC
    ## 329  TGTGCCAGTAGTATTGGACAGGGGACATCCTACGAGCAGTACTTC
    ## 330  TGTGCCAGTAGTATTGGACAGGGGACATCCTACGAGCAGTACTTC
    ## 342           TGTGCCAGCAGCTTCGCTAGCAATGAGCAGTTCTTC
    ## 345        TGTGCCAGCAGCTACAGGGGCTCCTACGAGCAGTACTTC
    ## 428           TGTGCCAGCAGCCCTTCGGGAAATGAGCAGTTCTTC
    ## 542           TGTGCCAGCAGCTTCGGGTACAATGAGCAGTTCTTC
    ## 544           TGTGCCAGCAGCTTCGGGTACAATGAGCAGTTCTTC
    ## 560           TGTGCCAGCAGCCCTACGGGGGGTGAAGCTTTCTTT
    ## 567        TGTGCCAGCAGTTACCGGGGAAGCACTGAAGCTTTCTTT
    ## 615           TGTGCCAGCAGTTTATACGGAGAGACCCAGTACTTC
    ## 635        TGTGCCAGTCAGCTTATCAACACCGGGGAGCTGTTTTTT
    ## 659           TGTGCCAGCAGCTCCTTTACCTACGAGCAGTACTTC
    ## 667  TGTGCCAGTAGTATTGGACAGGGACTCTCCTACGAGCAGTACTTC
    ## 687           TGTGCCAGCAGTCCGGGGAACACTGAAGCTTTCTTT
    ## 688           TGTGCCAGCAGTCCGGGGAACACTGAAGCTTTCTTT
    ## 728  TGTGCCAGCGGGACTAGCGGGGGGGCCCTAGATACGCAGTATTTT
    ## 760           TGTGCCAGCAGTTTTAGGGGCTACGAGCAGTACTTC
    ## 782           TGTGCCAGCAGTCTTAGTGGAGAGACCCAGTACTTC
    ## 794           TGTGCCAGCAGTTTTGCGGGTTACGAGCAGTACTTC
    ## 818           TGTGCCAGCAGCTTATATACCTACGAGCAGTACTTC
    ## 919           TGTGCCAGCAGCTCGGATTCCTACGAGCAGTACTTC
    ## 929        TGTGCCAGCAGCTCGGGGACCTCCTACGAGCAGTACTTC
    ## 930        TGTGCCAGCAGCTCGGGGACCTCCTACGAGCAGTACTTC
    ## 971  TGTGCCAGTAGTCTCGGGACAGGGTTCTCCTACGAGCAGTACTTC
    ## 1019 TGTGCCAGTAGTAAGGGACAGGGCCTCTCCTACGAGCAGTACTTC
    ## 1031 TGTGCCAGCAGGGCTAGCGGGGGGGCTCCAGATACGCAGTATTTT
    ## 1035             TGTGCCAGCAGCCTCGGGGGTGAAGCTTTCTTT
    ##                                             target distance
    ## 130           TGTGCCAGCAGTCCGGGGAACACTGAAGCTTTCTTT        8
    ## 159  TGTGCCAGCAGCTCTACAGGGGTGAGCACAGATACGCAGTATTTT        4
    ## 165  TGTGCCAGCAGCTCTAGCGGGGTGAGCACAGATACGCAGTATTTT        4
    ## 203  TGTGCCAGTAGTCTCGGGACAGGGTTCTCCTACGAGCAGTACTTC        8
    ## 204  TGTGCCAGTAGTGTCGGACAGGGCCTCTCCTACGAGCAGTACTTC        4
    ## 258        TGTGCCAGTCAGCTTATCAACACCGGGGAGCTGTTTTTT        8
    ## 286  TGTGCCAGTAGTAAGGGACAGGGCCTCTCCTACGAGCAGTACTTC        8
    ## 287  TGTGCCAGTAGTCTCGGACAGGGCCTCTCCTACGAGCAGTACTTC        4
    ## 292        TGTGCCAGCAGCTCGGGGACCTCCTACGAGCAGTACTTC        8
    ## 293           TGTGCCAGCAGCTCGGATTCCTACGAGCAGTACTTC        6
    ## 313  TGTGCCAGTAGTATTGGACAGGGGACATCCTACGAGCAGTACTTC        8
    ## 329  TGTGCCAGTAGTATTGGACAGGGACTCTCCTACGAGCAGTACTTC        8
    ## 330  TGTGCCAGTAGTATGGGACAGGGGCTATCCTACGAGCAGTACTTC        8
    ## 342           TGTGCCAGCAGCTTCGGGTACAATGAGCAGTTCTTC        8
    ## 345        TGTGCCAGCAGCTCGGGGACCTCCTACGAGCAGTACTTC        8
    ## 428           TGTGCCAGCAGCTTCGGGTACAATGAGCAGTTCTTC        8
    ## 542           TGTGCCAGCAGCCCTTCGGGAAATGAGCAGTTCTTC        8
    ## 544           TGTGCCAGCAGCTTCGCTAGCAATGAGCAGTTCTTC        8
    ## 560              TGTGCCAGCAGCCTCGGGGGTGAAGCTTTCTTT        6
    ## 567           TGTGCCAGCAGTCCGGGGAACACTGAAGCTTTCTTT        6
    ## 615           TGTGCCAGCAGTCTTAGTGGAGAGACCCAGTACTTC        8
    ## 635        TGTGCCAGCAGTTATTCGAACACCGGGGAGCTGTTTTTT        8
    ## 659           TGTGCCAGCAGCTTATATACCTACGAGCAGTACTTC        8
    ## 667  TGTGCCAGTAGTATTGGACAGGGGACATCCTACGAGCAGTACTTC        8
    ## 687        TGTGCCAGCAGTTACCGGGGAAGCACTGAAGCTTTCTTT        6
    ## 688           TGTGCCAGCAGTGGGGGGAACACTGAAGCTTTCTTT        8
    ## 728  TGTGCCAGCAGGGCTAGCGGGGGGGCTCCAGATACGCAGTATTTT        8
    ## 760           TGTGCCAGCAGTTTTGCGGGTTACGAGCAGTACTTC        8
    ## 782           TGTGCCAGCAGTTTATACGGAGAGACCCAGTACTTC        8
    ## 794           TGTGCCAGCAGTTTTAGGGGCTACGAGCAGTACTTC        8
    ## 818           TGTGCCAGCAGCTCCTTTACCTACGAGCAGTACTTC        8
    ## 919        TGTGCCAGCAGCTCGGGACTGTCCTACGAGCAGTACTTC        6
    ## 929        TGTGCCAGCAGCTACAGGGGCTCCTACGAGCAGTACTTC        8
    ## 930        TGTGCCAGCAGCTCGGGACTGTCCTACGAGCAGTACTTC        8
    ## 971  TGTGCCAGTAGTCTCGGACAGGGCCTCTCCTACGAGCAGTACTTC        8
    ## 1019 TGTGCCAGTAGTGTCGGACAGGGCCTCTCCTACGAGCAGTACTTC        8
    ## 1031 TGTGCCAGCGGGACTAGCGGGGGGGCCCTAGATACGCAGTATTTT        8
    ## 1035          TGTGCCAGCAGCCCTACGGGGGGTGAAGCTTTCTTT        6

### StarTree for fixed DNA self-joins

`star_tree` is a fixed DNA-only search structure based on the Starcode
all-pairs search strategy described by Zorita, Cuscó, and Filion (2015).
It is a modified version of that strategy, adapted to operate over a
radix trie rather than being a direct reimplementation.

This class follows two Starcode ideas:

- A lossless prefilter: a k-mer hash check can exclude the vast majority
  of negative matches up front, allowing only plausible matches to
  search the tree.
- A caching strategy: sorted queries with shared prefixes can reuse part
  of the tree search state instead of starting from scratch.

The Starcode algorithm remains one of the fastest methods for
self-similarity joins at small edit distances on DNA.

Unlike the other two classes, the sequences and search parameters are
supplied at construction time, the self-join runs immediately, and the
object does not support later insertion or deletion.

``` r
st <- star_tree(c("ACGT", "ACGA", "AAAA", "AAAT"),
                max_distance = 1,
                mismatch_cost = 1,
                gap_cost = 1,
                nthreads = 2)
result(st)
```

    ##   query target distance
    ## 1  AAAT   AAAA        1
    ## 2  ACGT   ACGA        1

``` r
# Search another query set using the same fixed costs and threshold.
align_search(st, c("ACGT", "AAAC"))
```

    ##   query target distance
    ## 1  AAAC   AAAA        1
    ## 2  AAAC   AAAT        1
    ## 3  ACGT   ACGA        1
    ## 4  ACGT   ACGT        0

``` r
# The same path is available through dist_search().
dist_search(c("ACGT", "ACGA", "AAAA", "AAAT"),
            max_distance = 1,
            tree_class = "StarTree")
```

    ##   query target distance
    ## 1  AAAT   AAAA        1
    ## 2  ACGT   ACGA        1

StarTree supports global/Levenshtein-style, anchored, and Hamming DNA
alignment with `A`, `C`, `G`, `T`, and `N`. It does not support custom
substitution matrices, affine gaps, or per-query distance thresholds.
Alignment parameters are fixed when the tree is built; `align_search()`
only accepts `query`, `nthreads`, and `show_progress` for `star_tree`
objects.

### StarTree anchored mode for fixed anchored self-joins

`star_tree(..., mode = "anchored")` is the fixed-tree equivalent for
anchored alignment. It is intended for DNA sets where every sequence
starts at the same biological position but may end at different
positions. The construction-time self-join returns each unordered pair
once, and additional query sets can be searched against the same fixed
targets.

``` r
ast <- star_tree(c("ACGT", "ACG", "AAAA", "AA"),
                 max_distance = 1,
                 mode = "anchored",
                 mismatch_cost = 1,
                 gap_cost = 1,
                 nthreads = 2)
result(ast)
```

    ##   query target distance query_size target_size
    ## 1  AAAA     AA        0          2           2
    ## 2   ACG     AA        1          1           2
    ## 3  ACGT     AA        1          1           2
    ## 4  ACGT    ACG        0          3           3

``` r
align_search(ast, c("ACGT", "AA"))
```

    ##   query target distance query_size target_size
    ## 1    AA     AA        0          2           2
    ## 2    AA   AAAA        0          2           2
    ## 3    AA    ACG        1          2           1
    ## 4    AA   ACGT        1          2           1
    ## 5  ACGT     AA        1          1           2
    ## 6  ACGT    ACG        0          3           3
    ## 7  ACGT   ACGT        0          4           4

``` r
dist_search(c("ACGT", "ACG", "AAAA", "AA"),
            max_distance = 1,
            mode = "anchored",
            tree_class = "StarTree")
```

    ##   query target distance query_size target_size
    ## 1  AAAA     AA        0          2           2
    ## 2   ACG     AA        1          1           2
    ## 3  ACGT     AA        1          1           2
    ## 4  ACGT    ACG        0          3           3

Anchored StarTree mode supports scalar mismatch and linear gap costs
through `mismatch_cost` and `gap_cost`. It does not support custom
substitution matrices, affine gaps, or per-query distance thresholds.

### StarTree Hamming mode for fixed same-length joins

`star_tree(..., mode = "hamming")` is a substitution-only fixed join:
only sequences of the same length can match, and `max_distance` is the
maximum number of mismatched positions. Because no insertions or
deletions are considered, `mismatch_cost` and `gap_cost` do not apply,
and it is typically much faster than global mode on the same data.

``` r
hst <- star_tree(c("ACGT", "ACGA", "TCGT", "ACG"),
                 max_distance = 1,
                 mode = "hamming",
                 nthreads = 2)
result(hst)
```

    ##   query target distance
    ## 1  ACGT   ACGA        1
    ## 2  TCGT   ACGT        1

``` r
align_search(hst, c("ACGT", "TTGT"))
```

    ##   query target distance
    ## 1  ACGT   ACGA        1
    ## 2  ACGT   ACGT        0
    ## 3  ACGT   TCGT        1
    ## 4  TTGT   TCGT        1

``` r
dist_search(c("ACGT", "ACGA", "TCGT", "ACG"),
            max_distance = 1,
            mode = "hamming",
            tree_class = "StarTree")
```

    ##   query target distance
    ## 1  ACGT   ACGA        1
    ## 2  TCGT   ACGT        1

`N` is treated as a regular base that mismatches every base, including
another `N`, so two aligned `N` positions still count as a mismatch.

### References

- Radix tree. Wikipedia. <https://en.wikipedia.org/wiki/Radix_tree>
- Nolan, Sean; Vignali, Marissa; Klinger, Mark; Dines, Jennifer N.; et
  al. “A large-scale database of T-cell receptor beta (TCRB) sequences
  and binding associations from natural and synthetic exposure to
  SARS-CoV-2.” Research Square. 2020 Aug 4:rs.3.rs-51964.
  <https://doi.org/10.21203/rs.3.rs-51964/v1>
- Schulz, Klaus U.; and Mihov, Stoyan. “Fast string correction with
  Levenshtein automata.” International Journal on Document Analysis and
  Recognition 5, no. 1 (2002): 67-85.
  <https://doi.org/10.1007/s10032-002-0082-8>
- Li, Heng; and Homer, Nils. “A survey of sequence alignment algorithms
  for next-generation sequencing.” Briefings in Bioinformatics 11, no. 5
  (2010): 473-483. <https://doi.org/10.1093/bib/bbq015>
- Hanov, Steve. “Fast and Easy Levenshtein distance using a Trie.” 2011.
  <https://stevehanov.ca/blog/index.php?id=114>
- Zorita, Eduard; Cuscó, Pol; and Filion, Guillaume J. “Starcode:
  sequence clustering based on all-pairs search.” Bioinformatics 31, no.
  12 (2015): 1913-1919. <https://doi.org/10.1093/bioinformatics/btv053>
- Buse-Dragomir, Alexandru; Popescu, Paul Stefan; and Mihaescu, Marian
  Cristian. “Spell Checker Application Based on Levenshtein Automaton.”
  In Intelligent Data Engineering and Automated Learning - IDEAL 2021,
  Lecture Notes in Computer Science 13113, 45-53. Springer, Cham, 2021.
  <https://doi.org/10.1007/978-3-030-91608-4_5>
