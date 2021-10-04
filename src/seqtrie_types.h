#ifndef seqtrie_TYPES_H
#define seqtrie_TYPES_H

#include <Rcpp.h>

// #include <seqtrie/prefixmap.h>
// #include <seqtrie/radixarray.h>
#include <seqtrie/radixmap.h>

using namespace Rcpp;

// design approach:
// this is a realization of Radix Tree for R, so whenever possible use explicit type labels
// rather than dependent type names

using cspan = nonstd::span<const char>;

template <typename T> class rtree;
using RadixTree = rtree<seqtrie::RadixMap<char, std::map, trqwe::small_array, size_t>>;

template <typename T> class rtree {
public:
  typedef typename T::self_type node_type;
  typedef typename T::pointer_type pointer_type;
  typedef typename T::search_context search_context;
  typedef typename T::path path;
  // typedef trqwe::nullable_array<char> ncstring; // nullable string to store sequences; NOT null terminated
  // typedef std::vector<ncstring> seqmap_type;
  static constexpr size_t nullidx = T::nullidx;
private:
  node_type root;
  size_t current_idx;
  // seqmap_type sequence_map;
public:
  rtree() : current_idx(0) {}
  double size() const;
  NumericVector insert(CharacterVector sequences);
  NumericVector erase(CharacterVector sequences);
  NumericVector find(CharacterVector sequences) const;
  SEXP levenshtein_search(CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress) const;
  SEXP hamming_search(CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress) const;
  SEXP prefix_search(CharacterVector sequences) const;
  std::string print() const;
  SEXP graph(const double depth) const;
  SEXP to_dataframe() const;
  bool validate() const;
  // cspan index_to_sequence(const size_t idx) const;
};

#endif // include guard