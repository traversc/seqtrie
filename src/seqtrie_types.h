#ifndef seqtrie_TYPES_H
#define seqtrie_TYPES_H

#include <Rcpp.h>
#include <RcppParallel.h>

#include <cstring>
#include <utility>
#include "seqtrie/radixmap.h"
#include "simple_array/small_array.h"

using namespace Rcpp;
using namespace RcppParallel;

// design approach:
// this is a realization of Radix Tree for R, so whenever possible use explicit type labels
// rather than dependent type names

using cspan = nonstd::span<const char>;

template <typename T> class RadixTreeContext;
using RadixTreeCtxForR = RadixTreeContext<seqtrie::RadixMap<char, std::map, trqwe::small_array, size_t>>;

template <typename T> class RadixTreeContext {
public:
  typedef typename T::self_type node_type;
  typedef typename T::pointer_type pointer_type;
  typedef typename T::search_context search_context;
  typedef typename T::path path;
  static constexpr size_t nullidx = T::nullidx; // the sequence does not exist in tree
  static constexpr size_t posidx = 1; // the sequence exists in tree
private:
  node_type root;
public:
  RadixTreeContext() {}
  double size() const;
  LogicalVector insert(CharacterVector sequences);
  LogicalVector erase(CharacterVector sequences);
  LogicalVector find(CharacterVector sequences) const;
  SEXP levenshtein_search(CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress) const;
  SEXP hamming_search(CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress) const;
  SEXP anchored_search(CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress) const;
  SEXP prefix_search(CharacterVector sequences) const;
  std::string print() const;
  SEXP graph(const double depth) const;
  SEXP to_vector() const;
  bool validate() const;
};

#endif // include guard