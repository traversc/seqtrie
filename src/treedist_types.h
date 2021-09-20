#ifndef TREEDIST_TYPES_H
#define TREEDIST_TYPES_H

#include <Rcpp.h>

#include <treedist/prefixmap.h>
#include <treedist/radixarray.h>
#include <treedist/radixmap.h>

using namespace Rcpp;

using cspan = treedist::cspan;
using uspan = treedist::uspan;

template <typename T> class rtree;
using DNATree = rtree<treedist::RadixArray<boost::mpl::string<'A','C','G','T'>, trqwe::small_array, size_t>>;
using RadixTree = rtree<treedist::RadixMap<std::unordered_map, trqwe::small_array, size_t>>;
using PrefixTree = rtree<treedist::PrefixMap<std::unordered_map, size_t>>;

template <typename T> class rtree {
public:
  typedef typename T::value_type value_type;
  typedef typename T::pointer_type pointer_type;
  typedef typename T::index_type index_type;
  typedef typename T::Levenshtein levenshtein_type;
  typedef typename T::Hamming hamming_type;
  typedef trqwe::nullable_array<char> ncstring; // nullable string to store sequences; NOT null terminated
  typedef std::vector<ncstring> seqmap_type;
  typedef typename seqmap_type::value_type seqmap_value_type;
private:
  pointer_type root;
  seqmap_type sequence_map;
public:
  rtree();
  index_type size() const;
  NumericVector insert(CharacterVector sequences);
  NumericVector erase(CharacterVector sequences);
  NumericVector find(CharacterVector sequences) const;
  SEXP find_prefix(CharacterVector sequences) const;
  SEXP levenshtein_search(CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress) const;
  SEXP hamming_search(CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress) const;
  SEXP prefix_search(CharacterVector sequences, const int nthreads, const bool show_progress) const;
  cspan index_to_sequence(const index_type idx) const;
  std::string print() const;
  SEXP to_dataframe() const;
};

#endif // include guard