#ifndef TREEDIST_TYPES_H
#define TREEDIST_TYPES_H

#include <Rcpp.h>

#include <treedist/prefixmap.h>
#include <treedist/radixarray.h>
#include <treedist/radixmap.h>

// #include <boost/bimap.hpp>
// #include <boost/bimap/unordered_set_of.hpp>

using cspan = treedist::cspan;
using uspan = treedist::uspan;

template <typename T> class rtree;
using DNATree = rtree<treedist::RadixArray<4>>;
using RadixTree = rtree<treedist::RadixMap<>>;
using PrefixTree = rtree<treedist::PrefixMap<>>;

template <typename T> class rtree {
public:
  typedef typename T::value_type value_type;
  typedef typename T::pointer_type pointer_type;
  typedef typename T::index_type index_type;
  typedef typename T::template Levenshtein<std::vector> levenshtein_type;
  typedef typename T::template Hamming<std::vector> hamming_type;
  typedef trqwe::nullable_array<char> ncstring; // nullable string to store sequences; NOT null terminated
  typedef std::vector<ncstring> seqmap_type;
  typedef typename seqmap_type::value_type seqmap_value_type;
private:
  pointer_type root;
  seqmap_type sequence_map;
public:
  rtree();
  index_type size() const;
  index_type insert(const cspan sequence);
  index_type erase(const cspan sequence);
  index_type find(const cspan sequence) const;
  auto levenshtein(const cspan sequence, const int max_distance) const;
  auto hamming(const cspan sequence, const int max_distance) const;
  cspan index_to_sequence(const index_type idx) const;
  std::string print() const;
  SEXP to_dataframe() const;
};

#endif // include guard