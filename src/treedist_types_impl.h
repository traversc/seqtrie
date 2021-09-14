#ifndef TREEDIST_TYPES_IMPL_H
#define TREEDIST_TYPES_IMPL_H


#include <Rcpp.h>

#include <cstring>
#include <utility>

#include <treedist/prefixmap.h>
#include <treedist/radixarray.h>
#include <treedist/radixmap.h>

#include "treedist_types.h"

using namespace Rcpp;

// using DNATree = rtree<treedist::RadixArray<4>>;
// using RadixTree = rtree<treedist::RadixMap<>>;
// using PrefixTree = rtree<treedist::PrefixMap<>>;

template <typename T> struct dispatch {};
template<> struct dispatch<RadixTree::value_type>{
  static inline uspan sequence_convert(const cspan sequence) { 
    return uspan{reinterpret_cast<const uint8_t*>(sequence.data()), sequence.size()};
  };
};
template<> struct dispatch<DNATree::value_type>{
  static inline trqwe::simple_array<uint8_t> sequence_convert(const cspan sequence) {
    trqwe::simple_array<uint8_t> result(sequence.size());
    for(size_t i=0; i<sequence.size(); ++i) {
      switch(sequence[i]) {
      case 'A':
        result[i] = 0;
        break;
      case 'C':
        result[i] = 1;
        break;
      case 'G':
        result[i] = 2;
        break;
      case 'T':
        result[i] = 3;
        break;
      default:
        throw std::runtime_error("sequence must have only A, C, G or T");
      }
    }
    return result;
  }
};

template<> struct dispatch<PrefixTree::value_type>{
  static inline uspan sequence_convert(const cspan sequence) { 
    return uspan{reinterpret_cast<const uint8_t*>(sequence.data()), sequence.size()};
  };
};

static std::string cspan_to_string(const cspan x) {
  return std::string(x.data(), x.size());
}

inline SEXP cspan_to_charsxp(const cspan x) {
  return Rf_mkCharLen(x.data(), x.size());
}

template <typename T> rtree<T>::rtree() : root(new value_type) {}

template <typename T> typename rtree<T>::index_type rtree<T>::size() const {
  index_type result = 0;
  for(size_t i=0; i<sequence_map.size(); ++i) {
    if(!sequence_map[i].is_null()) result++;
  }
  return result;
}

template <typename T> typename rtree<T>::index_type rtree<T>::insert(const cspan sequence) {
  index_type next_idx = static_cast<index_type>(sequence_map.size());
  auto usequence = dispatch<typename T::value_type>::sequence_convert(sequence);
  index_type idx = value_type::insert(root, usequence, next_idx);
  // std::cout << "insert " <<  (void*)sequence.data() << " " << sequence.size() << " " << idx << std::endl; 
  if(idx == value_type::nullidx) {
    sequence_map.push_back(ncstring(sequence.data(), sequence.size()));
    // std::cout << "insert " << (void*)sequence_map[sequence_map.size() - 1].data() << " " << sequence_map[sequence_map.size() - 1].size() << std::endl;
  }
  return idx;
}

template <typename T> typename rtree<T>::index_type rtree<T>::erase(const cspan sequence) {
  auto usequence = dispatch<typename T::value_type>::sequence_convert(sequence);
  index_type idx = value_type::erase(root, usequence);
  if(idx != value_type::nullidx) {
    sequence_map[idx].nullify();
  }
  return idx;
}

template <typename T> typename rtree<T>::index_type rtree<T>::find(const cspan sequence) const {
  auto usequence = dispatch<typename T::value_type>::sequence_convert(sequence);
  return value_type::find(root, usequence);
}

template <typename T> auto rtree<T>::levenshtein(const cspan sequence, const int max_distance) const {
  auto usequence = dispatch<typename T::value_type>::sequence_convert(sequence);
  return levenshtein_type{root, usequence, max_distance}.search();
}

template <typename T> auto rtree<T>::hamming(const cspan sequence, const int max_distance) const {
  auto usequence = dispatch<typename T::value_type>::sequence_convert(sequence);
  return hamming_type{root, usequence, max_distance}.search();
}

template <typename T> cspan rtree<T>::index_to_sequence(const index_type idx) const {
  // std::cout << (void*)sequence_map[idx].data() << " " << sequence_map[idx].size() << std::endl;
  return cspan(sequence_map[idx].data(), sequence_map[idx].size()); // must not be nullptr
}

template <typename T> std::string rtree<T>::print() const {
  return value_type::print(root, 0);
}
template <> std::string DNATree::print() const {
  return value_type::print(root, 0, "ACGT");
}

template <typename T> SEXP rtree<T>::to_dataframe() const {
  size_t n = this->size();
  if(n == 0) return R_NilValue;
  CharacterVector sequence(n);
  NumericVector index(n);
  double * indexp = REAL(index);
  size_t q = 0;
  for(size_t i=0; i<sequence_map.size(); ++i) {
    if(!sequence_map[i].is_null()) {
      SET_STRING_ELT(sequence, q, Rf_mkCharLen(sequence_map[i].data(), sequence_map[i].size()));
      indexp[q] = static_cast<double>(i);
    }
  }
  return DataFrame::create(_["sequence"] = sequence, _["index"] = index, _["stringsAsFactors"] = false);
}

#endif // include guard
