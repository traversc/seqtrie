#ifndef seqtrie_TYPES_H
#define seqtrie_TYPES_H

#include <Rcpp.h>
#include <RcppParallel.h>

#include <unordered_map>
#include <map>
#include <cstring>
#include <utility>
#include <set>
#include <memory>
#include <tuple>

#include "seqtrie/radixmap.h"
#include "simple_array/small_array.h"
#include "simple_progress/simple_progress.h"

using namespace Rcpp;
using namespace RcppParallel;

// constants and types
using pairchar = std::pair<char, char>;
using pairchar_map = std::unordered_map<pairchar, int>;
using cspan = nonstd::span<const char>;

// defined in utils.cpp
// Convert a string to a SEXP
// Be careful about R protection / GC
SEXP to_charsxp(const trqwe::small_array<char> & x);

// defined in utils.cpp
// Define a span of const char from a SEXP
cspan charsxp_to_cspan(SEXP x);

// defined in utils.cpp
// Input: cost_matrix
// a NxN matrix where column/row names are the characters to use for pairchar_map keys
// The special column "gap" is recoded as '\0'
// Output: pairchar_map
pairchar_map convert_cost_matrix(IntegerMatrix cost_matrix);

// defined in utils.cpp
// counts chars for strings
// keys returns the chars that are non-zero as a CharacterVector
struct CharCounter {
  std::unordered_map<char, size_t> counts;
  void add(cspan s);
  void subtract(cspan s);
  CharacterVector keys(); 
};

namespace SeqTrie {
  using RadixTreeR = seqtrie::RadixMap<char, std::map, trqwe::small_array, size_t>;
  using RadixForestR = std::unordered_map<size_t, RadixTreeR>;
  typedef typename RadixTreeR::self_type node_type;
  typedef typename RadixTreeR::pointer_type pointer_type;
  typedef typename RadixTreeR::search_context search_context;
  typedef typename RadixTreeR::path path;
  static constexpr size_t nullidx = RadixTreeR::nullidx; // the sequence does not exist in tree
  static constexpr size_t posidx = 1; // the sequence exists in tree
}
using RadixTreeRXPtr = Rcpp::XPtr<SeqTrie::RadixTreeR>;
using RadixForestRXPtr = Rcpp::XPtr<SeqTrie::RadixForestR>;
using CharCounterXPtr = Rcpp::XPtr<CharCounter>;


#endif
