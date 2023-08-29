#ifndef seqtrie_TYPES_H
#define seqtrie_TYPES_H

#include <Rcpp.h>

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

// constants and types
using pairchar_type = std::pair<char, char>;
using pairchar_map_type = std::unordered_map<pairchar_type, int>;
using cspan = nonstd::span<const char>;
constexpr char GAP_CHAR = '\0';                                     // '\0' any gap cost for non-affine
constexpr char GAP_OPEN_CHAR = std::numeric_limits<char>::min();    // '\255' gap open cost for affine
constexpr char GAP_EXTN_CHAR = '\0';    

// used in utils.cpp, a map for counting chars, to make sure input cost_matrix contains all chars in a trie
using CharCounter = std::unordered_map<char, size_t>;
using CharCounterXPtr = Rcpp::XPtr<CharCounter>;

// defined in utils.cpp
// Convert a string to a SEXP
// Be careful about R protection / GC
SEXP to_charsxp(const trqwe::small_array<char> & x);

// defined in utils.cpp
// Define a span of const char from a SEXP
cspan charsxp_to_cspan(SEXP x);

// defined in utils.cpp
// CharacterVector to cspan vector
std::vector<cspan> strsxp_to_cspan(CharacterVector x);

// defined in utils.cpp
// Input: cost_matrix
// a NxN matrix where column/row names are the characters to use for pairchar_map_type keys
// The special column "gap" is recoded as '\0'
// Output: pairchar_map_type
pairchar_map_type convert_cost_matrix(IntegerMatrix cost_matrix);

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
