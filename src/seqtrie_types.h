#ifndef seqtrie_TYPES_H
#define seqtrie_TYPES_H

#include <Rcpp.h>

#include <unordered_map>
#include <map>
#include <cstring>
#include <utility>

#include "seqtrie/radixmap.h"
#include "simple_array/small_array.h"

using namespace Rcpp;

using pairchar = std::pair<char, char>;
using pairchar_map = std::unordered_map<pairchar, int>;
using cspan = nonstd::span<const char>;

// Convert a string to a SEXP
// Be careful about R protection / GC
SEXP to_charsxp(const trqwe::small_array<char> & x);

// Define a span of const char from a SEXP
cspan charsxp_to_cspan(SEXP x);

// Input: cost_matrix
// a NxN matrix where column/row names are the characters to use for pairchar_map keys
// The special column "gap" is recoded as '\0'
// Output: pairchar_map
pairchar_map convert_cost_matrix(IntegerMatrix cost_matrix);

namespace SeqTrie {
  using RadixTreeR = seqtrie::RadixMap<char, std::map, trqwe::small_array, size_t>;
  using RadixForestR = std::unordered_map<size_t, RadixTreeR>;
  typedef typename RadixTreeR::self_type node_type;
  typedef typename RadixTreeR::pointer_type pointer_type;
  typedef typename RadixTreeR::search_context search_context;
  typedef typename RadixTreeR::path path;
  static constexpr size_t nullidx = RadixTreeR::nullidx; // the sequence does not exist in tree
  static constexpr size_t posidx = 1; // the sequence exists in tree
};

using RadixTreeRXPtr = Rcpp::XPtr<SeqTrie::RadixTreeR>;
using RadixForestRXPtr = Rcpp::XPtr<SeqTrie::RadixForestR>;

#endif

// double size() const;
// LogicalVector insert(CharacterVector sequences);
// LogicalVector erase(CharacterVector sequences);
// LogicalVector find(CharacterVector sequences) const;
// SEXP levenshtein_search(CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress) const;
// SEXP hamming_search(CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress) const;
// SEXP anchored_search(CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress) const;
// SEXP prefix_search(CharacterVector sequences) const;
// SEXP search(CharacterVector sequences, IntegerVector max_distance, IntegerMatrix substitution_distance, const int gap_distance, const int nthreads, const bool show_progress) const;
// std::string print() const;
// SEXP graph(const double depth) const;
// SEXP to_vector() const;
// bool validate() const;