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

#if __cplusplus >= 201703L
#define USE_ANKERL 1
#include "ankerl/unordered_dense.h"
#else
#include <boost/functional/hash.hpp>
#endif

#include "seqtrie/radixmap.h"
#include "simple_array/small_array.h"

// defined in some headers in windows and Mac, conflicts with R headers
#ifdef TRUE
#undef TRUE
#endif
#ifdef FALSE
#undef FALSE
#endif

#if (R_VERSION < R_Version(3, 5, 0))
#define STRING_PTR_RO STRING_PTR
#endif


#define USE_SEQTRIE_SMALL_ARRAY_SIZE SEQTRIE_SMALL_ARRAY_SIZE

using namespace Rcpp;
using namespace RcppParallel;

// constants and types
using pairchar_type = std::pair<char, char>;
#ifdef USE_ANKERL
using pairchar_map_type = ankerl::unordered_dense::map<pairchar_type, int>;
#else
using pairchar_map_type = std::unordered_map<pairchar_type, int, boost::hash<pairchar_type>>;
#endif
using cspan = nonstd::span<const char>;
constexpr char GAP_CHAR = '\0';                                     // '\0' any gap cost for non-affine
constexpr char GAP_OPEN_CHAR = std::numeric_limits<char>::min();    // '\255' gap open cost for affine
constexpr char GAP_EXTN_CHAR = '\0';    

// used in utils.cpp, a map for counting chars, to make sure input cost_matrix contains all chars in a trie
#ifdef USE_ANKERL
using CharCounter = ankerl::unordered_dense::map<char, size_t>;
#else
using CharCounter = std::unordered_map<char, size_t>;
#endif

using CharCounterXPtr = Rcpp::XPtr<CharCounter>;

// inline
// Convert a string to a SEXP
// Be careful about R protection / GC
// SEXP to_charsxp(const SeqTrie::array_r<char> & x);

// inline
// Define a span of const char from a SEXP
// cspan charsxp_to_cspan(SEXP x);

// inline
// CharacterVector to cspan vector
// std::vector<cspan> strsxp_to_cspan(CharacterVector x);

// Input: cost_matrix
// a NxN matrix where column/row names are the characters to use for pairchar_map_type keys
// The special column "gap" is recoded as '\0'
// Output: pairchar_map_type
// pairchar_map_type convert_cost_matrix(IntegerMatrix cost_matrix);

namespace SeqTrie {
  // If using ankerl, unordered_dense map does not have stable addresses on insert/delete
  // due to data being stored in a plain std::vector. Cannot use for RadixForest.
  // Not necessary for RadixTree as all nodes use unique_ptr
  // To do: convert RadixForest to use unique_ptr
#if USE_SEQTRIE_SMALL_ARRAY_SIZE > 0
  template <typename T> using array_r = trqwe::small_array<T, std::allocator<T>, size_t, std::integral_constant<size_t, USE_SEQTRIE_SMALL_ARRAY_SIZE>>;
#else
  template <typename T> using array_r = std::vector<T>;
#endif
#ifdef USE_ANKERL
  using RadixTreeR = seqtrie::RadixMap<char, ankerl::unordered_dense::map, array_r, size_t>;
#else
  using RadixTreeR = seqtrie::RadixMap<char, std::map, array_r, size_t>;
#endif
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

// parallel for helper function
// accepts functors (lambdas) with std::size_t begin, std::size_t end
template <typename Func>
struct DoParallelFor : public RcppParallel::Worker {
  Func f;
  DoParallelFor(Func f) : f(f) {}
  void operator()(std::size_t begin, std::size_t end) {
    f(begin, end);
  }
  
};
template <typename Func>
inline void do_parallel_for(Func f, std::size_t begin, std::size_t end, std::size_t grainSize = 1, int numThreads = -1) {
  DoParallelFor<Func> w(f);
  parallelFor(begin, end, w, grainSize, numThreads);
}

// Convert a string to a SEXP
// Be careful about R protection / GC
inline SEXP to_charsxp(const SeqTrie::array_r<char> & x) {
  return Rf_mkCharLen(x.data(), x.size());
}

// Define a span of const char from a SEXP
inline cspan charsxp_to_cspan(SEXP x) {
  return cspan(CHAR(x), Rf_xlength(x));
}

inline std::vector<cspan> strsxp_to_cspan(CharacterVector x) {
  size_t n = Rf_xlength(x);
  const SEXP * xp = STRING_PTR_RO(x);
  std::vector<cspan> out(n);
  for(size_t i=0; i<n; ++i) {
    out[i] = charsxp_to_cspan(xp[i]);
  }
  return out;
}

// Input: cost_matrix
// a NxN matrix where column/row names are the characters to use for pairchar_map_type keys
// The special column "gap" is recoded as '\0'
// Output: pairchar_map_type
inline pairchar_map_type convert_cost_matrix(IntegerMatrix cost_matrix) {
  pairchar_map_type cost_map;
  std::vector<char> map_elements;
  {
    List dimnames = cost_matrix.attr("dimnames");
    CharacterVector rownames = dimnames[0];
    map_elements.resize(rownames.size());
    for(size_t i=0; i<map_elements.size(); ++i) {
      if(rownames[i] == "gap") {
        map_elements[i] = GAP_CHAR; // '\0' same symbol as GAP_EXTN_CHAR
      } else if(rownames[i] == "gap_open") {
        map_elements[i] = GAP_OPEN_CHAR; // '\255'
      } else {
        Rcpp::String s = rownames[i];
        map_elements[i] = s.get_cstring()[0];
      }
    }
  }
  size_t N = map_elements.size();
  for(size_t i=0; i<N; ++i) {
    for(size_t j=0; j<N; ++j) {
      // skip any combination of gap_open, gap
      if((map_elements[i] == GAP_CHAR || map_elements[i] == GAP_OPEN_CHAR) &&
         (map_elements[j] == GAP_CHAR || map_elements[j] == GAP_OPEN_CHAR)) continue;
      cost_map[pairchar_type(map_elements[i], map_elements[j])] = cost_matrix(i,j);
    }
  }
  return cost_map;
}

// converts seqtrie search results to a data.frame
inline DataFrame seqtrie_results_to_dataframe(CharacterVector query, std::vector<SeqTrie::search_context> & output) {
  size_t nresults = 0;
  size_t nseqs = output.size();
  for(size_t i=0; i<nseqs; ++i) { nresults += output[i].match.size(); }
  CharacterVector query_results(nresults);
  CharacterVector target_results(nresults);
  IntegerVector distance_results(nresults);
  int * distance_results_ptr = INTEGER(distance_results);
  size_t q = 0;
  for(size_t i=0; i<nseqs; ++i) {
    auto & targets = output[i].match;
    auto & distances = output[i].distance;
    for(size_t j=0; j<targets.size(); ++j) {
      SET_STRING_ELT(query_results, q, STRING_ELT(query, i));
      auto s = targets[j]->template sequence<SeqTrie::array_r<char>>();
      SET_STRING_ELT(target_results, q, to_charsxp(s));
      distance_results_ptr[q] = distances[j];
      q++;
    }
  }
  return DataFrame::create(_["query"] = query_results,
                           _["target"] = target_results,
                           _["distance"] = distance_results,
                           _["stringsAsFactors"] = false);
}

#endif
