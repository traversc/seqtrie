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

#include "ankerl/unordered_dense.h"

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

// basic pairchar type aligned with RadixMap
using pairchar_type     = seqtrie::RadixMap::pairchar_type;
using pairchar_map_type = ankerl::unordered_dense::map<pairchar_type, int>;
using cspan             = nonstd::span<const char>;

// alias constants from RadixMap
static constexpr char GAP_CHAR      = seqtrie::RadixMap::GAP_CHAR;
static constexpr char GAP_OPEN_CHAR = seqtrie::RadixMap::GAP_OPEN_CHAR;
static constexpr char GAP_EXTN_CHAR = seqtrie::RadixMap::GAP_EXTN_CHAR;

// char counter map type using ankerl
using CharCounter      = ankerl::unordered_dense::map<char, size_t>;
using CharCounterXPtr  = Rcpp::XPtr<CharCounter>;

namespace SeqTrie {
#if USE_SEQTRIE_SMALL_ARRAY_SIZE > 0
  template <typename T>
  using array_r = trqwe::small_array<T, std::allocator<T>, size_t,
                      std::integral_constant<size_t, USE_SEQTRIE_SMALL_ARRAY_SIZE>>;
#else
  template <typename T>
  using array_r = std::vector<T>;
#endif

  using RadixTreeR   = seqtrie::RadixMap;
  // Note: ankerl::unordered_dense::map does not guarantee pointer stability on insert/delete;
  // therefore RadixForest must use std::unordered_map for stable addresses.
  using RadixForestR = std::unordered_map<size_t, RadixTreeR>;
  using node_type    = RadixTreeR::self_type;
  using pointer_type = RadixTreeR::pointer_type;
  using search_context = RadixTreeR::search_context;
  using path           = RadixTreeR::path;

  static constexpr size_t nullidx = RadixTreeR::nullidx;
  static constexpr size_t posidx  = 1;
}

using RadixTreeRXPtr   = Rcpp::XPtr<SeqTrie::RadixTreeR>;
using RadixForestRXPtr = Rcpp::XPtr<SeqTrie::RadixForestR>;
using CharCounterXPtr  = Rcpp::XPtr<CharCounter>;

// parallel-for helper

template <typename Func>
struct DoParallelFor : public RcppParallel::Worker {
  Func f;
  DoParallelFor(Func f) : f(f) {}
  void operator()(std::size_t begin, std::size_t end) { f(begin, end); }
};

template <typename Func>
inline void do_parallel_for(Func f, std::size_t begin, std::size_t end,
                            std::size_t grainSize = 1, int numThreads = -1) {
  DoParallelFor<Func> w(f);
  parallelFor(begin, end, w, grainSize, numThreads);
}

// Convert a string to a SEXP
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
  for(size_t i = 0; i < n; ++i) {
    out[i] = charsxp_to_cspan(xp[i]);
  }
  return out;
}

// convert cost matrix to map
inline pairchar_map_type convert_cost_matrix(IntegerMatrix cost_matrix) {
  pairchar_map_type cost_map;
  std::vector<char> map_elements;
  List dimnames = cost_matrix.attr("dimnames");
  CharacterVector rownames = dimnames[0];
  map_elements.resize(rownames.size());
  for(size_t i = 0; i < rownames.size(); ++i) {
    if(rownames[i] == "gap") {
      map_elements[i] = GAP_CHAR;
    } else if(rownames[i] == "gap_open") {
      map_elements[i] = GAP_OPEN_CHAR;
    } else {
      Rcpp::String s = rownames[i];
      map_elements[i] = s.get_cstring()[0];
    }
  }
  size_t N = map_elements.size();
  for(size_t i = 0; i < N; ++i) {
    for(size_t j = 0; j < N; ++j) {
      if((map_elements[i] == GAP_CHAR || map_elements[i] == GAP_OPEN_CHAR) &&
         (map_elements[j] == GAP_CHAR || map_elements[j] == GAP_OPEN_CHAR))
        continue;
      cost_map[{map_elements[i], map_elements[j]}] = cost_matrix(i, j);
    }
  }
  return cost_map;
}

// convert search results to DataFrame
inline DataFrame seqtrie_results_to_dataframe(CharacterVector query,
                                              std::vector<SeqTrie::search_context> & output) {
  size_t nresults = 0;
  for(auto &ctx : output) nresults += ctx.match.size();
  CharacterVector query_results(nresults);
  CharacterVector target_results(nresults);
  IntegerVector distance_results(nresults);
  int *dist_ptr = INTEGER(distance_results);
  size_t q = 0;
  for(size_t i = 0; i < output.size(); ++i) {
    auto &ctx = output[i];
    for(size_t j = 0; j < ctx.match.size(); ++j) {
      SET_STRING_ELT(query_results, q, STRING_ELT(query, i));
      auto s = ctx.match[j]->template sequence<SeqTrie::array_r<char>>();
      SET_STRING_ELT(target_results, q, to_charsxp(s));
      dist_ptr[q++] = ctx.distance[j];
    }
  }
  return DataFrame::create(_["query"]    = query_results,
                           _["target"]   = target_results,
                           _["distance"] = distance_results,
                           _["stringsAsFactors"] = false);
}

#endif // seqtrie_TYPES_H
