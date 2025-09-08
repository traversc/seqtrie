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
#include "seqtrie/utility.h"
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

// basic types and spans
using pairchar_type = seqtrie::pairchar_type;
using cspan        = nonstd::span<const char>;
using CostMap      = seqtrie::CostMap;

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

// ===== Centralized alignment algorithm selection =====
enum class AlignmentAlgo {
  Hamming,
  GlobalUnit,
  AnchoredUnit,
  GlobalLinear,
  AnchoredLinear,
  GlobalAffine,
  AnchoredAffine
};

inline bool is_unit_substitution(const Rcpp::Nullable<Rcpp::IntegerMatrix>& cm) {
  if (cm.isNull()) return true; // treat NULL as unit substitution (match=0, mismatch=1)
  Rcpp::IntegerMatrix m = cm.get();
  if (m.nrow() != m.ncol()) return false;
  for (int i = 0; i < m.nrow(); ++i) {
    for (int j = 0; j < m.ncol(); ++j) {
      const int expected = (i == j) ? 0 : 1;
      if (m(i, j) != expected) return false;
    }
  }
  return true;
}

inline AlignmentAlgo decide_alignment_algo(std::string mode,
                                           const Rcpp::Nullable<Rcpp::IntegerMatrix>& cost_matrix,
                                           int gap_cost,
                                           int gap_open_cost) {
  // normalize mode
  for (auto& c : mode) c = static_cast<char>(std::tolower(c));
  if (mode == "hm") mode = "hamming";
  else if (mode == "gb" || mode == "lv" || mode == "levenshtein") mode = "global";
  else if (mode == "an" || mode == "en" || mode == "extension") mode = "anchored";

  if (mode == "hamming") return AlignmentAlgo::Hamming;

  // If substitution matrix is NULL, always use unit substitution and unit gaps
  // ignoring gap_cost and gap_open_cost entirely, per refined logic.
  if (cost_matrix.isNull()) {
    return (mode == "global") ? AlignmentAlgo::GlobalUnit : AlignmentAlgo::AnchoredUnit;
  }

  const bool affine = gap_open_cost > 0;      // non-zero gap_open means affine
  const bool unit_subs = is_unit_substitution(cost_matrix);
  const bool unit_gaps = (gap_cost == 1) && !affine; // unit gaps only when linear with gap=1

  // Treat boolean (0/1) matrices with unit gaps as Unit, so callers can route to Myers
  bool boolean_subs = unit_subs;
  if (!unit_subs) {
    // Check explicitly for 0/1 matrices (not necessarily identity)
    Rcpp::IntegerMatrix m = cost_matrix.get();
    if (m.nrow() == m.ncol()) {
      boolean_subs = true;
      for (int i = 0; i < m.nrow() && boolean_subs; ++i) {
        for (int j = 0; j < m.ncol(); ++j) {
          int v = m(i, j);
          if (v != 0 && v != 1) { boolean_subs = false; break; }
        }
      }
    }
  }

  if (mode == "global") {
    if (affine) return AlignmentAlgo::GlobalAffine;
    if (unit_gaps && (unit_subs || boolean_subs)) return AlignmentAlgo::GlobalUnit;
    return AlignmentAlgo::GlobalLinear;
  } else { // anchored
    if (affine) return AlignmentAlgo::AnchoredAffine;
    if (unit_gaps && (unit_subs || boolean_subs)) return AlignmentAlgo::AnchoredUnit;
    return AlignmentAlgo::AnchoredLinear;
  }
}

// convert cost matrix to map
inline CostMap convert_cost_matrix(IntegerMatrix cost_matrix, int gap_cost, int gap_open_cost) {
  CostMap cm;
  std::vector<char> map_elements;
  List dimnames = cost_matrix.attr("dimnames");
  CharacterVector rownames = dimnames[0];
  map_elements.resize(rownames.size());
  for(size_t i = 0; i < rownames.size(); ++i) {
    if(rownames[i] == "gap" || rownames[i] == "gap_open") {
      // special tokens are ignored for substitution table
      map_elements[i] = '\0';
    } else {
      Rcpp::String s = rownames[i];
      map_elements[i] = s.get_cstring()[0];
    }
  }
  // Collect alphabet (exclude special rows/cols)
  std::vector<char> alphabet;
  for(size_t i = 0; i < rownames.size(); ++i) {
    if(rownames[i] != "gap" && rownames[i] != "gap_open") {
      alphabet.push_back(map_elements[i]);
    }
  }
  // Fill substitution cost map using original indices for non-special rows/cols
  std::vector<int> alpha_idx;
  alpha_idx.reserve(alphabet.size());
  for(int i = 0; i < rownames.size(); ++i) {
    if(rownames[i] != "gap" && rownames[i] != "gap_open") alpha_idx.push_back(i);
  }
  for(size_t ai = 0; ai < alpha_idx.size(); ++ai) {
    for(size_t aj = 0; aj < alpha_idx.size(); ++aj) {
      char ci = map_elements[alpha_idx[ai]];
      char cj = map_elements[alpha_idx[aj]];
      cm.char_cost_map[{ci, cj}] = cost_matrix(alpha_idx[ai], alpha_idx[aj]);
    }
  }
  cm.gap_cost = gap_cost;
  cm.gap_open_cost = gap_open_cost;
  return cm;
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
