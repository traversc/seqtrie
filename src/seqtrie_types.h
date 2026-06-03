#ifndef seqtrie_TYPES_H
#define seqtrie_TYPES_H

#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>
#include <array>
#include <cstring>
#include <utility>
#include <memory>
#include <tuple>
#include <limits>
#include <stdexcept>
#include <string>

#include "ankerl/unordered_dense.h"

#include "seqtrie/radix_forest.h"
#include "seqtrie/utility.h"
#include "simple_array/small_array.h"
#include "startree/startree.h"

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

// SEQTRIE_SMALL_ARRAY_SIZE default is provided by seqtrie/radixmap.h (included
// above), so it is defined here. 0 selects the std::vector fallback below.
#define USE_SEQTRIE_SMALL_ARRAY_SIZE SEQTRIE_SMALL_ARRAY_SIZE

using namespace Rcpp;
using namespace RcppParallel;

// basic types and spans
using pairchar_type = seqtrie::pairchar_type;
using cspan        = nonstd::span<const char>;
using CostMap      = seqtrie::CostMap;
using QueryCostCache = seqtrie::QueryCostCache;

// Fixed byte-alphabet character counter.  Character values are indexed as
// unsigned bytes so all possible char values, including negative signed-char
// representations, map safely into 0..255.
struct CharCounter {
  static constexpr size_t alphabet_size = 256;
  std::array<size_t, alphabet_size> counts{};

  static size_t index(char c) noexcept {
    return static_cast<size_t>(static_cast<unsigned char>(c));
  }

  void add(char c) noexcept {
    ++counts[index(c)];
  }

  void subtract(char c) noexcept {
    size_t & value = counts[index(c)];
    if(value > 0) --value;
  }

  size_t count(char c) const noexcept {
    return counts[index(c)];
  }

  size_t size() const noexcept {
    size_t n = 0;
    for(size_t value : counts) {
      if(value != 0) ++n;
    }
    return n;
  }
};
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

  using RadixMapCore = seqtrie::RadixMap;
  using RadixTreeR = seqtrie::RadixTree;
  using RadixForestR = seqtrie::RadixForest;
  using node_type = RadixMapCore;
  using pointer_type = RadixMapCore::pointer_type;
  using search_result = RadixMapCore::search_result;
  using search_context = RadixMapCore::search_context;
  using search_context_cached = RadixMapCore::search_context_cached;
  using path = RadixMapCore::path;

  static constexpr size_t nullidx = RadixMapCore::nullidx;
  static constexpr size_t posidx  = 1;
}

struct SearchFilterOptions {
  bool lower_triangle = false;
  bool best_only = false;

  bool needs_hook() const noexcept {
    return lower_triangle || best_only;
  }
};

inline SearchFilterOptions make_search_filter_options(const bool lower_triangle,
                                                      const std::string & match_mode) {
  SearchFilterOptions options;
  options.lower_triangle = lower_triangle;
  if(match_mode == "best") {
    options.best_only = true;
  } else if(match_mode != "all") {
    throw std::runtime_error("Internal error: match_mode must be 'all' or 'best'");
  }
  return options;
}

struct SearchResultFilterHook {
  SeqTrie::RadixMapCore::index_type query_index;
  bool lower_triangle;
  bool best_only;
  int best_distance = SeqTrie::RadixMapCore::NO_ALIGN;

  SearchResultFilterHook(SeqTrie::RadixMapCore::index_type query_index_value,
                         bool lower_triangle_value,
                         bool best_only_value) noexcept
    : query_index(query_index_value),
      lower_triangle(lower_triangle_value),
      best_only(best_only_value) {}

  template <typename Context>
  SeqTrie::RadixMapCore::SearchHookAction operator()(Context & ctx,
                                                     SeqTrie::RadixMapCore::path candidate,
                                                     int distance) noexcept {
    using Action = SeqTrie::RadixMapCore::SearchHookAction;
    if(lower_triangle) {
      const auto target_index = candidate->get_terminal_idx();
      if(target_index == SeqTrie::nullidx || query_index <= target_index) {
        return Action::skip_continue;
      }
    }

    if(best_only) {
      if(distance < best_distance) {
        best_distance = distance;
        ctx.match.clear();
        ctx.distance.clear();
        return Action::add_continue;
      }
      if(distance == best_distance) {
        return Action::add_continue;
      }
      return Action::skip_continue;
    }

    return Action::add_continue;
  }
};

inline SearchResultFilterHook make_search_result_filter_hook(const SearchFilterOptions & options,
                                                             const int query_index) {
  return SearchResultFilterHook(
    static_cast<SeqTrie::RadixMapCore::index_type>(query_index),
    options.lower_triangle,
    options.best_only
  );
}

template <typename DefaultSearch, typename HookedSearch>
inline SeqTrie::search_result run_search_with_optional_hook(const SearchFilterOptions & options,
                                                            const int query_index,
                                                            DefaultSearch default_search,
                                                            HookedSearch hooked_search) {
  if(options.needs_hook()) {
    auto hook = make_search_result_filter_hook(options, query_index);
    return hooked_search(hook);
  }
  return default_search();
}

using RadixTreeRXPtr   = Rcpp::XPtr<SeqTrie::RadixTreeR>;
using RadixForestRXPtr = Rcpp::XPtr<SeqTrie::RadixForestR>;

struct StarTreeR {
  startree::InputData data;
  std::vector<std::string> sequences;
  startree::SearchParams params;
  std::vector<startree::PairRecord> self_pairs;
  int nthreads = 1;
  bool hamming = false;  // mode == "hamming": substitution-only, equal length
};
using StarTreeRXPtr = Rcpp::XPtr<StarTreeR>;

struct AnchoredStarTreeR {
  startree::AnchoredInputData data;
  std::vector<std::string> sequences;
  startree::SearchParams params;
  std::vector<startree::PairRecord> self_pairs;
  int nthreads = 1;
};
using AnchoredStarTreeRXPtr = Rcpp::XPtr<AnchoredStarTreeR>;

inline int checked_r_len(const size_t n) {
  if(n > static_cast<size_t>(std::numeric_limits<int>::max())) {
    throw std::length_error("object is too large for an R vector");
  }
  return static_cast<int>(n);
}

inline size_t checked_size_from_r_xlen(const R_xlen_t n) {
  if(n < 0) {
    throw std::length_error("negative R vector length");
  }
  return static_cast<size_t>(n);
}

inline R_xlen_t checked_r_xlen(const size_t n) {
  if(n > static_cast<size_t>(std::numeric_limits<R_xlen_t>::max())) {
    throw std::length_error("object is too large for an R vector");
  }
  return static_cast<R_xlen_t>(n);
}

// parallel-for helper

template <typename Func>
struct DoParallelFor : public RcppParallel::Worker {
  Func f;
  explicit DoParallelFor(Func func) : f(func) {}
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
  return Rf_mkCharLen(x.data(), checked_r_len(x.size()));
}

// Define a span of const char from a SEXP
inline cspan charsxp_to_cspan(SEXP x) {
  return cspan(CHAR(x), checked_size_from_r_xlen(Rf_xlength(x)));
}

inline std::vector<cspan> strsxp_to_cspan(CharacterVector x) {
  size_t n = checked_size_from_r_xlen(Rf_xlength(x));
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

inline bool is_global_algo(const AlignmentAlgo algo) noexcept {
  return algo == AlignmentAlgo::GlobalUnit ||
         algo == AlignmentAlgo::GlobalLinear ||
         algo == AlignmentAlgo::GlobalAffine;
}

inline bool is_anchored_algo(const AlignmentAlgo algo) noexcept {
  return algo == AlignmentAlgo::AnchoredUnit ||
         algo == AlignmentAlgo::AnchoredLinear ||
         algo == AlignmentAlgo::AnchoredAffine;
}

inline bool is_unit_substitution(const Rcpp::Nullable<Rcpp::IntegerMatrix>& cm) {
  if (cm.isNull()) return true; // treat NULL as unit substitution (match=0, mismatch=1)
  Rcpp::IntegerMatrix m = cm.get();
  if (m.nrow() != m.ncol()) return false;
  const size_t n = static_cast<size_t>(m.nrow());
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      const int expected = (i == j) ? 0 : 1;
      if (m(i, j) != expected) return false;
    }
  }
  return true;
}

inline AlignmentAlgo decide_alignment_algo(std::string mode,
                                           const Rcpp::Nullable<Rcpp::IntegerMatrix>& cost_matrix,
                                           int & gap_cost,
                                           int & gap_open_cost) {
  // normalize mode
  for (auto& c : mode) c = static_cast<char>(std::tolower(c));
  if (mode == "hm") mode = "hamming";
  else if (mode == "gb" || mode == "lv" || mode == "levenshtein") mode = "global";
  else if (mode == "an" || mode == "en" || mode == "extension") mode = "anchored";

  if (mode == "hamming") return AlignmentAlgo::Hamming;

  const bool gap_cost_provided = gap_cost != NA_INTEGER;
  bool gap_open_provided = gap_open_cost != NA_INTEGER;

  // If substitution matrix is NULL, always use unit substitution and unit gaps
  // ignoring gap_cost and gap_open_cost entirely, per refined logic.
  if (cost_matrix.isNull()) {
    if (gap_cost_provided || gap_open_provided) {
      Rcpp::warning("gap_cost and gap_open_cost are ignored when cost_matrix is NULL; provide a cost_matrix to enable custom gap penalties");
    }
    return (mode == "global") ? AlignmentAlgo::GlobalUnit : AlignmentAlgo::AnchoredUnit;
  }

  if (!gap_cost_provided && gap_open_provided) {
    Rcpp::warning("gap_open_cost is defined but gap_cost is NA; ignoring gap_open_cost");
    gap_open_provided = false;
  }

  if (!gap_cost_provided) {
    gap_cost = 1;
  }
  if (!gap_open_provided) {
    gap_open_cost = 0;
  }

  // New semantics: gap_open_cost includes the first extension.
  // We consider the model affine only when opening (incl. first extension) exceeds per-step gap cost.
  const bool affine = gap_open_cost > gap_cost;
  const bool unit_subs = is_unit_substitution(cost_matrix);
  const bool unit_gaps = (gap_cost == 1) && !affine; // unit gaps only when linear with gap=1

  if (mode == "global") {
    if (affine) return AlignmentAlgo::GlobalAffine;
    if (unit_gaps && unit_subs) return AlignmentAlgo::GlobalUnit;
    return AlignmentAlgo::GlobalLinear;
  } else { // anchored
    if (affine) return AlignmentAlgo::AnchoredAffine;
    if (unit_gaps && unit_subs) return AlignmentAlgo::AnchoredUnit;
    return AlignmentAlgo::AnchoredLinear;
  }
}

// convert cost matrix to map
inline CostMap convert_cost_matrix(IntegerMatrix cost_matrix, int gap_cost, int gap_open_cost) {
  CostMap cm;
  std::vector<char> map_elements;
  List dimnames = Rcpp::as<List>(cost_matrix.attr("dimnames"));
  CharacterVector rownames = Rcpp::as<CharacterVector>(dimnames[0]);
  const size_t nrow = checked_size_from_r_xlen(rownames.size());
  map_elements.resize(nrow);
  for(size_t i = 0; i < nrow; ++i) {
    const R_xlen_t ri = checked_r_xlen(i);
    if(rownames[ri] == "gap" || rownames[ri] == "gap_open") {
      // special tokens are ignored for substitution table
      map_elements[i] = '\0';
    } else {
      Rcpp::String s = rownames[ri];
      map_elements[i] = s.get_cstring()[0];
    }
  }
  // Fill substitution cost map using original indices for non-special rows/cols.
  std::vector<size_t> alpha_idx;
  alpha_idx.reserve(nrow);
  for(size_t i = 0; i < nrow; ++i) {
    const R_xlen_t ri = checked_r_xlen(i);
    if(rownames[ri] != "gap" && rownames[ri] != "gap_open") alpha_idx.push_back(i);
  }
  for(size_t ai = 0; ai < alpha_idx.size(); ++ai) {
    for(size_t aj = 0; aj < alpha_idx.size(); ++aj) {
      const size_t row_idx = alpha_idx[ai];
      const size_t col_idx = alpha_idx[aj];
      char ci = map_elements[row_idx];
      char cj = map_elements[col_idx];
      cm.set_subst(ci, cj, cost_matrix(row_idx, col_idx));
    }
  }
  cm.gap_cost = gap_cost;
  cm.gap_open_including_first_extension = gap_open_cost;
  return cm;
}

// convert search results to DataFrame
// The query column is reconstructed from the original R query vector, so C++
// search result containers only need to carry matches and distances.
template <typename SearchResult>
inline DataFrame seqtrie_results_to_dataframe(CharacterVector query,
                                              std::vector<SearchResult> & output) {
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
      SET_STRING_ELT(query_results, checked_r_xlen(q), STRING_ELT(query, checked_r_xlen(i)));
      auto s = ctx.match[j]->template sequence<SeqTrie::array_r<char>>();
      SET_STRING_ELT(target_results, checked_r_xlen(q), to_charsxp(s));
      dist_ptr[checked_r_xlen(q++)] = ctx.distance[j];
    }
  }
  return DataFrame::create(_["query"]    = query_results,
                           _["target"]   = target_results,
                           _["distance"] = distance_results,
                           _["stringsAsFactors"] = false);
}

#endif // seqtrie_TYPES_H
