#ifndef seqtrie_UTILS_H
#define seqtrie_UTILS_H

#include <Rcpp.h>
#include <RcppParallel.h>
#include "seqtrie_types.h"
#include "seqtrie/radixmap.h"
#include "simple_array/small_array.h"

using namespace Rcpp;
using namespace RcppParallel;

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
inline SEXP to_charsxp(const trqwe::small_array<char> & x) {
    return Rf_mkCharLen(x.data(), x.size());
}

// Define a span of const char from a SEXP
inline cspan charsxp_to_cspan(SEXP x) {
    return cspan(CHAR(x), Rf_xlength(x));
}

inline std::vector<cspan> strsxp_to_cspan(CharacterVector x) {
    size_t n = Rf_xlength(x);
    SEXP * xp = STRING_PTR(x);
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
      auto s = targets[j]->template sequence<trqwe::small_array<char>>();
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