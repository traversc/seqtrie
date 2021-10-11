#ifndef seqtrie_TYPES_IMPL_H
#define seqtrie_TYPES_IMPL_H


#include <Rcpp.h>
#include <RcppParallel.h>

#include <cstring>
#include <utility>

// #include <seqtrie/prefixmap.h>
// #include <seqtrie/radixarray.h>
#include <seqtrie/radixmap.h>

#include "seqtrie_types.h"
#include "simple_progress/simple_progress.h"

using namespace Rcpp;
using namespace RcppParallel;

template <class T> struct LevenshteinWorker : public Worker {
  const T & root;
  const std::vector<cspan> & query;
  int const * const max_distance_ptr;
  std::vector<typename T::search_context> & output;
  trqwe::simple_progress & progress_bar;
  LevenshteinWorker(const T & root,
                const std::vector<cspan> & query,
                int const * const max_distance_ptr,
                std::vector<typename T::search_context> & output,
                trqwe::simple_progress & progress_bar) : 
    root(root), query(query), max_distance_ptr(max_distance_ptr), output(output), progress_bar(progress_bar) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t i=begin; i<end; ++i) {
      output[i] = root.levenshtein_search(query[i], max_distance_ptr[i]);
      progress_bar.increment();
    }
  }
};

template <class T> struct HammingWorker : public Worker {
  const T & root;
  const std::vector<cspan> & query;
  int const * const max_distance_ptr;
  std::vector<typename T::search_context> & output;
  trqwe::simple_progress & progress_bar;
  HammingWorker(const T & root,
                const std::vector<cspan> & query,
                int const * const max_distance_ptr,
                std::vector<typename T::search_context> & output,
                trqwe::simple_progress & progress_bar) : 
    root(root), query(query), max_distance_ptr(max_distance_ptr), output(output), progress_bar(progress_bar) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t i=begin; i<end; ++i) {
      output[i] = root.hamming_search(query[i], max_distance_ptr[i]);
      progress_bar.increment();
    }
  }
};

template <typename T> SEXP to_charsxp(const T & x) { return Rf_mkCharLen(x.data(), x.size()); }
cspan charsxp_to_cspan(SEXP x) { return cspan(CHAR(x), Rf_xlength(x)); }

////////////////////////////////////////////////////////////////////////////////
// rtree class definitions 

template <typename T> double rtree<T>::size() const {
  return static_cast<double>(root.size());
}

template <typename T> NumericVector rtree<T>::insert(CharacterVector sequences) {
  SEXP * sequence_ptr = STRING_PTR(sequences);
  size_t nseqs = Rf_xlength(sequences);
  NumericVector result(nseqs);
  double * result_ptr = REAL(result);
  for(size_t i=0; i<nseqs; ++i) {
    cspan sequence = charsxp_to_cspan(sequence_ptr[i]);
    size_t idx = root.insert(sequence, current_idx);
    if(idx == nullidx) {
      result_ptr[i] = NA_REAL;
      current_idx++;
    } else {
      result_ptr[i] = static_cast<double>(idx);
    }
  }
  return result;
}

template <typename T> NumericVector rtree<T>::erase(CharacterVector sequences) {
  SEXP * sequence_ptr = STRING_PTR(sequences);
  size_t nseqs = Rf_xlength(sequences);
  NumericVector result(nseqs);
  double * result_ptr = REAL(result);
  for(size_t i=0; i<nseqs; ++i) {
    cspan sequence = charsxp_to_cspan(sequence_ptr[i]);
    size_t idx = root.erase(sequence);
    if(idx == nullidx) {
      result_ptr[i] = NA_REAL;
    } else {
      result_ptr[i] = static_cast<double>(idx);
    }
  }
  return result;
}

template <typename T> NumericVector rtree<T>::find(CharacterVector sequences) const {
  SEXP * sequence_ptr = STRING_PTR(sequences);
  size_t nseqs = Rf_xlength(sequences);
  NumericVector result(nseqs);
  double * result_ptr = REAL(result);
  for(size_t i=0; i<nseqs; ++i) {
    cspan sequence = charsxp_to_cspan(sequence_ptr[i]);
    size_t idx = root.find(sequence);
    result_ptr[i] = idx == nullidx ? NA_REAL : static_cast<double>(idx);
  }
  return result;
}

template <typename T> SEXP rtree<T>::levenshtein_search(CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress) const {
  size_t nseqs = Rf_xlength(sequences);
  SEXP * sequence_ptr = STRING_PTR(sequences);
  int * max_distance_ptr = INTEGER(max_distance);
  
  if(nseqs == 0) return R_NilValue;
  
  std::vector<cspan> query(nseqs);
  for(size_t i=0; i<nseqs; ++i) { query[i] = charsxp_to_cspan(sequence_ptr[i]); }
  std::vector<search_context> output(nseqs);
  trqwe::simple_progress progress_bar(nseqs, show_progress);
  
  if(nthreads == 1) {
    for(size_t i=0; i<nseqs; ++i) {
      output[i] = root.levenshtein_search(query[i], max_distance_ptr[i]);
      progress_bar.increment();
    }
  } else {
    LevenshteinWorker<node_type> w(root, query, max_distance_ptr, output, progress_bar);
    parallelFor(0, nseqs, w, 1, nthreads);
  }
  
  size_t nresults = 0;
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
      SET_STRING_ELT(query_results, q, STRING_ELT(sequences, i));
      auto s = targets[j]->template sequence<trqwe::small_array<char>>();
      SET_STRING_ELT(target_results, q, to_charsxp(s));
      distance_results_ptr[q] = distances[j];
      q++;
    }
  }
  return DataFrame::create(_["query"] = query_results, _["target"] = target_results, _["distance"] = distance_results, _["stringsAsFactors"] = false);
}

template <typename T> SEXP rtree<T>::hamming_search(CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress) const {
  size_t nseqs = Rf_xlength(sequences);
  SEXP * sequence_ptr = STRING_PTR(sequences);
  int * max_distance_ptr = INTEGER(max_distance);
  
  if(nseqs == 0) return R_NilValue;
  
  std::vector<cspan> query(nseqs);
  for(size_t i=0; i<nseqs; ++i) { query[i] = charsxp_to_cspan(sequence_ptr[i]); }
  std::vector<search_context> output(nseqs);
  trqwe::simple_progress progress_bar(nseqs, show_progress);
  
  if(nthreads == 1) {
    for(size_t i=0; i<nseqs; ++i) {
      output[i] = root.hamming_search(query[i], max_distance_ptr[i]);
      progress_bar.increment();
    }
  } else {
    HammingWorker<node_type> w(root, query, max_distance_ptr, output, progress_bar);
    parallelFor(0, nseqs, w, 1, nthreads);
  }
  
  size_t nresults = 0;
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
      SET_STRING_ELT(query_results, q, STRING_ELT(sequences, i));
      auto s = targets[j]->template sequence<trqwe::small_array<char>>();
      SET_STRING_ELT(target_results, q, to_charsxp(s));
      distance_results_ptr[q] = distances[j];
      q++;
    }
  }
  return DataFrame::create(_["query"] = query_results, _["target"] = target_results, _["distance"] = distance_results, _["stringsAsFactors"] = false);
}

template <typename T> SEXP rtree<T>::prefix_search(CharacterVector sequences) const {
  SEXP * sequence_ptr = STRING_PTR(sequences);
  size_t nseqs = Rf_xlength(sequences);
  std::vector<std::vector<path>> output(nseqs);
  
  if(nseqs == 0) return R_NilValue;
  
  for(size_t i=0; i<nseqs; ++i) {
    cspan sequence = charsxp_to_cspan(sequence_ptr[i]);
    output[i] = root.prefix_search(sequence);
  }
  
  size_t nresults = 0;
  for(size_t i=0; i<nseqs; ++i) { nresults += output[i].size(); }
  CharacterVector query_results(nresults);
  CharacterVector target_results(nresults);
  size_t q = 0;
  for(size_t i=0; i<nseqs; ++i) {
    auto & targets = output[i];
    for(size_t j=0; j<targets.size(); ++j) {
      SET_STRING_ELT(query_results, q, STRING_ELT(sequences, i));
      auto s = targets[j]->template sequence<trqwe::small_array<char>>();
      SET_STRING_ELT(target_results, q, to_charsxp(s));
      q++;
    }
  }
  return DataFrame::create(_["query"] = query_results, _["target"] = target_results, _["stringsAsFactors"] = false);
}

template <typename T> std::string rtree<T>::print() const {
  return root.print();
}

template <typename T> SEXP rtree<T>::graph(const double max_depth) const {
  size_t depth2;
  if(max_depth < 0) {
    depth2 = -1;
  } else if(max_depth >= std::numeric_limits<size_t>::max()) {
    depth2 = -1;
  } else {
    depth2 = static_cast<size_t>(max_depth);
  }
  auto seqs = root.graph(depth2);
  if(seqs.first.size() == 0) return R_NilValue;
  CharacterVector parent(seqs.first.size());
  CharacterVector child(seqs.first.size());
  for(size_t i=0; i<seqs.first.size(); ++i) {
    SET_STRING_ELT(parent, i, to_charsxp(seqs.first[i]->get_branch()));
    SET_STRING_ELT(child, i, to_charsxp(seqs.second[i]->get_branch()));
  }
  return DataFrame::create(_["parent"] = parent, _["child"] = child, _["stringsAsFactors"] = false);
}

template <typename T> SEXP rtree<T>::to_dataframe() const {
  auto seqs = root.all();
  if(seqs.size() == 0) return R_NilValue;
  CharacterVector sequence(seqs.size());
  NumericVector index(seqs.size());
  double * indexp = REAL(index);
  for(size_t i=0; i<seqs.size(); ++i) {
    auto s = seqs[i]->template sequence<trqwe::small_array<char>>();
    SET_STRING_ELT(sequence, i, to_charsxp(s));
    size_t idx = seqs[i]->get_terminal_idx();
    indexp[i] = idx == nullidx ? NA_REAL : static_cast<double>(idx);
  }
  return DataFrame::create(_["sequence"] = sequence, _["index"] = index, _["stringsAsFactors"] = false);
}

template <typename T> bool rtree<T>::validate() const {
  return root.validate();
}

// template <typename T> cspan rtree<T>::index_to_sequence(const size_t idx) const {
//   return cspan(sequence_map[idx].data(), sequence_map[idx].size()); // must not be nullptr
// }

#endif // include guard
