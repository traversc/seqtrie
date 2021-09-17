#ifndef TREEDIST_TYPES_IMPL_H
#define TREEDIST_TYPES_IMPL_H


#include <Rcpp.h>
#include <RcppParallel.h>

#include <cstring>
#include <utility>

#include <treedist/prefixmap.h>
#include <treedist/radixarray.h>
#include <treedist/radixmap.h>

#include "treedist_types.h"
#include "simple_progress/simple_progress.h"

using namespace Rcpp;
using namespace RcppParallel;

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


template <class T> struct LevenshteinWorker : public Worker {
  const std::unique_ptr<T> & root;
  const std::vector<cspan> & query;
  int const * const max_distance_ptr;
  std::vector<typename T::Levenshtein::result_type> & output;
  trqwe::simple_progress & progress_bar;
  LevenshteinWorker(const std::unique_ptr<T> & root,
                    const std::vector<cspan> & query,
                    int const * const max_distance_ptr,
                    std::vector<std::pair<std::vector<typename T::index_type>, std::vector<int>>> & output,
                    trqwe::simple_progress & progress_bar) : 
    root(root), query(query), max_distance_ptr(max_distance_ptr), output(output), progress_bar(progress_bar) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t i=begin; i<end; ++i) {
      auto usequence = dispatch<T>::sequence_convert(query[i]);
      output[i] = typename T::Levenshtein(root, usequence, max_distance_ptr[i]).search();
      progress_bar.increment();
    }
  }
};

template <class T> struct HammingWorker : public Worker {
  const std::unique_ptr<T> & root;
  const std::vector<cspan> & query;
  int const * const max_distance_ptr;
  std::vector<typename T::Levenshtein::result_type> & output;
  trqwe::simple_progress & progress_bar;
  HammingWorker(const std::unique_ptr<T> & root,
                const std::vector<cspan> & query,
                int const * const max_distance_ptr,
                std::vector<typename T::Levenshtein::result_type> & output,
                trqwe::simple_progress & progress_bar) : 
    root(root), query(query), max_distance_ptr(max_distance_ptr), output(output), progress_bar(progress_bar) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t i=begin; i<end; ++i) {
      auto usequence = dispatch<T>::sequence_convert(query[i]);
      output[i] = typename T::Hamming(root, usequence, max_distance_ptr[i]).search();
      progress_bar.increment();
    }
  }
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

template <typename T> NumericVector rtree<T>::insert(CharacterVector sequences) {
  SEXP * sequence_ptr = STRING_PTR(sequences);
  size_t nseqs = Rf_xlength(sequences);
  NumericVector result(nseqs);
  double * result_ptr = REAL(result);
  for(size_t i=0; i<nseqs; ++i) {
    cspan sequence(CHAR(sequence_ptr[i]), Rf_xlength(sequence_ptr[i]));
    auto usequence = dispatch<typename T::value_type>::sequence_convert(sequence);
    index_type idx = value_type::insert(root, usequence, sequence_map.size());
    if(idx == value_type::nullidx) {
      sequence_map.push_back(ncstring{sequence.data(), sequence.size()});
      result_ptr[i] = NA_REAL;
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
    cspan sequence(CHAR(sequence_ptr[i]), Rf_xlength(sequence_ptr[i]));
    auto usequence = dispatch<typename T::value_type>::sequence_convert(sequence);
    index_type idx = value_type::erase(root, usequence);
    if(idx == value_type::nullidx) {
      result_ptr[i] = NA_REAL;
    } else {
      sequence_map[idx].nullify();
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
    cspan sequence(CHAR(sequence_ptr[i]), Rf_xlength(sequence_ptr[i]));
    auto usequence = dispatch<typename T::value_type>::sequence_convert(sequence);
    index_type idx = value_type::find(root, usequence);
    result_ptr[i] = idx == value_type::nullidx ? NA_REAL : static_cast<double>(idx);
  }
  return result;
}

template <typename T> SEXP rtree<T>::find_prefix(CharacterVector sequences) const {
  SEXP * sequence_ptr = STRING_PTR(sequences);
  size_t nseqs = Rf_xlength(sequences);
  std::vector< std::vector<index_type> > output(nseqs);
  
  if(nseqs == 0) return R_NilValue;
  
  for(size_t i=0; i<nseqs; ++i) {
    cspan sequence(CHAR(sequence_ptr[i]), Rf_xlength(sequence_ptr[i]));
    auto usequence = dispatch<typename T::value_type>::sequence_convert(sequence);
    output[i] = value_type::find_prefix(root, usequence);
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
      cspan tj = index_to_sequence(targets[j]);
      SET_STRING_ELT(target_results, q, Rf_mkCharLen(tj.data(), tj.size()));
      q++;
    }
  }
  return DataFrame::create(_["query"] = query_results, _["target"] = target_results, _["stringsAsFactors"] = false);
}

template <typename T> SEXP rtree<T>::levenshtein_search(CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress) const {
  size_t nseqs = Rf_xlength(sequences);
  if(Rf_xlength(max_distance) != nseqs) throw std::runtime_error("sequences and max_distance must be same length (or length 1)");
  SEXP * sp = STRING_PTR(sequences);
  int * max_distance_ptr = INTEGER(max_distance);
  
  if(nseqs == 0) return R_NilValue;
  
  std::vector<cspan> query(nseqs);
  for(size_t i=0; i<nseqs; ++i) { query[i] = cspan(CHAR(sp[i]), Rf_xlength(sp[i])); }
  std::vector<typename levenshtein_type::result_type> output(nseqs);
  trqwe::simple_progress progress_bar(nseqs, show_progress);
  
  if(nthreads == 1) {
    for(size_t i=0; i<nseqs; ++i) {
      auto usequence = dispatch<value_type>::sequence_convert(query[i]);
      output[i] = levenshtein_type(root, usequence, max_distance_ptr[i]).search();
      progress_bar.increment();
    }
  } else {
    LevenshteinWorker<value_type> w(root, query, max_distance_ptr, output, progress_bar);
    parallelFor(0, nseqs, w, 1, nthreads);
  }
  
  size_t nresults = 0;
  for(size_t i=0; i<nseqs; ++i) { nresults += output[i].first.size(); }
  CharacterVector query_results(nresults);
  CharacterVector target_results(nresults);
  IntegerVector distance_results(nresults);
  int * distance_results_ptr = INTEGER(distance_results);
  size_t q = 0;
  for(size_t i=0; i<nseqs; ++i) {
    auto & targets = output[i].first;
    auto & distances = output[i].second;
    for(size_t j=0; j<targets.size(); ++j) {
      SET_STRING_ELT(query_results, q, STRING_ELT(sequences, i));
      cspan tj = index_to_sequence(targets[j]);
      SET_STRING_ELT(target_results, q, Rf_mkCharLen(tj.data(), tj.size()));
      distance_results_ptr[q] = distances[j];
      q++;
    }
  }
  return DataFrame::create(_["query"] = query_results, _["target"] = target_results, _["distance"] = distance_results, _["stringsAsFactors"] = false);
}

template <typename T> SEXP rtree<T>::hamming_search(CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress) const {
  size_t nseqs = Rf_xlength(sequences);
  if(Rf_xlength(max_distance) != nseqs) throw std::runtime_error("sequences and max_distance must be same length (or length 1)");
  SEXP * sp = STRING_PTR(sequences);
  int * max_distance_ptr = INTEGER(max_distance);
  
  if(nseqs == 0) return R_NilValue;
  
  std::vector<cspan> query(nseqs);
  for(size_t i=0; i<nseqs; ++i) { query[i] = cspan(CHAR(sp[i]), Rf_xlength(sp[i])); }
  std::vector<typename hamming_type::result_type> output(nseqs);
  trqwe::simple_progress progress_bar(nseqs, show_progress);
  
  if(nthreads == 1) {
    for(size_t i=0; i<nseqs; ++i) {
      auto usequence = dispatch<value_type>::sequence_convert(query[i]);
      output[i] = hamming_type(root, usequence, max_distance_ptr[i]).search();
      progress_bar.increment();
    }
  } else {
    HammingWorker<value_type> w(root, query, max_distance_ptr, output, progress_bar);
    parallelFor(0, nseqs, w, 1, nthreads);
  }
  
  size_t nresults = 0;
  for(size_t i=0; i<nseqs; ++i) { nresults += output[i].first.size(); }
  CharacterVector query_results(nresults);
  CharacterVector target_results(nresults);
  IntegerVector distance_results(nresults);
  int * distance_results_ptr = INTEGER(distance_results);
  size_t q = 0;
  for(size_t i=0; i<nseqs; ++i) {
    auto & targets = output[i].first;
    auto & distances = output[i].second;
    for(size_t j=0; j<targets.size(); ++j) {
      SET_STRING_ELT(query_results, q, STRING_ELT(sequences, i));
      cspan tj = index_to_sequence(targets[j]);
      SET_STRING_ELT(target_results, q, Rf_mkCharLen(tj.data(), tj.size()));
      distance_results_ptr[q] = distances[j];
      q++;
    }
  }
  return DataFrame::create(_["query"] = query_results, _["target"] = target_results, _["distance"] = distance_results, _["stringsAsFactors"] = false);
}

template <typename T> cspan rtree<T>::index_to_sequence(const index_type idx) const {
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
