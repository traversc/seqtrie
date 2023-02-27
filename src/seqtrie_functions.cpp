#include <Rcpp.h>
#include <RcppParallel.h>

#include <set>
#include <memory>
#include <tuple>
#include "seqtrie_types.h"

#include "simple_progress/simple_progress.h"

#include <boost/numeric/ublas/matrix.hpp>
using IMatrix = boost::numeric::ublas::matrix<int>;

using namespace Rcpp;
using namespace RcppParallel;


template <typename T> SEXP to_charsxp(const T & x) { return Rf_mkCharLen(x.data(), x.size()); }
cspan charsxp_to_cspan(SEXP x) { return cspan(CHAR(x), Rf_xlength(x)); }


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

template <class T> struct AnchoredWorker : public Worker {
  const T & root;
  const std::vector<cspan> & query;
  int const * const max_distance_ptr;
  std::vector<typename T::search_context> & output;
  trqwe::simple_progress & progress_bar;
  AnchoredWorker(const T & root,
                const std::vector<cspan> & query,
                int const * const max_distance_ptr,
                std::vector<typename T::search_context> & output,
                trqwe::simple_progress & progress_bar) : 
    root(root), query(query), max_distance_ptr(max_distance_ptr), output(output), progress_bar(progress_bar) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t i=begin; i<end; ++i) {
      output[i] = root.anchored_search(query[i], max_distance_ptr[i]);
      progress_bar.increment();
    }
  }
};

////////////////////////////////////////////////////////////////////////////////
// RadixTreeContext class definitions 

template <typename T> double RadixTreeContext<T>::size() const {
  return static_cast<double>(root.size());
}

template <typename T> LogicalVector RadixTreeContext<T>::insert(CharacterVector sequences) {
  SEXP * sequence_ptr = STRING_PTR(sequences);
  size_t nseqs = Rf_xlength(sequences);
  LogicalVector result(nseqs);
  int * result_ptr = LOGICAL(result);
  for(size_t i=0; i<nseqs; ++i) {
    cspan sequence = charsxp_to_cspan(sequence_ptr[i]);
    size_t idx = root.insert(sequence, posidx);
    result_ptr[i] = idx == nullidx ? 1 : 0; // nullidx means it was successfully inserted
  }
  return result;
}

template <typename T> LogicalVector RadixTreeContext<T>::erase(CharacterVector sequences) {
  SEXP * sequence_ptr = STRING_PTR(sequences);
  size_t nseqs = Rf_xlength(sequences);
  LogicalVector result(nseqs);
  int * result_ptr = LOGICAL(result);
  for(size_t i=0; i<nseqs; ++i) {
    cspan sequence = charsxp_to_cspan(sequence_ptr[i]);
    size_t idx = root.erase(sequence);
    result_ptr[i] = idx == nullidx ? 0 : 1; // nullidx means sequence did not exist, erase was not succesful
  }
  return result;
}

template <typename T> LogicalVector RadixTreeContext<T>::find(CharacterVector sequences) const {
  SEXP * sequence_ptr = STRING_PTR(sequences);
  size_t nseqs = Rf_xlength(sequences);
  LogicalVector result(nseqs);
  int * result_ptr = LOGICAL(result);
  for(size_t i=0; i<nseqs; ++i) {
    cspan sequence = charsxp_to_cspan(sequence_ptr[i]);
    size_t idx = root.find(sequence);
    result_ptr[i] = idx == nullidx ? 0 : 1; // nullidx means sequence was not found
  }
  return result;
}

template <typename T> SEXP RadixTreeContext<T>::levenshtein_search(CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress) const {
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

template <typename T> SEXP RadixTreeContext<T>::hamming_search(CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress) const {
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

template <typename T> SEXP RadixTreeContext<T>::anchored_search(CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress) const {
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
      output[i] = root.anchored_search(query[i], max_distance_ptr[i]);
      progress_bar.increment();
    }
  } else {
    AnchoredWorker<node_type> w(root, query, max_distance_ptr, output, progress_bar);
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

template <typename T> SEXP RadixTreeContext<T>::prefix_search(CharacterVector sequences) const {
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

template <typename T> std::string RadixTreeContext<T>::print() const {
  return root.print();
}

template <typename T> SEXP RadixTreeContext<T>::graph(const double max_depth) const {
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

template <typename T> SEXP RadixTreeContext<T>::to_vector() const {
  auto seqs = root.all();
  if(seqs.size() == 0) return R_NilValue;
  CharacterVector sequence(seqs.size());
  for(size_t i=0; i<seqs.size(); ++i) {
    auto s = seqs[i]->template sequence<trqwe::small_array<char>>();
    SET_STRING_ELT(sequence, i, to_charsxp(s));
  }
  return sequence;
}

template <typename T> bool RadixTreeContext<T>::validate() const {
  return root.validate();
}


// [[Rcpp::export(rng = false)]]
SEXP RadixTree_create() { return Rcpp::XPtr<RadixTreeCtxForR>(new RadixTreeCtxForR, true); }

// [[Rcpp::export(rng = false)]]
double RadixTree_size(Rcpp::XPtr<RadixTreeCtxForR> xp) { return xp->size(); }

// [[Rcpp::export(rng = false)]]
std::string RadixTree_print(Rcpp::XPtr<RadixTreeCtxForR> xp) { return xp->print(); }

// [[Rcpp::export(rng = false)]]
SEXP RadixTree_graph(Rcpp::XPtr<RadixTreeCtxForR> xp, const double max_depth) { return xp->graph(max_depth); }

// [[Rcpp::export(rng = false)]]
LogicalVector RadixTree_insert(Rcpp::XPtr<RadixTreeCtxForR> xp, CharacterVector sequences) { return xp->insert(sequences); }

// [[Rcpp::export(rng = false)]]
LogicalVector RadixTree_erase(Rcpp::XPtr<RadixTreeCtxForR> xp, CharacterVector sequences) { return xp->erase(sequences); }

// [[Rcpp::export(rng = false)]]
LogicalVector RadixTree_find(Rcpp::XPtr<RadixTreeCtxForR> xp, CharacterVector sequences) { return xp->find(sequences); }

// [[Rcpp::export(rng = false)]]
SEXP RadixTree_to_vector(Rcpp::XPtr<RadixTreeCtxForR> xp) { return xp->to_vector(); }

// [[Rcpp::export(rng = false)]]
SEXP RadixTree_levenshtein_search(Rcpp::XPtr<RadixTreeCtxForR> xp, CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress) {
  return xp->levenshtein_search(sequences, max_distance, nthreads, show_progress);
}

// [[Rcpp::export(rng = false)]]
SEXP RadixTree_hamming_search(Rcpp::XPtr<RadixTreeCtxForR> xp, CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress) {
  return xp->hamming_search(sequences, max_distance, nthreads, show_progress);
}

// [[Rcpp::export(rng = false)]]
SEXP RadixTree_anchored_search(Rcpp::XPtr<RadixTreeCtxForR> xp, CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress) {
  return xp->anchored_search(sequences, max_distance, nthreads, show_progress);
}

// [[Rcpp::export(rng = false)]]
SEXP RadixTree_prefix_search(Rcpp::XPtr<RadixTreeCtxForR> xp, CharacterVector sequences) { return xp->prefix_search(sequences); }

// [[Rcpp::export(rng = false)]]
bool RadixTree_validate(Rcpp::XPtr<RadixTreeCtxForR> xp) { return xp->validate(); }

////////////////////////////////////////////////////////////////////////////////
// distance matrix functions

IMatrix get_dprog_matrix(cspan query, cspan target) {
  IMatrix mat(query.size()+1, target.size()+1);
  // boundaries are the same for both levenshtein and anchored
  for(size_t i=0; i<mat.size1(); ++i) mat(i,0) = i;
  for(size_t j=1; j<mat.size2(); ++j) mat(0,j) = j;

  // fill it in
  for(size_t i=1; i<mat.size1(); ++i) {
    for(size_t j=1; j<mat.size2(); ++j) {
      int match_cost  = mat(i-1, j-1) + (query[i-1] == target[j-1] ? 0 : 1);
      int insert_cost = mat(i, j-1) + 1;
      int delete_cost = mat(i-1, j) + 1;
      mat(i,j) = std::min({match_cost, insert_cost, delete_cost});
    }
  }
  return mat;
}

int hamming_distance(cspan query, cspan target) {
  if(query.size() != target.size()) return NA_INTEGER;
  int distance = 0;
  for(size_t i=0; i<query.size(); ++i) {
    if(query[i] != target[i]) distance++;
  }
  return distance;
}

int levenshtein_distance(cspan query, cspan target) {
  IMatrix mat = get_dprog_matrix(query, target);
  return mat(mat.size1()-1, mat.size2()-1);
}

int anchored_distance(cspan query, cspan target) {
  IMatrix mat = get_dprog_matrix(query, target);
  int distance = std::numeric_limits<int>::max();
  for(size_t i=0; i<mat.size1(); ++i) distance = std::min(distance, mat(i, mat.size2()-1));
  for(size_t j=0; j<mat.size2(); ++j) distance = std::min(distance, mat(mat.size1()-1, j));
  return distance;
}

// if we implement more distance metrics, might want to use if constexpr (require C++17)
// the output matrix is parallelized over target because R matrices are column major

struct HammingDistanceWorker : public Worker {
  const std::vector<cspan> & query;
  const std::vector<cspan> & target;
  int * output; // IntegerMatrix
  trqwe::simple_progress & progress_bar;
  HammingDistanceWorker(const std::vector<cspan> & query,
                        const std::vector<cspan> & target,
                        int * output,
                        trqwe::simple_progress & progress_bar) : 
  query(query), target(target), output(output), progress_bar(progress_bar) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t j=begin; j<end; ++j) {
      for(size_t i=0; i<query.size(); ++i) {
        output[i + j*query.size()] = hamming_distance(query[i], target[j]);
      }
    }
  }
};

struct LevenshteinDistanceWorker : public Worker {
  const std::vector<cspan> & query;
  const std::vector<cspan> & target;
  int * output; // IntegerMatrix
  trqwe::simple_progress & progress_bar;
  LevenshteinDistanceWorker(const std::vector<cspan> & query,
                        const std::vector<cspan> & target,
                        int * output,
                        trqwe::simple_progress & progress_bar) : 
  query(query), target(target), output(output), progress_bar(progress_bar) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t j=begin; j<end; ++j) {
      for(size_t i=0; i<query.size(); ++i) {
        output[i + j*query.size()] = levenshtein_distance(query[i], target[j]);
      }
    }
  }
};

struct AnchoredDistanceWorker : public Worker {
  const std::vector<cspan> & query;
  const std::vector<cspan> & target;
  int * output; // IntegerMatrix
  trqwe::simple_progress & progress_bar;
  AnchoredDistanceWorker(const std::vector<cspan> & query,
                        const std::vector<cspan> & target,
                        int * output,
                        trqwe::simple_progress & progress_bar) : 
  query(query), target(target), output(output), progress_bar(progress_bar) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t j=begin; j<end; ++j) {
      for(size_t i=0; i<query.size(); ++i) {
        output[i + j*query.size()] = anchored_distance(query[i], target[j]);
      }
    }
  }
};

// [[Rcpp::export(rng = false)]]
IntegerMatrix distance_matrix(CharacterVector query, CharacterVector target, 
                              const std::string mode = "levenshtein", 
                              const int nthreads = 1, const bool show_progress = false) {
  if((mode != "levenshtein") && (mode != "hamming") && (mode != "anchored")) {
    throw std::runtime_error("Metric must be one of levenshtein, hamming or anchored");
  }
  size_t query_len = Rf_xlength(query);
  SEXP * query_ptr = STRING_PTR(query);
  size_t target_len = Rf_xlength(target);
  SEXP * target_ptr = STRING_PTR(target);
  std::vector<cspan> query_span(query_len);
  for(size_t i=0; i<query_len; ++i) query_span[i] = charsxp_to_cspan(query_ptr[i]);
  std::vector<cspan> target_span(target_len);
  for(size_t i=0; i<target_len; ++i) target_span[i] = charsxp_to_cspan(target_ptr[i]);

  IntegerMatrix output(query_len, target_len);
  int * output_ptr = INTEGER(output);
  trqwe::simple_progress progress_bar(target_len, show_progress);
  
  if(mode == "levenshtein") {
    LevenshteinDistanceWorker w(query_span, target_span, output_ptr, progress_bar);
    parallelFor(0, target_len, w, 1, nthreads);
  } else if(mode == "hamming") {
    HammingDistanceWorker w(query_span, target_span, output_ptr, progress_bar);
    parallelFor(0, target_len, w, 1, nthreads);
  } else if(mode == "anchored") {
    AnchoredDistanceWorker w(query_span, target_span, output_ptr, progress_bar);
    parallelFor(0, target_len, w, 1, nthreads);
  }
  return output;
}
