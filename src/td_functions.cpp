#include <Rcpp.h>
#include <set>
#include <memory>
#include <tuple>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

#include "treedist_types_impl.h"

using namespace Rcpp;
using namespace RcppParallel;

////////////////////////////////////////////////////////////////////////////////
// create
// [[Rcpp::export(rng = false)]]
SEXP DNATree_create() {
  return Rcpp::XPtr<DNATree>(new DNATree, true);
}

// [[Rcpp::export(rng = false)]]
SEXP RadixTree_create() {
  return Rcpp::XPtr<RadixTree>(new RadixTree, true);
}

// [[Rcpp::export(rng = false)]]
SEXP PrefixTree_create() {
  return Rcpp::XPtr<PrefixTree>(new PrefixTree, true);
}

////////////////////////////////////////////////////////////////////////////////
// size
// [[Rcpp::export(rng = false)]]
double DNATree_size(Rcpp::XPtr<DNATree> xp) {
  return static_cast<double>(xp->size());
}

// [[Rcpp::export(rng = false)]]
double RadixTree_size(Rcpp::XPtr<RadixTree> xp) {
  return static_cast<double>(xp->size());
}

// [[Rcpp::export(rng = false)]]
double PrefixTree_size(Rcpp::XPtr<PrefixTree> xp) {
  return static_cast<double>(xp->size());
}

////////////////////////////////////////////////////////////////////////////////
// print
// [[Rcpp::export(rng = false)]]
std::string DNATree_print(Rcpp::XPtr<DNATree> xp) {
  return xp->print();
}

// [[Rcpp::export(rng = false)]]
std::string RadixTree_print(Rcpp::XPtr<RadixTree> xp) {
  return xp->print();
}

// [[Rcpp::export(rng = false)]]
std::string PrefixTree_print(Rcpp::XPtr<PrefixTree> xp) {
  return xp->print();
}

////////////////////////////////////////////////////////////////////////////////
// insert
// [[Rcpp::export(rng = false)]]
NumericVector DNATree_insert(Rcpp::XPtr<DNATree> xp, CharacterVector sequences) {
  using tree_type = DNATree;
  auto & tree = *xp.get();
  SEXP * sp = STRING_PTR(sequences);
  size_t nseqs = Rf_xlength(sequences);
  NumericVector result(nseqs);
  double * rp = REAL(result);
  for(size_t i=0; i<nseqs; ++i) {
    auto xi = tree.insert(cspan(CHAR(sp[i]), Rf_xlength(sp[i])));
    rp[i] = xi == tree_type::value_type::nullidx ? NA_REAL : static_cast<double>(xi);
  }
  return result;
}

// [[Rcpp::export(rng = false)]]
NumericVector RadixTree_insert(Rcpp::XPtr<RadixTree> xp, CharacterVector sequences) {
  using tree_type = RadixTree;
  auto & tree = *xp.get();
  SEXP * sp = STRING_PTR(sequences);
  size_t nseqs = Rf_xlength(sequences);
  NumericVector result(nseqs);
  double * rp = REAL(result);
  for(size_t i=0; i<nseqs; ++i) {
    auto xi = tree.insert(cspan(CHAR(sp[i]), Rf_xlength(sp[i])));
    rp[i] = xi == tree_type::value_type::nullidx ? NA_REAL : static_cast<double>(xi);
  }
  return result;
}

// [[Rcpp::export(rng = false)]]
NumericVector PrefixTree_insert(Rcpp::XPtr<PrefixTree> xp, CharacterVector sequences) {
  using tree_type = PrefixTree;
  auto & tree = *xp.get();
  SEXP * sp = STRING_PTR(sequences);
  size_t nseqs = Rf_xlength(sequences);
  NumericVector result(nseqs);
  double * rp = REAL(result);
  for(size_t i=0; i<nseqs; ++i) {
    auto xi = tree.insert(cspan(CHAR(sp[i]), Rf_xlength(sp[i])));
    rp[i] = xi == tree_type::value_type::nullidx ? NA_REAL : static_cast<double>(xi);
  }
  return result;
}

////////////////////////////////////////////////////////////////////////////////
// erase

// [[Rcpp::export(rng = false)]]
NumericVector DNATree_erase(Rcpp::XPtr<DNATree> xp, CharacterVector sequences) {
  using tree_type = DNATree;
  auto & tree = *xp.get();
  size_t nseqs = Rf_xlength(sequences);
  SEXP * sp = STRING_PTR(sequences);
  NumericVector result(nseqs);
  double * rp = REAL(result);
  for(size_t i=0; i<nseqs; ++i) {
    auto xi = tree.erase(cspan(CHAR(sp[i]), Rf_xlength(sp[i])));
    rp[i] = xi == tree_type::value_type::get_null_idx() ? NA_REAL : static_cast<double>(xi);
  }
  return result;
}

// [[Rcpp::export(rng = false)]]
NumericVector RadixTree_erase(Rcpp::XPtr<RadixTree> xp, CharacterVector sequences) {
  using tree_type = RadixTree;
  auto & tree = *xp.get();
  size_t nseqs = Rf_xlength(sequences);
  SEXP * sp = STRING_PTR(sequences);
  NumericVector result(nseqs);
  double * rp = REAL(result);
  for(size_t i=0; i<nseqs; ++i) {
    auto xi = tree.erase(cspan(CHAR(sp[i]), Rf_xlength(sp[i])));
    rp[i] = xi == tree_type::value_type::get_null_idx() ? NA_REAL : static_cast<double>(xi);
  }
  return result;
}

// [[Rcpp::export(rng = false)]]
NumericVector PrefixTree_erase(Rcpp::XPtr<PrefixTree> xp, CharacterVector sequences) {
  using tree_type = PrefixTree;
  auto & tree = *xp.get();
  size_t nseqs = Rf_xlength(sequences);
  SEXP * sp = STRING_PTR(sequences);
  NumericVector result(nseqs);
  double * rp = REAL(result);
  for(size_t i=0; i<nseqs; ++i) {
    auto xi = tree.erase(cspan(CHAR(sp[i]), Rf_xlength(sp[i])));
    rp[i] = xi == tree_type::value_type::get_null_idx() ? NA_REAL : static_cast<double>(xi);
  }
  return result;
}


////////////////////////////////////////////////////////////////////////////////
// find

// [[Rcpp::export(rng = false)]]
NumericVector DNATree_find(Rcpp::XPtr<DNATree> xp, CharacterVector sequences) {
  using tree_type = DNATree;
  auto & tree = *xp.get();
  size_t nseqs = Rf_xlength(sequences);
  SEXP * sp = STRING_PTR(sequences);
  NumericVector result(nseqs);
  double * rp = REAL(result);
  for(size_t i=0; i<nseqs; ++i) {
    auto xi = tree.find(cspan(CHAR(sp[i]), Rf_xlength(sp[i])));
    rp[i] = xi == tree_type::value_type::get_null_idx() ? NA_REAL : static_cast<double>(xi);
  }
  return result;
}

// [[Rcpp::export(rng = false)]]
NumericVector RadixTree_find(Rcpp::XPtr<RadixTree> xp, CharacterVector sequences) {
  using tree_type = RadixTree;
  auto & tree = *xp.get();
  size_t nseqs = Rf_xlength(sequences);
  SEXP * sp = STRING_PTR(sequences);
  NumericVector result(nseqs);
  double * rp = REAL(result);
  for(size_t i=0; i<nseqs; ++i) {
    auto xi = tree.find(cspan(CHAR(sp[i]), Rf_xlength(sp[i])));
    rp[i] = xi == tree_type::value_type::get_null_idx() ? NA_REAL : static_cast<double>(xi);
  }
  return result;
}

// [[Rcpp::export(rng = false)]]
NumericVector PrefixTree_find(Rcpp::XPtr<PrefixTree> xp, CharacterVector sequences) {
  using tree_type = PrefixTree;
  auto & tree = *xp.get();
  size_t nseqs = Rf_xlength(sequences);
  SEXP * sp = STRING_PTR(sequences);
  NumericVector result(nseqs);
  double * rp = REAL(result);
  for(size_t i=0; i<nseqs; ++i) {
    auto xi = tree.find(cspan(CHAR(sp[i]), Rf_xlength(sp[i])));
    rp[i] = xi == tree_type::value_type::get_null_idx() ? NA_REAL : static_cast<double>(xi);
  }
  return result;
}




////////////////////////////////////////////////////////////////////////////////
// to_dataframe

// [[Rcpp::export(rng = false)]]
SEXP DNATree_to_dataframe(Rcpp::XPtr<DNATree> xp) {
  return xp->to_dataframe();
}

// [[Rcpp::export(rng = false)]]
SEXP RadixTree_to_dataframe(Rcpp::XPtr<RadixTree> xp) {
  return xp->to_dataframe();
}

// [[Rcpp::export(rng = false)]]
SEXP PrefixTree_to_dataframe(Rcpp::XPtr<PrefixTree> xp) {
  return xp->to_dataframe();
}

////////////////////////////////////////////////////////////////////////////////
// levenshtein

template <class T> struct LevenshteinWorker : public Worker {
  const T & tree;
  const std::vector<cspan> & query;
  int const * const max_distance_ptr;
  std::vector<std::pair<std::vector<uint64_t>, std::vector<int>>> & output;
  LevenshteinWorker(const T & tree,
                    const std::vector<cspan> & query,
                    int const * const max_distance_ptr,
                    std::vector<std::pair<std::vector<uint64_t>, std::vector<int>>> & output) : 
    tree(tree), query(query), max_distance_ptr(max_distance_ptr), output(output) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t i=begin; i<end; ++i) {
      output[i] = tree.levenshtein(query[i], max_distance_ptr[i]);
      // progress_bar.increment();
    }
  }
};

// [[Rcpp::export(rng = false)]]
SEXP DNATree_levenshtein(Rcpp::XPtr<DNATree> xp, CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool display_progress) {
  using tree_type = DNATree;
  auto & tree = *xp.get();
  size_t nseqs = Rf_xlength(sequences);
  if(Rf_xlength(max_distance) != nseqs) throw std::runtime_error("sequences and max_distance must be same length (or length 1)");
  SEXP * sp = STRING_PTR(sequences);
  int * max_distance_ptr = INTEGER(max_distance);
  
  if(nseqs == 0) return R_NilValue;
  
  std::vector<cspan> query(nseqs);
  for(size_t i=0; i<nseqs; ++i) { query[i] = cspan(CHAR(sp[i]), Rf_xlength(sp[i])); }
  std::vector<std::pair<std::vector<uint64_t>, std::vector<int>>> output(nseqs);
  if(nthreads == 1) {
    Progress progress_bar(nseqs, display_progress);
    for(size_t i=0; i<nseqs; ++i) {
      output[i] = tree.levenshtein(query[i], max_distance_ptr[i]);
      progress_bar.increment();
    }
  } else {
    LevenshteinWorker<tree_type> w(tree, query, max_distance_ptr, output);
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
      cspan tj = tree.index_to_sequence(targets[j]);
      SET_STRING_ELT(target_results, q, Rf_mkCharLen(tj.data(), tj.size()));
      distance_results_ptr[q] = distances[j];
      q++;
    }
  }
  return DataFrame::create(_["query"] = query_results, _["target"] = target_results, _["distance"] = distance_results);
}

// [[Rcpp::export(rng = false)]]
SEXP RadixTree_levenshtein(Rcpp::XPtr<RadixTree> xp, CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool display_progress) {
  using tree_type = RadixTree;
  auto & tree = *xp.get();
  size_t nseqs = Rf_xlength(sequences);
  if(Rf_xlength(max_distance) != nseqs) throw std::runtime_error("sequences and max_distance must be same length (or length 1)");
  SEXP * sp = STRING_PTR(sequences);
  int * max_distance_ptr = INTEGER(max_distance);
  
  if(nseqs == 0) return R_NilValue;
  
  std::vector<cspan> query(nseqs);
  for(size_t i=0; i<nseqs; ++i) { query[i] = cspan(CHAR(sp[i]), Rf_xlength(sp[i])); }
  std::vector<std::pair<std::vector<uint64_t>, std::vector<int>>> output(nseqs);
  if(nthreads == 1) {
    Progress progress_bar(nseqs, display_progress);
    for(size_t i=0; i<nseqs; ++i) {
      output[i] = tree.levenshtein(query[i], max_distance_ptr[i]);
      progress_bar.increment();
    }
  } else {
    LevenshteinWorker<tree_type> w(tree, query, max_distance_ptr, output);
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
      cspan tj = tree.index_to_sequence(targets[j]);
      SET_STRING_ELT(target_results, q, Rf_mkCharLen(tj.data(), tj.size()));
      distance_results_ptr[q] = distances[j];
      q++;
    }
  }
  return DataFrame::create(_["query"] = query_results, _["target"] = target_results, _["distance"] = distance_results);
}

// [[Rcpp::export(rng = false)]]
SEXP PrefixTree_levenshtein(Rcpp::XPtr<PrefixTree> xp, CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool display_progress) {
  using tree_type = PrefixTree;
  auto & tree = *xp.get();
  size_t nseqs = Rf_xlength(sequences);
  if(Rf_xlength(max_distance) != nseqs) throw std::runtime_error("sequences and max_distance must be same length (or length 1)");
  SEXP * sp = STRING_PTR(sequences);
  int * max_distance_ptr = INTEGER(max_distance);
  
  if(nseqs == 0) return R_NilValue;
  
  std::vector<cspan> query(nseqs);
  for(size_t i=0; i<nseqs; ++i) { query[i] = cspan(CHAR(sp[i]), Rf_xlength(sp[i])); }
  std::vector<std::pair<std::vector<uint64_t>, std::vector<int>>> output(nseqs);
  if(nthreads == 1) {
    Progress progress_bar(nseqs, display_progress);
    for(size_t i=0; i<nseqs; ++i) {
      output[i] = tree.levenshtein(query[i], max_distance_ptr[i]);
      progress_bar.increment();
    }
  } else {
    LevenshteinWorker<tree_type> w(tree, query, max_distance_ptr, output);
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
      cspan tj = tree.index_to_sequence(targets[j]);
      SET_STRING_ELT(target_results, q, Rf_mkCharLen(tj.data(), tj.size()));
      distance_results_ptr[q] = distances[j];
      q++;
    }
  }
  return DataFrame::create(_["query"] = query_results, _["target"] = target_results, _["distance"] = distance_results);
}

////////////////////////////////////////////////////////////////////////////////
// hamming
template <class T> struct HammingWorker : public Worker {
  const T & tree;
  const std::vector<cspan> & query;
  int const * const max_distance_ptr;
  std::vector<std::pair<std::vector<uint64_t>, std::vector<int>>> & output;
  HammingWorker(const T & tree,
                    const std::vector<cspan> & query,
                    int const * const max_distance_ptr,
                    std::vector<std::pair<std::vector<uint64_t>, std::vector<int>>> & output) : 
    tree(tree), query(query), max_distance_ptr(max_distance_ptr), output(output) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t i=begin; i<end; ++i) {
      output[i] = tree.hamming(query[i], max_distance_ptr[i]);
      // progress_bar.increment();
    }
  }
};

// [[Rcpp::export(rng = false)]]
SEXP DNATree_hamming(Rcpp::XPtr<DNATree> xp, CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool display_progress) {
  using tree_type = DNATree;
  auto & tree = *xp.get();
  size_t nseqs = Rf_xlength(sequences);
  if(Rf_xlength(max_distance) != nseqs) throw std::runtime_error("sequences and max_distance must be same length (or length 1)");
  SEXP * sp = STRING_PTR(sequences);
  int * max_distance_ptr = INTEGER(max_distance);
  
  if(nseqs == 0) return R_NilValue;
  
  std::vector<cspan> query(nseqs);
  for(size_t i=0; i<nseqs; ++i) { query[i] = cspan(CHAR(sp[i]), Rf_xlength(sp[i])); }
  std::vector<std::pair<std::vector<uint64_t>, std::vector<int>>> output(nseqs);
  if(nthreads == 1) {
    Progress progress_bar(nseqs, display_progress);
    for(size_t i=0; i<nseqs; ++i) {
      output[i] = tree.hamming(query[i], max_distance_ptr[i]);
      progress_bar.increment();
    }
  } else {
    HammingWorker<tree_type> w(tree, query, max_distance_ptr, output);
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
      cspan tj = tree.index_to_sequence(targets[j]);
      SET_STRING_ELT(target_results, q, Rf_mkCharLen(tj.data(), tj.size()));
      distance_results_ptr[q] = distances[j];
      q++;
    }
  }
  return DataFrame::create(_["query"] = query_results, _["target"] = target_results, _["distance"] = distance_results);
}

// [[Rcpp::export(rng = false)]]
SEXP RadixTree_hamming(Rcpp::XPtr<RadixTree> xp, CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool display_progress) {
  using tree_type = RadixTree;
  auto & tree = *xp.get();
  size_t nseqs = Rf_xlength(sequences);
  if(Rf_xlength(max_distance) != nseqs) throw std::runtime_error("sequences and max_distance must be same length (or length 1)");
  SEXP * sp = STRING_PTR(sequences);
  int * max_distance_ptr = INTEGER(max_distance);
  
  if(nseqs == 0) return R_NilValue;
  
  std::vector<cspan> query(nseqs);
  for(size_t i=0; i<nseqs; ++i) { query[i] = cspan(CHAR(sp[i]), Rf_xlength(sp[i])); }
  std::vector<std::pair<std::vector<uint64_t>, std::vector<int>>> output(nseqs);
  if(nthreads == 1) {
    Progress progress_bar(nseqs, display_progress);
    for(size_t i=0; i<nseqs; ++i) {
      output[i] = tree.hamming(query[i], max_distance_ptr[i]);
      progress_bar.increment();
    }
  } else {
    HammingWorker<tree_type> w(tree, query, max_distance_ptr, output);
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
      cspan tj = tree.index_to_sequence(targets[j]);
      SET_STRING_ELT(target_results, q, Rf_mkCharLen(tj.data(), tj.size()));
      distance_results_ptr[q] = distances[j];
      q++;
    }
  }
  return DataFrame::create(_["query"] = query_results, _["target"] = target_results, _["distance"] = distance_results);
}

// [[Rcpp::export(rng = false)]]
SEXP PrefixTree_hamming(Rcpp::XPtr<PrefixTree> xp, CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool display_progress) {
  using tree_type = PrefixTree;
  auto & tree = *xp.get();
  size_t nseqs = Rf_xlength(sequences);
  if(Rf_xlength(max_distance) != nseqs) throw std::runtime_error("sequences and max_distance must be same length (or length 1)");
  SEXP * sp = STRING_PTR(sequences);
  int * max_distance_ptr = INTEGER(max_distance);
  
  if(nseqs == 0) return R_NilValue;
  
  std::vector<cspan> query(nseqs);
  for(size_t i=0; i<nseqs; ++i) { query[i] = cspan(CHAR(sp[i]), Rf_xlength(sp[i])); }
  std::vector<std::pair<std::vector<uint64_t>, std::vector<int>>> output(nseqs);
  if(nthreads == 1) {
    Progress progress_bar(nseqs, display_progress);
    for(size_t i=0; i<nseqs; ++i) {
      output[i] = tree.hamming(query[i], max_distance_ptr[i]);
      progress_bar.increment();
    }
  } else {
    HammingWorker<tree_type> w(tree, query, max_distance_ptr, output);
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
      cspan tj = tree.index_to_sequence(targets[j]);
      SET_STRING_ELT(target_results, q, Rf_mkCharLen(tj.data(), tj.size()));
      distance_results_ptr[q] = distances[j];
      q++;
    }
  }
  return DataFrame::create(_["query"] = query_results, _["target"] = target_results, _["distance"] = distance_results);
}




