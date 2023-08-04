#include <Rcpp.h>
#include <RcppParallel.h>

#include <set>
#include <memory>
#include <tuple>
#include "seqtrie_types.h"

#include "simple_progress/simple_progress.h"

using namespace Rcpp;
using namespace RcppParallel;

struct HammingWorker : public Worker {
  const SeqTrie::RadixTreeR & root;
  const std::vector<cspan> & query;
  int const * const max_distance_ptr;
  std::vector<SeqTrie::search_context> & output;
  trqwe::simple_progress & progress_bar;
  HammingWorker(const SeqTrie::RadixTreeR & root,
                const std::vector<cspan> & query,
                int const * const max_distance_ptr,
                std::vector<SeqTrie::search_context> & output,
                trqwe::simple_progress & progress_bar) : 
    root(root), query(query), max_distance_ptr(max_distance_ptr), output(output), progress_bar(progress_bar) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t i=begin; i<end; ++i) {
      output[i] = root.hamming_search(query[i], max_distance_ptr[i]);
      progress_bar.increment();
    }
  }
};

struct LevenshteinWorker : public Worker {
  const SeqTrie::RadixTreeR & root;
  const std::vector<cspan> & query;
  int const * const max_distance_ptr;
  std::vector<SeqTrie::search_context> & output;
  trqwe::simple_progress & progress_bar;
  LevenshteinWorker(const SeqTrie::RadixTreeR & root,
                const std::vector<cspan> & query,
                int const * const max_distance_ptr,
                std::vector<SeqTrie::search_context> & output,
                trqwe::simple_progress & progress_bar) : 
    root(root), query(query), max_distance_ptr(max_distance_ptr), output(output), progress_bar(progress_bar) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t i=begin; i<end; ++i) {
      output[i] = root.levenshtein_search(query[i], max_distance_ptr[i]);
      progress_bar.increment();
    }
  }
};

struct AnchoredWorker : public Worker {
  const SeqTrie::RadixTreeR & root;
  const std::vector<cspan> & query;
  int const * const max_distance_ptr;
  std::vector<SeqTrie::search_context> & output;
  trqwe::simple_progress & progress_bar;
  AnchoredWorker(const SeqTrie::RadixTreeR & root,
                const std::vector<cspan> & query,
                int const * const max_distance_ptr,
                std::vector<SeqTrie::search_context> & output,
                trqwe::simple_progress & progress_bar) : 
    root(root), query(query), max_distance_ptr(max_distance_ptr), output(output), progress_bar(progress_bar) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t i=begin; i<end; ++i) {
      output[i] = root.anchored_search(query[i], max_distance_ptr[i]);
      progress_bar.increment();
    }
  }
};

struct LevenshteinWorkerWithCostMap : public Worker {
  const SeqTrie::RadixTreeR & root;
  const std::vector<cspan> & query;
  int const * const max_distance_ptr;
  std::vector<SeqTrie::search_context> & output;
  pairchar_map & cost_map;
  trqwe::simple_progress & progress_bar;
  LevenshteinWorkerWithCostMap(const SeqTrie::RadixTreeR & root,
                const std::vector<cspan> & query,
                int const * const max_distance_ptr,
                std::vector<SeqTrie::search_context> & output,
                pairchar_map & cost_map,
                trqwe::simple_progress & progress_bar) : 
    root(root), query(query), max_distance_ptr(max_distance_ptr), output(output), 
    cost_map(cost_map), progress_bar(progress_bar) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t i=begin; i<end; ++i) {
      output[i] = root.levenshtein_search(query[i], max_distance_ptr[i], cost_map);
      progress_bar.increment();
    }
  }
};

struct AnchoredWorkerWithCostMap : public Worker {
  const SeqTrie::RadixTreeR & root;
  const std::vector<cspan> & query;
  int const * const max_distance_ptr;
  std::vector<SeqTrie::search_context> & output;
  pairchar_map & cost_map;
  trqwe::simple_progress & progress_bar;
  AnchoredWorkerWithCostMap(const SeqTrie::RadixTreeR & root,
                const std::vector<cspan> & query,
                int const * const max_distance_ptr,
                std::vector<SeqTrie::search_context> & output,
                pairchar_map & cost_map,
                trqwe::simple_progress & progress_bar) : 
    root(root), query(query), max_distance_ptr(max_distance_ptr), output(output), 
    cost_map(cost_map), progress_bar(progress_bar) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t i=begin; i<end; ++i) {
      output[i] = root.anchored_search(query[i], max_distance_ptr[i], cost_map);
      progress_bar.increment();
    }
  }
};

////////////////////////////////////////////////////////////////////////////////
// RadixTreeContext class definitions 

// [[Rcpp::export(rng = false)]]
double RadixTree_size(RadixTreeRXPtr xp) { return static_cast<double>(xp->size()); }

// [[Rcpp::export(rng = false)]]
LogicalVector RadixTree_insert(RadixTreeRXPtr xp, CharacterVector sequences) {
  auto & root = *xp;
  SEXP * sequence_ptr = STRING_PTR(sequences);
  size_t nseqs = Rf_xlength(sequences);
  LogicalVector result(nseqs);
  int * result_ptr = LOGICAL(result);
  for(size_t i=0; i<nseqs; ++i) {
    cspan sequence = charsxp_to_cspan(sequence_ptr[i]);
    size_t idx = root.insert(sequence, SeqTrie::posidx);
    result_ptr[i] = idx == SeqTrie::nullidx ? 1 : 0; // nullidx means it was successfully inserted
  }
  return result;
}

// [[Rcpp::export(rng = false)]]
LogicalVector RadixTree_erase(RadixTreeRXPtr xp, CharacterVector sequences) {
  auto & root = *xp;
  SEXP * sequence_ptr = STRING_PTR(sequences);
  size_t nseqs = Rf_xlength(sequences);
  LogicalVector result(nseqs);
  int * result_ptr = LOGICAL(result);
  for(size_t i=0; i<nseqs; ++i) {
    cspan sequence = charsxp_to_cspan(sequence_ptr[i]);
    size_t idx = root.erase(sequence);
    result_ptr[i] = idx == SeqTrie::nullidx ? 0 : 1; // nullidx means sequence did not exist, erase was not succesful
  }
  return result;
}

// [[Rcpp::export(rng = false)]]
LogicalVector RadixTree_find(RadixTreeRXPtr xp, CharacterVector sequences) {
  auto & root = *xp;
  SEXP * sequence_ptr = STRING_PTR(sequences);
  size_t nseqs = Rf_xlength(sequences);
  LogicalVector result(nseqs);
  int * result_ptr = LOGICAL(result);
  for(size_t i=0; i<nseqs; ++i) {
    cspan sequence = charsxp_to_cspan(sequence_ptr[i]);
    size_t idx = root.find(sequence);
    result_ptr[i] = idx == SeqTrie::nullidx ? 0 : 1; // nullidx means sequence was not found
  }
  return result;
}

// [[Rcpp::export(rng = false)]]
DataFrame RadixTree_hamming_search(RadixTreeRXPtr xp, CharacterVector sequences, IntegerVector max_distance, Rcpp::Nullable<IntegerMatrix> cost_matrix, const int nthreads, const bool show_progress) {
  auto & root = *xp;
  size_t nseqs = Rf_xlength(sequences);
  SEXP * sequence_ptr = STRING_PTR(sequences);
  int * max_distance_ptr = INTEGER(max_distance);
  
  if(nseqs == 0) {
    return DataFrame::create(_["query"] = CharacterVector(), _["target"] = CharacterVector(), _["distance"] = IntegerVector(), _["stringsAsFactors"] = false);
  }
  
  std::vector<cspan> query(nseqs);
  for(size_t i=0; i<nseqs; ++i) { query[i] = charsxp_to_cspan(sequence_ptr[i]); }
  std::vector<SeqTrie::search_context> output(nseqs);
  trqwe::simple_progress progress_bar(nseqs, show_progress);
  
  if(nthreads == 1) {
    for(size_t i=0; i<nseqs; ++i) {
      output[i] = root.hamming_search(query[i], max_distance_ptr[i]);
      progress_bar.increment();
    }
  } else {
    HammingWorker w(root, query, max_distance_ptr, output, progress_bar);
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


// [[Rcpp::export(rng = false)]]
DataFrame RadixTree_levenshtein_search(RadixTreeRXPtr xp, CharacterVector sequences, IntegerVector max_distance, Rcpp::Nullable<IntegerMatrix> cost_matrix, const int nthreads, const bool show_progress) {
  auto & root = *xp;
  size_t nseqs = Rf_xlength(sequences);
  SEXP * sequence_ptr = STRING_PTR(sequences);
  int * max_distance_ptr = INTEGER(max_distance);
  
  if(nseqs == 0) {
    return DataFrame::create(_["query"] = CharacterVector(), _["target"] = CharacterVector(), _["distance"] = IntegerVector(), _["stringsAsFactors"] = false);
  }
  
  std::vector<cspan> query(nseqs);
  for(size_t i=0; i<nseqs; ++i) { query[i] = charsxp_to_cspan(sequence_ptr[i]); }
  std::vector<SeqTrie::search_context> output(nseqs);
  trqwe::simple_progress progress_bar(nseqs, show_progress);
  
  if(cost_matrix.isNotNull()) {
    IntegerMatrix cost_matrix_(cost_matrix);
    pairchar_map cost_map = convert_cost_matrix(cost_matrix_);
    LevenshteinWorkerWithCostMap w(root, query, max_distance_ptr, output, cost_map, progress_bar);
    parallelFor(0, nseqs, w, 1, nthreads);
  } else {
    LevenshteinWorker w(root, query, max_distance_ptr, output, progress_bar);
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

// [[Rcpp::export(rng = false)]]
DataFrame RadixTree_anchored_search(RadixTreeRXPtr xp, CharacterVector sequences, IntegerVector max_distance, Rcpp::Nullable<IntegerMatrix> cost_matrix, const int nthreads, const bool show_progress) {
  auto & root = *xp;
  size_t nseqs = Rf_xlength(sequences);
  SEXP * sequence_ptr = STRING_PTR(sequences);
  int * max_distance_ptr = INTEGER(max_distance);
  
  if(nseqs == 0) {
    return DataFrame::create(_["query"] = CharacterVector(), _["target"] = CharacterVector(), _["distance"] = IntegerVector(), _["stringsAsFactors"] = false);
  }
  
  std::vector<cspan> query(nseqs);
  for(size_t i=0; i<nseqs; ++i) { query[i] = charsxp_to_cspan(sequence_ptr[i]); }
  std::vector<SeqTrie::search_context> output(nseqs);
  trqwe::simple_progress progress_bar(nseqs, show_progress);
  
  if(cost_matrix.isNotNull()) {
    IntegerMatrix cost_matrix_(cost_matrix);
    pairchar_map cost_map = convert_cost_matrix(cost_matrix_);
    AnchoredWorkerWithCostMap w(root, query, max_distance_ptr, output, cost_map, progress_bar);
    parallelFor(0, nseqs, w, 1, nthreads);
  } else {
    AnchoredWorker w(root, query, max_distance_ptr, output, progress_bar);
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

// [[Rcpp::export(rng = false)]]
DataFrame RadixTree_prefix_search(RadixTreeRXPtr xp, CharacterVector sequences) {
  auto & root = *xp;
  SEXP * sequence_ptr = STRING_PTR(sequences);
  size_t nseqs = Rf_xlength(sequences);
  std::vector<std::vector<SeqTrie::path>> output(nseqs);
  
  if(nseqs == 0) {
    return DataFrame::create(_["query"] = CharacterVector(), _["target"] = CharacterVector(), _["stringsAsFactors"] = false);
  }
  
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

// [[Rcpp::export(rng = false)]]
std::string RadixTree_print(RadixTreeRXPtr xp) {
  auto & root = *xp;
  return root.print();
}

// [[Rcpp::export(rng = false)]]
DataFrame RadixTree_graph(RadixTreeRXPtr xp, const double max_depth) {
  auto & root = *xp;

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

// [[Rcpp::export(rng = false)]]
CharacterVector RadixTree_to_vector(RadixTreeRXPtr xp) {
  auto & root = *xp;
  auto seqs = root.all();
  if(seqs.size() == 0) return R_NilValue;
  CharacterVector sequence(seqs.size());
  for(size_t i=0; i<seqs.size(); ++i) {
    auto s = seqs[i]->template sequence<trqwe::small_array<char>>();
    SET_STRING_ELT(sequence, i, to_charsxp(s));
  }
  return sequence;
}

// [[Rcpp::export(rng = false)]]
bool RadixTree_validate(RadixTreeRXPtr xp) {
  auto & root = *xp;
  return root.validate();
}

// [[Rcpp::export(rng = false)]]
RadixTreeRXPtr RadixTree_create() {
  return RadixTreeRXPtr(new SeqTrie::RadixTreeR, true);
}

