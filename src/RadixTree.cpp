#include <Rcpp.h>
#include <RcppParallel.h>

#include <set>
#include <memory>
#include <tuple>
#include "seqtrie_types.h"
#include "seqtrie_utils.h"

#include "simple_progress/simple_progress.h"

using namespace Rcpp;
using namespace RcppParallel;

////////////////////////////////////////////////////////////////////////////////
// RadixTree R functions

// [[Rcpp::export(rng = false)]]
double RadixTree_size(RadixTreeRXPtr xp) {
  return static_cast<double>(xp->size());
}

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

// All input parameters should be checked in R, so any error thrown here is an internal error
// [[Rcpp::export(rng = false)]]
DataFrame RadixTree_search(RadixTreeRXPtr xp,
                           CharacterVector query,
                           IntegerVector max_distance,
                           const std::string mode = "global", // global, anchored or hamming
                           const std::string gap_type = "linear", // linear or affine; ignored if hamming or if cost_matrix is NULL
                           Rcpp::Nullable<IntegerMatrix> cost_matrix = R_NilValue,
                           const int nthreads = 1, const bool show_progress = false) {
  auto & root = *xp;
  size_t nseqs = Rf_xlength(query);
  int * max_distance_ptr = INTEGER(max_distance);
  std::vector<cspan> query_span =  strsxp_to_cspan(query);
  std::vector<SeqTrie::search_context> output(nseqs);
  trqwe::simple_progress progress_bar(nseqs, show_progress);

  if(nseqs == 0) {
    return DataFrame::create(_["query"] = CharacterVector(), _["target"] = CharacterVector(), _["distance"] = IntegerVector(), _["stringsAsFactors"] = false);
  }

    if(mode == "hamming") {
    do_parallel_for([&root, &query_span, max_distance_ptr, &output, &progress_bar](size_t begin, size_t end) {
      for(size_t i=begin; i<end; ++i) {
        output[i] = root.hamming_search(query_span[i], max_distance_ptr[i]);
        progress_bar.increment();
      }
    }, 0, nseqs, 1, nthreads);
  } else if(mode == "global") {
    if(cost_matrix.isNull()) {
      do_parallel_for([&root, &query_span, max_distance_ptr, &output, &progress_bar](size_t begin, size_t end) {
        for(size_t i=begin; i<end; ++i) {
          output[i] = root.global_search(query_span[i], max_distance_ptr[i]);
          progress_bar.increment();
        }
      }, 0, nseqs, 1, nthreads);
    } else if(gap_type == "linear") {
      pairchar_map_type cost_map = convert_cost_matrix(cost_matrix.get());
      do_parallel_for([&root, &query_span, max_distance_ptr, &output, &cost_map, &progress_bar](size_t begin, size_t end) {
        for(size_t i=begin; i<end; ++i) {
          output[i] = root.template global_search_linear<pairchar_map_type>(query_span[i], max_distance_ptr[i], cost_map);
          progress_bar.increment();
        }
      }, 0, nseqs, 1, nthreads);
    } else if(gap_type == "affine") {
      pairchar_map_type cost_map = convert_cost_matrix(cost_matrix.get());
      do_parallel_for([&root, &query_span, max_distance_ptr, &output, &cost_map, &progress_bar](size_t begin, size_t end) {
        for(size_t i=begin; i<end; ++i) {
          output[i] = root.template global_search_affine<pairchar_map_type>(query_span[i], max_distance_ptr[i], cost_map);
          progress_bar.increment();
        }
      }, 0, nseqs, 1, nthreads);
    } else {
      throw std::runtime_error("Internal Error: gap_type must be one of linear or affine");
    }
  } else if(mode == "anchored") { // anchored
    if(cost_matrix.isNull()) {
      do_parallel_for([&root, &query_span, max_distance_ptr, &output, &progress_bar](size_t begin, size_t end) {
        for(size_t i=begin; i<end; ++i) {
          output[i] = root.anchored_search(query_span[i], max_distance_ptr[i]);
          progress_bar.increment();
        }
      }, 0, nseqs, 1, nthreads);
    } else if(gap_type == "linear") {
      pairchar_map_type cost_map = convert_cost_matrix(cost_matrix.get());
      do_parallel_for([&root, &query_span, max_distance_ptr, &output, &cost_map, &progress_bar](size_t begin, size_t end) {
        for(size_t i=begin; i<end; ++i) {
          output[i] = root.template anchored_search_linear<pairchar_map_type>(query_span[i], max_distance_ptr[i], cost_map);
          progress_bar.increment();
        }
      }, 0, nseqs, 1, nthreads);
    } else if(gap_type == "affine") {
      pairchar_map_type cost_map = convert_cost_matrix(cost_matrix.get());
      do_parallel_for([&root, &query_span, max_distance_ptr, &output, &cost_map, &progress_bar](size_t begin, size_t end) {
        for(size_t i=begin; i<end; ++i) {
          output[i] = root.template anchored_search_affine<pairchar_map_type>(query_span[i], max_distance_ptr[i], cost_map);
          progress_bar.increment();
        }
      }, 0, nseqs, 1, nthreads);
    } else {
      throw std::runtime_error("Internal Error: gap_type must be one of linear or affine");
    }
  } else {
    throw std::runtime_error("Internal Error: mode must be one of global, anchored or hamming");
  }
  return seqtrie_results_to_dataframe(query, output);
}
