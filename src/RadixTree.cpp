#include "seqtrie_types.h"
#include "simple_progress/simple_progress.h"

////////////////////////////////////////////////////////////////////////////////
// RadixTree R functions

// [[Rcpp::export(rng = false)]]
double RadixTree_size(RadixTreeRXPtr xp) {
  return static_cast<double>(xp->size());
}

// [[Rcpp::export(rng = false)]]
LogicalVector RadixTree_insert(RadixTreeRXPtr xp, CharacterVector sequences) {
  auto & root = *xp;
  const SEXP * sequence_ptr = STRING_PTR_RO(sequences);
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
  const SEXP * sequence_ptr = STRING_PTR_RO(sequences);
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
  const SEXP * sequence_ptr = STRING_PTR_RO(sequences);
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
  const SEXP * sequence_ptr = STRING_PTR_RO(sequences);
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
      auto s = targets[j]->template sequence<SeqTrie::array_r<char>>();
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
  } else if(max_depth >= static_cast<double>(std::numeric_limits<size_t>::max())) {
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
    auto s = seqs[i]->template sequence<SeqTrie::array_r<char>>();
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
                           Rcpp::Nullable<IntegerMatrix> cost_matrix = R_NilValue,
                           const int gap_cost = 1,
                           const int gap_open_cost = 0,
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
  auto algo = decide_alignment_algo(mode, cost_matrix, gap_cost, gap_open_cost);
  if(algo == AlignmentAlgo::Hamming) {
    do_parallel_for([&root, &query_span, max_distance_ptr, &output, &progress_bar](size_t begin, size_t end) {
      for(size_t i=begin; i<end; ++i) {
        output[i] = root.hamming_search(query_span[i], max_distance_ptr[i]);
        progress_bar.increment();
      }
    }, 0, nseqs, 1, nthreads);
  } else if(mode == "global" || mode == "gb" || mode == "lv" || mode == "levenshtein") {
    if(algo == AlignmentAlgo::GlobalUnit) {
      do_parallel_for([&root, &query_span, max_distance_ptr, &output, &progress_bar](size_t begin, size_t end) {
        for(size_t i=begin; i<end; ++i) {
          output[i] = root.global_search(query_span[i], max_distance_ptr[i]);
          progress_bar.increment();
        }
      }, 0, nseqs, 1, nthreads);
    } else if(algo == AlignmentAlgo::GlobalLinear) {
      CostMap cost_map = convert_cost_matrix(cost_matrix.get(), gap_cost, gap_open_cost);
      do_parallel_for([&root, &query_span, max_distance_ptr, &output, &cost_map, &progress_bar](size_t begin, size_t end) {
        for(size_t i=begin; i<end; ++i) {
          output[i] = root.global_search_linear(query_span[i], max_distance_ptr[i], cost_map);
          progress_bar.increment();
        }
      }, 0, nseqs, 1, nthreads);
    } else { // GlobalAffine
      CostMap cost_map = convert_cost_matrix(cost_matrix.get(), gap_cost, gap_open_cost);
      do_parallel_for([&root, &query_span, max_distance_ptr, &output, &cost_map, &progress_bar](size_t begin, size_t end) {
        for(size_t i=begin; i<end; ++i) {
          output[i] = root.global_search_affine(query_span[i], max_distance_ptr[i], cost_map);
          progress_bar.increment();
        }
      }, 0, nseqs, 1, nthreads);
    }
  } else { // anchored
    if(algo == AlignmentAlgo::AnchoredUnit) {
      do_parallel_for([&root, &query_span, max_distance_ptr, &output, &progress_bar](size_t begin, size_t end) {
        for(size_t i=begin; i<end; ++i) {
          output[i] = root.anchored_search(query_span[i], max_distance_ptr[i]);
          progress_bar.increment();
        }
      }, 0, nseqs, 1, nthreads);
    } else if(algo == AlignmentAlgo::AnchoredLinear) {
      CostMap cost_map = convert_cost_matrix(cost_matrix.get(), gap_cost, gap_open_cost);
      do_parallel_for([&root, &query_span, max_distance_ptr, &output, &cost_map, &progress_bar](size_t begin, size_t end) {
        for(size_t i=begin; i<end; ++i) {
          output[i] = root.anchored_search_linear(query_span[i], max_distance_ptr[i], cost_map);
          progress_bar.increment();
        }
      }, 0, nseqs, 1, nthreads);
    } else { // AnchoredAffine
      CostMap cost_map = convert_cost_matrix(cost_matrix.get(), gap_cost, gap_open_cost);
      do_parallel_for([&root, &query_span, max_distance_ptr, &output, &cost_map, &progress_bar](size_t begin, size_t end) {
        for(size_t i=begin; i<end; ++i) {
          output[i] = root.anchored_search_affine(query_span[i], max_distance_ptr[i], cost_map);
          progress_bar.increment();
        }
      }, 0, nseqs, 1, nthreads);
    }
  }
  return seqtrie_results_to_dataframe(query, output);
}
