#include "seqtrie_types.h"
#include "simple_progress/simple_progress.h"

#include <utility>

////////////////////////////////////////////////////////////////////////////////
// RadixTree R functions

// [[Rcpp::export(rng = false)]]
double RadixTree_size(RadixTreeRXPtr xp) {
  return static_cast<double>(xp->root.size());
}

// [[Rcpp::export(rng = false)]]
LogicalVector RadixTree_insert(RadixTreeRXPtr xp, CharacterVector sequences) {
  auto & tree = *xp;
  const SEXP * sequence_ptr = STRING_PTR_RO(sequences);
  size_t nseqs = checked_size_from_r_xlen(Rf_xlength(sequences));
  LogicalVector result(nseqs);
  int * result_ptr = LOGICAL(result);
  for(size_t i=0; i<nseqs; ++i) {
    cspan sequence = charsxp_to_cspan(sequence_ptr[i]);
    const auto terminal_idx = tree.next_terminal_idx++;
    SeqTrie::path p = tree.insert(sequence, terminal_idx);
    const bool inserted = !p.m;
    result_ptr[i] = inserted ? 1 : 0;
  }
  tree.compact();
  return result;
}

// [[Rcpp::export(rng = false)]]
LogicalVector RadixTree_erase(RadixTreeRXPtr xp, CharacterVector sequences) {
  auto & tree = *xp;
  const SEXP * sequence_ptr = STRING_PTR_RO(sequences);
  size_t nseqs = checked_size_from_r_xlen(Rf_xlength(sequences));
  LogicalVector result(nseqs);
  int * result_ptr = LOGICAL(result);
  for(size_t i=0; i<nseqs; ++i) {
    cspan sequence = charsxp_to_cspan(sequence_ptr[i]);
    size_t idx = tree.erase(sequence);
    result_ptr[i] = idx == SeqTrie::nullidx ? 0 : 1; // nullidx means sequence did not exist, erase was not succesful
  }
  tree.compact();
  return result;
}

// [[Rcpp::export(rng = false)]]
LogicalVector RadixTree_find(RadixTreeRXPtr xp, CharacterVector sequences) {
  auto & root = xp->root;
  const SEXP * sequence_ptr = STRING_PTR_RO(sequences);
  size_t nseqs = checked_size_from_r_xlen(Rf_xlength(sequences));
  LogicalVector result(nseqs);
  int * result_ptr = LOGICAL(result);
  for(size_t i=0; i<nseqs; ++i) {
    cspan sequence = charsxp_to_cspan(sequence_ptr[i]);
    SeqTrie::path p = root.find(sequence);
    size_t idx = p.m ? p->get_terminal_idx() : SeqTrie::nullidx;
    result_ptr[i] = idx == SeqTrie::nullidx ? 0 : 1; // nullidx means sequence was not found
  }
  return result;
}


// [[Rcpp::export(rng = false)]]
DataFrame RadixTree_prefix_search(RadixTreeRXPtr xp, CharacterVector sequences) {
  auto & root = xp->root;
  const SEXP * sequence_ptr = STRING_PTR_RO(sequences);
  size_t nseqs = checked_size_from_r_xlen(Rf_xlength(sequences));
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
      SET_STRING_ELT(query_results, checked_r_xlen(q), STRING_ELT(sequences, checked_r_xlen(i)));
      auto s = targets[j]->template sequence<SeqTrie::array_r<char>>();
      SET_STRING_ELT(target_results, checked_r_xlen(q), to_charsxp(s));
      q++;
    }
  }
  return DataFrame::create(_["query"] = query_results, _["target"] = target_results, _["stringsAsFactors"] = false);
}

// [[Rcpp::export(rng = false)]]
std::string RadixTree_print(RadixTreeRXPtr xp) {
  auto & root = xp->root;
  return root.print();
}

// [[Rcpp::export(rng = false)]]
DataFrame RadixTree_graph(RadixTreeRXPtr xp, const double max_depth) {
  auto & root = xp->root;

  size_t depth2;
  if(max_depth < 0) {
    depth2 = std::numeric_limits<size_t>::max();
  } else if(max_depth >= static_cast<double>(std::numeric_limits<size_t>::max())) {
    depth2 = std::numeric_limits<size_t>::max();
  } else {
    depth2 = static_cast<size_t>(max_depth);
  }
  auto seqs = root.graph(depth2);
  if(seqs.first.size() == 0) return R_NilValue;
  CharacterVector parent(seqs.first.size());
  CharacterVector child(seqs.first.size());
  for(size_t i=0; i<seqs.first.size(); ++i) {
    SET_STRING_ELT(parent, checked_r_xlen(i), to_charsxp(seqs.first[i]->get_branch()));
    SET_STRING_ELT(child, checked_r_xlen(i), to_charsxp(seqs.second[i]->get_branch()));
  }
  return DataFrame::create(_["parent"] = parent, _["child"] = child, _["stringsAsFactors"] = false);
}

// [[Rcpp::export(rng = false)]]
CharacterVector RadixTree_to_vector(RadixTreeRXPtr xp) {
  auto & root = xp->root;
  auto seqs = root.all();
  CharacterVector sequence(seqs.size());
  for(size_t i=0; i<seqs.size(); ++i) {
    auto s = seqs[i]->template sequence<SeqTrie::array_r<char>>();
    SET_STRING_ELT(sequence, checked_r_xlen(i), to_charsxp(s));
  }
  return sequence;
}

// [[Rcpp::export(rng = false)]]
bool RadixTree_validate(RadixTreeRXPtr xp) {
  auto & root = xp->root;
  return root.validate();
}

// [[Rcpp::export(rng = false)]]
RadixTreeRXPtr RadixTree_create() {
  auto ptr = std::make_unique<SeqTrie::RadixTreeR>();
  ptr->compact();
  return RadixTreeRXPtr(ptr.release(), true);
}

// All input parameters should be checked in R, so any error thrown here is an internal error
// [[Rcpp::export(rng = false)]]
DataFrame RadixTree_search(RadixTreeRXPtr xp,
                           CharacterVector query,
                           IntegerVector max_distance,
                           const std::string mode = "global", // global, anchored or hamming
                           Rcpp::Nullable<IntegerMatrix> cost_matrix = R_NilValue,
                           int gap_cost = NA_INTEGER,
                           int gap_open_cost = NA_INTEGER,
                           const bool lower_triangle = false,
                           const std::string match_mode = "all",
                           IntegerVector query_index = IntegerVector(),
                           const int nthreads = 1, const bool show_progress = false) {
  auto & root = xp->root;
  size_t nseqs = checked_size_from_r_xlen(Rf_xlength(query));
  const SearchFilterOptions search_options = make_search_filter_options(lower_triangle, match_mode);
  if(checked_size_from_r_xlen(Rf_xlength(query_index)) == 0) {
    query_index = IntegerVector(checked_r_xlen(nseqs));
    for(size_t i = 0; i < nseqs; ++i) {
      query_index[checked_r_xlen(i)] = checked_r_len(SeqTrie::posidx + i);
    }
  } else if(checked_size_from_r_xlen(Rf_xlength(query_index)) != nseqs) {
    throw std::runtime_error("Internal error: query_index length must match query length");
  }
  int * query_index_ptr = INTEGER(query_index);
  int * max_distance_ptr = INTEGER(max_distance);
  std::vector<cspan> query_span =  strsxp_to_cspan(query);
  std::vector<SeqTrie::search_result> output(nseqs);
  trqwe::simple_progress progress_bar(nseqs, show_progress);

  if(nseqs == 0) {
    return DataFrame::create(_["query"] = CharacterVector(), _["target"] = CharacterVector(), _["distance"] = IntegerVector(), _["stringsAsFactors"] = false);
  }
  auto algo = decide_alignment_algo(mode, cost_matrix, gap_cost, gap_open_cost);
  if(algo == AlignmentAlgo::Hamming) {
    do_parallel_for([&root, &query_span, max_distance_ptr, query_index_ptr, &output, &progress_bar, &search_options](size_t begin, size_t end) {
      for(size_t i=begin; i<end; ++i) {
        output[i] = run_search_with_optional_hook(
          search_options,
          query_index_ptr[i],
          [&]() { return root.hamming_search(query_span[i], max_distance_ptr[i]).release_result(); },
          [&](SearchResultFilterHook & hook) { return root.hamming_search(query_span[i], max_distance_ptr[i], hook).release_result(); }
        );
        progress_bar.increment();
      }
    }, 0, nseqs, 1, nthreads);
  } else if(is_global_algo(algo)) {
    if(algo == AlignmentAlgo::GlobalUnit) {
      do_parallel_for([&root, &query_span, max_distance_ptr, query_index_ptr, &output, &progress_bar, &search_options](size_t begin, size_t end) {
        for(size_t i=begin; i<end; ++i) {
          output[i] = run_search_with_optional_hook(
            search_options,
            query_index_ptr[i],
            [&]() { return root.global_search(query_span[i], max_distance_ptr[i]).release_result(); },
            [&](SearchResultFilterHook & hook) { return root.global_search(query_span[i], max_distance_ptr[i], hook).release_result(); }
          );
          progress_bar.increment();
        }
      }, 0, nseqs, 1, nthreads);
    } else if(algo == AlignmentAlgo::GlobalLinear) {
      CostMap cost_map = convert_cost_matrix(cost_matrix.get(), gap_cost, gap_open_cost);
      do_parallel_for([&root, &query_span, max_distance_ptr, query_index_ptr, &output, &cost_map, &progress_bar, &search_options](size_t begin, size_t end) {
        for(size_t i=begin; i<end; ++i) {
          output[i] = run_search_with_optional_hook(
            search_options,
            query_index_ptr[i],
            [&]() { return root.global_search_linear(query_span[i], max_distance_ptr[i], cost_map); },
            [&](SearchResultFilterHook & hook) { return root.global_search_linear(query_span[i], max_distance_ptr[i], cost_map, hook); }
          );
          progress_bar.increment();
        }
      }, 0, nseqs, 1, nthreads);
    } else { // GlobalAffine
      CostMap cost_map = convert_cost_matrix(cost_matrix.get(), gap_cost, gap_open_cost);
      do_parallel_for([&root, &query_span, max_distance_ptr, query_index_ptr, &output, &cost_map, &progress_bar, &search_options](size_t begin, size_t end) {
        for(size_t i=begin; i<end; ++i) {
          output[i] = run_search_with_optional_hook(
            search_options,
            query_index_ptr[i],
            [&]() { return root.global_search_affine(query_span[i], max_distance_ptr[i], cost_map); },
            [&](SearchResultFilterHook & hook) { return root.global_search_affine(query_span[i], max_distance_ptr[i], cost_map, hook); }
          );
          progress_bar.increment();
        }
      }, 0, nseqs, 1, nthreads);
    }
  } else if(is_anchored_algo(algo)) { // anchored
    if(algo == AlignmentAlgo::AnchoredUnit) {
      do_parallel_for([&root, &query_span, max_distance_ptr, query_index_ptr, &output, &progress_bar, &search_options](size_t begin, size_t end) {
        for(size_t i=begin; i<end; ++i) {
          output[i] = run_search_with_optional_hook(
            search_options,
            query_index_ptr[i],
            [&]() { return root.anchored_search(query_span[i], max_distance_ptr[i]).release_result(); },
            [&](SearchResultFilterHook & hook) { return root.anchored_search(query_span[i], max_distance_ptr[i], hook).release_result(); }
          );
          progress_bar.increment();
        }
      }, 0, nseqs, 1, nthreads);
    } else if(algo == AlignmentAlgo::AnchoredLinear) {
      CostMap cost_map = convert_cost_matrix(cost_matrix.get(), gap_cost, gap_open_cost);
      do_parallel_for([&root, &query_span, max_distance_ptr, query_index_ptr, &output, &cost_map, &progress_bar, &search_options](size_t begin, size_t end) {
        for(size_t i=begin; i<end; ++i) {
          output[i] = run_search_with_optional_hook(
            search_options,
            query_index_ptr[i],
            [&]() { return root.anchored_search_linear(query_span[i], max_distance_ptr[i], cost_map); },
            [&](SearchResultFilterHook & hook) { return root.anchored_search_linear(query_span[i], max_distance_ptr[i], cost_map, hook); }
          );
          progress_bar.increment();
        }
      }, 0, nseqs, 1, nthreads);
    } else { // AnchoredAffine
      CostMap cost_map = convert_cost_matrix(cost_matrix.get(), gap_cost, gap_open_cost);
      do_parallel_for([&root, &query_span, max_distance_ptr, query_index_ptr, &output, &cost_map, &progress_bar, &search_options](size_t begin, size_t end) {
        for(size_t i=begin; i<end; ++i) {
          output[i] = run_search_with_optional_hook(
            search_options,
            query_index_ptr[i],
            [&]() { return root.anchored_search_affine(query_span[i], max_distance_ptr[i], cost_map); },
            [&](SearchResultFilterHook & hook) { return root.anchored_search_affine(query_span[i], max_distance_ptr[i], cost_map, hook); }
          );
          progress_bar.increment();
        }
      }, 0, nseqs, 1, nthreads);
    }
  } else {
    throw std::runtime_error("Internal error: unsupported alignment algorithm");
  }
  return seqtrie_results_to_dataframe(query, output);
}

// [[Rcpp::export(rng = false)]]
DataFrame RadixTree_single_gap_search(RadixTreeRXPtr xp,
                                      CharacterVector query,
                                      IntegerVector max_distance,
                                      const int gap_cost = 1,
                                      const int nthreads = 1,
                                      const bool show_progress = false) {
  auto & root = xp->root;
  const size_t nseqs = checked_size_from_r_xlen(Rf_xlength(query));
  if(nseqs == 0) {
    return DataFrame::create(_["query"] = CharacterVector(),
                             _["target"] = CharacterVector(),
                             _["distance"] = IntegerVector(),
                             _["stringsAsFactors"] = false);
  }

  int * max_distance_ptr = INTEGER(max_distance);
  std::vector<cspan> query_span = strsxp_to_cspan(query);
  std::vector<SeqTrie::search_result> output(nseqs);
  trqwe::simple_progress progress_bar(nseqs, show_progress);

  do_parallel_for([&root, &query_span, max_distance_ptr, &output, gap_cost, &progress_bar](size_t begin, size_t end) {
    for(size_t i = begin; i < end; ++i) {
      output[i] = root.single_gap_search(query_span[i], max_distance_ptr[i], gap_cost).release_result();
      progress_bar.increment();
    }
  }, 0, nseqs, 1, nthreads);

  return seqtrie_results_to_dataframe(query, output);
}
