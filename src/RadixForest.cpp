#include "seqtrie_types.h"
#include "simple_progress/simple_progress.h"
#include <utility>

////////////////////////////////////////////////////////////////////////////////
// RadixForest

// [[Rcpp::export(rng = false)]]
double RadixForest_size(RadixForestRXPtr xp) {
  auto & forest = xp->tree_map;
  size_t size = 0;
  for(const auto & x : forest) {
    size += x.second->root.size();
  }
  return static_cast<double>(size);
}

// [[Rcpp::export(rng = false)]]
LogicalVector RadixForest_insert(RadixForestRXPtr xp, CharacterVector sequences) {
  auto & forest_obj = *xp;
  auto & forest = forest_obj.tree_map;
  const SEXP * sequence_ptr = STRING_PTR_RO(sequences);
  size_t nseqs = checked_size_from_r_xlen(Rf_xlength(sequences));
  LogicalVector result(nseqs);
  int * result_ptr = LOGICAL(result);
  std::vector<SeqTrie::RadixTreeR *> touched;
  touched.reserve(nseqs);

  for(size_t i=0; i<nseqs; ++i) {
    cspan sequence = charsxp_to_cspan(sequence_ptr[i]);
    auto it = forest.find(sequence.size());
    if(it == forest.end()) {
      it = forest.emplace(sequence.size(), std::make_unique<SeqTrie::RadixTreeR>()).first;
    }
    SeqTrie::RadixTreeR * tree = it->second.get();
    if(std::find(touched.begin(), touched.end(), tree) == touched.end()) {
      touched.push_back(tree);
    }
    const auto terminal_idx = forest_obj.next_terminal_idx++;
    SeqTrie::path p = tree->insert(sequence, terminal_idx);
    const bool inserted = !p.m;
    result_ptr[i] = inserted ? 1 : 0;
  }

  for(auto * tree : touched) {
    tree->compact();
  }

  return result;
}

// [[Rcpp::export(rng = false)]]
LogicalVector RadixForest_erase(RadixForestRXPtr xp, CharacterVector sequences) {
  auto & forest_obj = *xp;
  const SEXP * sequence_ptr = STRING_PTR_RO(sequences);
  size_t nseqs = checked_size_from_r_xlen(Rf_xlength(sequences));
  LogicalVector result(nseqs);
  int * result_ptr = LOGICAL(result);
  for(size_t i=0; i<nseqs; ++i) {
    cspan sequence = charsxp_to_cspan(sequence_ptr[i]);
    size_t idx = forest_obj.erase(sequence);
    result_ptr[i] = idx == SeqTrie::nullidx ? 0 : 1; // nullidx means sequence did not exist, erase was not succesful
  }
  forest_obj.compact();
  return result;
}

// [[Rcpp::export(rng = false)]]
LogicalVector RadixForest_find(RadixForestRXPtr xp, CharacterVector sequences) {
  auto & forest = xp->tree_map;
  const SEXP * sequence_ptr = STRING_PTR_RO(sequences);
  size_t nseqs = checked_size_from_r_xlen(Rf_xlength(sequences));
  LogicalVector result(nseqs);
  int * result_ptr = LOGICAL(result);
  for(size_t i=0; i<nseqs; ++i) {
    cspan sequence = charsxp_to_cspan(sequence_ptr[i]);
    auto it = forest.find(sequence.size());
    if(it != forest.end()) {
      SeqTrie::path p = it->second->root.find(sequence);
      size_t idx = p.m ? p->get_terminal_idx() : SeqTrie::nullidx;
      result_ptr[i] = idx == SeqTrie::nullidx ? 0 : 1; // nullidx means sequence was not found
    } else {
      result_ptr[i] = 0;
    }
  }
  return result;
}

// [[Rcpp::export(rng = false)]]
DataFrame RadixForest_prefix_search(RadixForestRXPtr xp, CharacterVector sequences) {
  auto & forest = xp->tree_map;
  const SEXP * sequence_ptr = STRING_PTR_RO(sequences);
  size_t nseqs = checked_size_from_r_xlen(Rf_xlength(sequences));

  std::vector<size_t> queries;
  std::vector<std::vector<SeqTrie::path>> output;  
  for(size_t i=0; i<nseqs; ++i) {
    cspan sequence = charsxp_to_cspan(sequence_ptr[i]);
    for(auto & x : forest) {
      auto res = x.second->root.prefix_search(sequence);
      if(res.size() > 0) {
        queries.push_back(i);
        output.push_back(res);
      }
    }
  }

  size_t nresults = 0;
  for(auto & x: output) { nresults += x.size(); }
  CharacterVector query_results(nresults);
  CharacterVector target_results(nresults);
  size_t q = 0;
  for(size_t i=0; i<output.size(); ++i) {
    auto & targets = output[i];
    for(size_t j=0; j<targets.size(); ++j) {
      SET_STRING_ELT(query_results, checked_r_xlen(q), STRING_ELT(sequences, checked_r_xlen(queries[i])));
      auto s = targets[j]->template sequence<SeqTrie::array_r<char>>();
      SET_STRING_ELT(target_results, checked_r_xlen(q), to_charsxp(s));
      q++;
    }
  }
  return DataFrame::create(_["query"] = query_results, _["target"] = target_results, _["stringsAsFactors"] = false);
}

// [[Rcpp::export(rng = false)]]
std::vector<std::string> RadixForest_print(RadixForestRXPtr xp) {
  auto & forest = xp->tree_map;
  std::vector<std::string> output;
  for(auto & x : forest) {
    output.push_back(x.second->root.print());
  }
  return output;
}

// [[Rcpp::export(rng = false)]]
DataFrame RadixForest_graph(RadixForestRXPtr xp, const double max_depth) {
  auto & forest = xp->tree_map;

  size_t depth2;
  if(max_depth < 0) {
    depth2 = std::numeric_limits<size_t>::max();
  } else if(max_depth >= static_cast<double>(std::numeric_limits<size_t>::max())) {
    depth2 = std::numeric_limits<size_t>::max();
  } else {
    depth2 = static_cast<size_t>(max_depth);
  }

  std::vector<SeqTrie::path> parent_vec;
  std::vector<SeqTrie::path> child_vec;
  for(auto & x : forest) {
    auto seqs = x.second->root.graph(depth2);
    parent_vec.insert(parent_vec.end(), seqs.first.begin(), seqs.first.end());
    child_vec.insert(child_vec.end(), seqs.second.begin(), seqs.second.end());
  }
  CharacterVector parent(parent_vec.size());
  CharacterVector child(child_vec.size());
  for(size_t i=0; i<parent_vec.size(); ++i) {
    SET_STRING_ELT(parent, checked_r_xlen(i), to_charsxp(parent_vec[i]->get_branch()));
    SET_STRING_ELT(child, checked_r_xlen(i), to_charsxp(child_vec[i]->get_branch()));
  }
  return DataFrame::create(_["parent"] = parent, _["child"] = child, _["stringsAsFactors"] = false);
}

// [[Rcpp::export(rng = false)]]
CharacterVector RadixForest_to_vector(RadixForestRXPtr xp) {
  auto & forest = xp->tree_map;
  std::vector<SeqTrie::path> seqs;
  for(auto & x : forest) {
    auto s = x.second->root.all();
    seqs.insert(seqs.end(), s.begin(), s.end());
  }
  CharacterVector sequence(seqs.size());
  for(size_t i=0; i<seqs.size(); ++i) {
    auto s = seqs[i]->template sequence<SeqTrie::array_r<char>>();
    SET_STRING_ELT(sequence, checked_r_xlen(i), to_charsxp(s));
  }
  return sequence;
}

// [[Rcpp::export(rng = false)]]
bool RadixForest_validate(RadixForestRXPtr xp) {
  auto & forest = xp->tree_map;
  for(auto & x : forest) {
    if(!x.second->root.validate()) {
      return false;
    }
  }
  return true;
}

// [[Rcpp::export(rng = false)]]
RadixForestRXPtr RadixForest_create() {
  auto ptr = std::make_unique<SeqTrie::RadixForestR>();
  ptr->compact();
  return RadixForestRXPtr(ptr.release(), true);
}

// All input parameters should be checked in R, so any error thrown here is an internal error
// [[Rcpp::export(rng = false)]]
DataFrame RadixForest_search(RadixForestRXPtr xp,
                             CharacterVector query,
                             IntegerVector max_distance,
                             const std::string mode = "global", // global or hamming
                             Rcpp::Nullable<IntegerMatrix> cost_matrix = R_NilValue,
                             int gap_cost = NA_INTEGER,
                             int gap_open_cost = NA_INTEGER,
                             const bool lower_triangle = false,
                             const std::string match_mode = "all",
                             IntegerVector query_index = IntegerVector(),
                             const int nthreads = 1,
                             const bool show_progress = false) {

  auto & forest_obj = *xp;
  const size_t nseqs = checked_size_from_r_xlen(Rf_xlength(query));
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
  std::vector<cspan> query_span = strsxp_to_cspan(query);
  std::vector<SeqTrie::search_result> output(nseqs);
  trqwe::simple_progress progress_bar(nseqs, show_progress);

  if(nseqs == 0) {
    return DataFrame::create(_["query"] = CharacterVector(),
                             _["target"] = CharacterVector(),
                             _["distance"] = IntegerVector(),
                             _["stringsAsFactors"] = false);
  }

  AlignmentAlgo algo = decide_alignment_algo(mode, cost_matrix, gap_cost, gap_open_cost);

  if(algo == AlignmentAlgo::Hamming) {
    do_parallel_for([&forest_obj, &query_span, max_distance_ptr, query_index_ptr, &output, &progress_bar, &search_options](size_t begin, size_t end) {
      for(size_t i = begin; i < end; ++i) {
        output[i] = run_search_with_optional_hook(
          search_options,
          query_index_ptr[i],
          [&]() { return forest_obj.hamming_search(query_span[i], max_distance_ptr[i]); },
          [&](SearchResultFilterHook & hook) { return forest_obj.hamming_search(query_span[i], max_distance_ptr[i], hook); }
        );
        progress_bar.increment();
      }
    }, 0, nseqs, 1, nthreads);
  } else if(algo == AlignmentAlgo::GlobalUnit) {
    do_parallel_for([&forest_obj, &query_span, max_distance_ptr, query_index_ptr, &output, &progress_bar, &search_options](size_t begin, size_t end) {
      for(size_t i = begin; i < end; ++i) {
        output[i] = run_search_with_optional_hook(
          search_options,
          query_index_ptr[i],
          [&]() { return forest_obj.global_search(query_span[i], max_distance_ptr[i]); },
          [&](SearchResultFilterHook & hook) { return forest_obj.global_search(query_span[i], max_distance_ptr[i], hook); }
        );
        progress_bar.increment();
      }
    }, 0, nseqs, 1, nthreads);
  } else if(algo == AlignmentAlgo::GlobalLinear) {
    CostMap cost_map = convert_cost_matrix(cost_matrix.get(), gap_cost, gap_open_cost);
    do_parallel_for([&forest_obj, &query_span, max_distance_ptr, query_index_ptr, &output, &cost_map, &progress_bar, &search_options](size_t begin, size_t end) {
      for(size_t i = begin; i < end; ++i) {
        output[i] = run_search_with_optional_hook(
          search_options,
          query_index_ptr[i],
          [&]() { return forest_obj.global_search_linear(query_span[i], max_distance_ptr[i], cost_map); },
          [&](SearchResultFilterHook & hook) { return forest_obj.global_search_linear(query_span[i], max_distance_ptr[i], cost_map, hook); }
        );
        progress_bar.increment();
      }
    }, 0, nseqs, 1, nthreads);
  } else if(algo == AlignmentAlgo::GlobalAffine) {
    CostMap cost_map = convert_cost_matrix(cost_matrix.get(), gap_cost, gap_open_cost);
    do_parallel_for([&forest_obj, &query_span, max_distance_ptr, query_index_ptr, &output, &cost_map, &progress_bar, &search_options](size_t begin, size_t end) {
      for(size_t i = begin; i < end; ++i) {
        output[i] = run_search_with_optional_hook(
          search_options,
          query_index_ptr[i],
          [&]() { return forest_obj.global_search_affine(query_span[i], max_distance_ptr[i], cost_map); },
          [&](SearchResultFilterHook & hook) { return forest_obj.global_search_affine(query_span[i], max_distance_ptr[i], cost_map, hook); }
        );
        progress_bar.increment();
      }
    }, 0, nseqs, 1, nthreads);
  } else {
    throw std::runtime_error("Internal error: RadixForest only supports hamming and global alignment");
  }

  return seqtrie_results_to_dataframe(query, output);
}
