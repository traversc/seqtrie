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
// RadixForest

// [[Rcpp::export(rng = false)]]
double RadixForest_size(RadixForestRXPtr xp) {
  auto & forest = *xp;
  size_t size = 0;
  for(auto & x : forest) {
    size += x.second.size();
  }
  return static_cast<double>(size);
}

// [[Rcpp::export(rng = false)]]
LogicalVector RadixForest_insert(RadixForestRXPtr xp, CharacterVector sequences) {
  auto & forest = *xp;
  SEXP * sequence_ptr = STRING_PTR(sequences);
  size_t nseqs = Rf_xlength(sequences);
  LogicalVector result(nseqs);
  int * result_ptr = LOGICAL(result);
  for(size_t i=0; i<nseqs; ++i) {
    cspan sequence = charsxp_to_cspan(sequence_ptr[i]);
    auto it = forest.find(sequence.size());
    if(it == forest.end()) {
      // the insert has to happen after initialization (emplace)
      // otherwise the address of the root will change and the tree will be invalid because a child will point to the old root
      forest.emplace(sequence.size(), SeqTrie::RadixTreeR());
      forest[sequence.size()].insert(sequence, SeqTrie::posidx);
      result_ptr[i] = 1;
    } else {
      it->second.insert(sequence, SeqTrie::posidx);
      result_ptr[i] = 0;
    }
  }
  return result;
}

// [[Rcpp::export(rng = false)]]
LogicalVector RadixForest_erase(RadixForestRXPtr xp, CharacterVector sequences) {
  auto & forest = *xp;
  SEXP * sequence_ptr = STRING_PTR(sequences);
  size_t nseqs = Rf_xlength(sequences);
  LogicalVector result(nseqs);
  int * result_ptr = LOGICAL(result);
  for(size_t i=0; i<nseqs; ++i) {
    cspan sequence = charsxp_to_cspan(sequence_ptr[i]);
    auto it = forest.find(sequence.size());
    if(it != forest.end()) {
      size_t idx = it->second.erase(sequence);
      result_ptr[i] = idx == SeqTrie::nullidx ? 0 : 1; // nullidx means sequence did not exist, erase was not succesful
      // erase a tree if it is empty
      if(it->second.get_child_nodes().size() == 0) {
        forest.erase(it);
      }
    } else {
      result_ptr[i] = 0;
    }
  }
  return result;
}

// [[Rcpp::export(rng = false)]]
LogicalVector RadixForest_find(RadixForestRXPtr xp, CharacterVector sequences) {
  auto & forest = *xp;
  SEXP * sequence_ptr = STRING_PTR(sequences);
  size_t nseqs = Rf_xlength(sequences);
  LogicalVector result(nseqs);
  int * result_ptr = LOGICAL(result);
  for(size_t i=0; i<nseqs; ++i) {
    cspan sequence = charsxp_to_cspan(sequence_ptr[i]);
    auto it = forest.find(sequence.size());
    if(it != forest.end()) {
      size_t idx = it->second.find(sequence);
      result_ptr[i] = idx == SeqTrie::nullidx ? 0 : 1; // nullidx means sequence was not found
    } else {
      result_ptr[i] = 0;
    }
  }
  return result;
}

// [[Rcpp::export(rng = false)]]
DataFrame RadixForest_prefix_search(RadixForestRXPtr xp, CharacterVector sequences) {
  auto & forest = *xp;
  SEXP * sequence_ptr = STRING_PTR(sequences);
  size_t nseqs = Rf_xlength(sequences);

  std::vector<size_t> queries;
  std::vector<std::vector<SeqTrie::path>> output;  
  for(size_t i=0; i<nseqs; ++i) {
    cspan sequence = charsxp_to_cspan(sequence_ptr[i]);
    for(auto & x : forest) {
      auto res = x.second.prefix_search(sequence);
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
      SET_STRING_ELT(query_results, q, STRING_ELT(sequences, queries[i]));
      auto s = targets[j]->template sequence<trqwe::small_array<char>>();
      SET_STRING_ELT(target_results, q, to_charsxp(s));
      q++;
    }
  }
  return DataFrame::create(_["query"] = query_results, _["target"] = target_results, _["stringsAsFactors"] = false);
}

// [[Rcpp::export(rng = false)]]
std::vector<std::string> RadixForest_print(RadixForestRXPtr xp) {
  auto & forest = *xp;
  std::vector<std::string> output;
  for(auto & x : forest) {
    output.push_back(x.second.print());
  }
  return output;
}

// [[Rcpp::export(rng = false)]]
DataFrame RadixForest_graph(RadixForestRXPtr xp, const double max_depth) {
  auto & forest = *xp;

  size_t depth2;
  if(max_depth < 0) {
    depth2 = -1;
  } else if(max_depth >= static_cast<double>(std::numeric_limits<size_t>::max())) {
    depth2 = -1;
  } else {
    depth2 = static_cast<size_t>(max_depth);
  }

  std::vector<SeqTrie::path> parent_vec;
  std::vector<SeqTrie::path> child_vec;
  for(auto & x : forest) {
    auto seqs = x.second.graph(depth2);
    parent_vec.insert(parent_vec.end(), seqs.first.begin(), seqs.first.end());
    child_vec.insert(child_vec.end(), seqs.second.begin(), seqs.second.end());
  }
  CharacterVector parent(parent_vec.size());
  CharacterVector child(child_vec.size());
  for(size_t i=0; i<parent_vec.size(); ++i) {
    SET_STRING_ELT(parent, i, to_charsxp(parent_vec[i]->get_branch()));
    SET_STRING_ELT(child, i, to_charsxp(child_vec[i]->get_branch()));
  }
  return DataFrame::create(_["parent"] = parent, _["child"] = child, _["stringsAsFactors"] = false);
}

// [[Rcpp::export(rng = false)]]
CharacterVector RadixForest_to_vector(RadixForestRXPtr xp) {
  auto & forest = *xp;
  std::vector<SeqTrie::path> seqs;
  for(auto & x : forest) {
    auto s = x.second.all();
    seqs.insert(seqs.end(), s.begin(), s.end());
  }
  CharacterVector sequence(seqs.size());
  for(size_t i=0; i<seqs.size(); ++i) {
    auto s = seqs[i]->template sequence<trqwe::small_array<char>>();
    SET_STRING_ELT(sequence, i, to_charsxp(s));
  }
  return sequence;
}

// [[Rcpp::export(rng = false)]]
bool RadixForest_validate(RadixForestRXPtr xp) {
  auto & forest = *xp;
  for(auto & x : forest) {
    if(!x.second.validate()) {
      return false;
    }
  }
  return true;
}

// [[Rcpp::export(rng = false)]]
RadixForestRXPtr RadixForest_create() {
  return RadixForestRXPtr(new SeqTrie::RadixForestR, true);
}

// All input parameters should be checked in R, so any error thrown here is an internal error
// [[Rcpp::export(rng = false)]]
DataFrame RadixForest_search(RadixForestRXPtr xp,
                           CharacterVector query,
                           IntegerVector max_distance,
                           const std::string mode = "global", // global or hamming
                           const int nthreads = 1, const bool show_progress = false) {

  auto & forest = *xp;
  size_t nseqs = Rf_xlength(query);
  int * max_distance_ptr = INTEGER(max_distance);
  std::vector<cspan> query_span =  strsxp_to_cspan(query);
  std::vector<SeqTrie::search_context> output(nseqs);
  trqwe::simple_progress progress_bar(nseqs, show_progress);

  if(nseqs == 0) {
    return DataFrame::create(_["query"] = CharacterVector(), _["target"] = CharacterVector(), _["distance"] = IntegerVector(), _["stringsAsFactors"] = false);
  }

  if(mode == "hamming") {
    do_parallel_for([&forest, &query_span, max_distance_ptr, &output, &progress_bar](size_t begin, size_t end) {
      for(size_t i=begin; i<end; ++i) {
        size_t len = query_span[i].size();
        auto it = forest.find(len);
        if(it != forest.end()) {
          output[i] = it->second.hamming_search(query_span[i], max_distance_ptr[i]);
        }
        progress_bar.increment();
      }
    }, 0, nseqs, 1, nthreads);
  } else if(mode == "global") {
    do_parallel_for([&forest, &query_span, max_distance_ptr, &output, &progress_bar](size_t begin, size_t end) {
      for(size_t i=begin; i<end; ++i) {
        size_t len = query_span[i].size();
        size_t min_search_len = len > static_cast<size_t>(max_distance_ptr[i]) ? len - static_cast<size_t>(max_distance_ptr[i]) : 0;
        size_t max_search_len = len + static_cast<size_t>(max_distance_ptr[i]);
        for(size_t j=min_search_len; j<=max_search_len; ++j) {
          auto it = forest.find(j);
          if(it != forest.end()) {
            SeqTrie::search_context res = it->second.global_search(query_span[i], max_distance_ptr[i]);
            output[i].append(res);
          }
        }
        progress_bar.increment();
      }
    }, 0, nseqs, 1, nthreads);
  }
  return seqtrie_results_to_dataframe(query, output);
}
