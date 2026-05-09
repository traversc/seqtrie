#include "seqtrie_types.h"
#include "simple_progress/simple_progress.h"

#include <memory>
#include <string>
#include <unordered_map>

namespace {

struct SplitParts {
  std::string left;
  std::string right;
};

struct SecondaryIndex {
  SeqTrie::RadixTreeR tree;
  std::vector<std::vector<size_t>> target_indices;
};

struct SplitIndex {
  bool primary_is_left = true;
  SeqTrie::RadixTreeR primary_tree;
  std::vector<std::unique_ptr<SecondaryIndex>> secondary_indices;
};

struct SplitMatch {
  size_t target_index;
  int distance;
};

inline cspan string_span(const std::string & x) {
  return cspan(x.data(), x.size());
}

SplitParts split_sequence(SEXP sequence, int split, int edge_trim) {
  const char * ptr = CHAR(sequence);
  const size_t len = static_cast<size_t>(Rf_xlength(sequence));
  const size_t trim = static_cast<size_t>(edge_trim);
  const size_t split_pos = static_cast<size_t>(split);
  const size_t right_end = len - trim;

  SplitParts out;
  out.left.reserve(split_pos - trim);
  for(size_t i = split_pos; i > trim; --i) {
    out.left.push_back(ptr[i - 1]);
  }
  out.right.assign(ptr + split_pos, ptr + right_end);
  return out;
}

SeqTrie::search_context anchored_search_dispatch(const SeqTrie::RadixTreeR & tree,
                                                 const cspan query,
                                                 const int max_distance,
                                                 const AlignmentAlgo algo,
                                                 const CostMap & cost_map) {
  if(algo == AlignmentAlgo::AnchoredUnit) {
    return tree.anchored_search(query, max_distance);
  } else if(algo == AlignmentAlgo::AnchoredLinear) {
    return tree.anchored_search_linear(query, max_distance, cost_map);
  } else if(algo == AlignmentAlgo::AnchoredAffine) {
    return tree.anchored_search_affine(query, max_distance, cost_map);
  }
  throw std::runtime_error("Internal error: split_search only supports anchored alignment");
}

void insert_target(SplitIndex & index,
                   const SplitParts & parts,
                   const size_t target_index) {
  const std::string & primary = index.primary_is_left ? parts.left : parts.right;
  const std::string & secondary = index.primary_is_left ? parts.right : parts.left;

  const size_t next_primary_id = index.secondary_indices.size();
  SeqTrie::path primary_path = index.primary_tree.insert_get_path(string_span(primary), next_primary_id);
  const size_t primary_id = primary_path->get_terminal_idx();

  if(primary_id == next_primary_id) {
    index.secondary_indices.push_back(std::make_unique<SecondaryIndex>());
  }

  SecondaryIndex & secondary_index = *index.secondary_indices[primary_id];
  const size_t next_secondary_id = secondary_index.target_indices.size();
  SeqTrie::path secondary_path = secondary_index.tree.insert_get_path(string_span(secondary), next_secondary_id);
  const size_t secondary_id = secondary_path->get_terminal_idx();

  if(secondary_id == next_secondary_id) {
    secondary_index.target_indices.emplace_back();
  }
  secondary_index.target_indices[secondary_id].push_back(target_index);
}

std::vector<SplitMatch> search_query(const SplitIndex & index,
                                     const SplitParts & parts,
                                     const int max_distance,
                                     const AlignmentAlgo algo,
                                     const CostMap & cost_map) {
  const std::string & primary = index.primary_is_left ? parts.left : parts.right;
  const std::string & secondary = index.primary_is_left ? parts.right : parts.left;
  std::vector<SplitMatch> out;

  SeqTrie::search_context primary_matches = anchored_search_dispatch(
    index.primary_tree,
    string_span(primary),
    max_distance,
    algo,
    cost_map
  );

  for(size_t i = 0; i < primary_matches.match.size(); ++i) {
    const int primary_distance = primary_matches.distance[i];
    const int remaining_distance = max_distance - primary_distance;
    if(remaining_distance < 0) continue;

    const size_t primary_id = primary_matches.match[i]->get_terminal_idx();
    if(primary_id >= index.secondary_indices.size()) {
      throw std::runtime_error("Internal error: invalid split_search primary index");
    }

    const SecondaryIndex & secondary_index = *index.secondary_indices[primary_id];
    SeqTrie::search_context secondary_matches = anchored_search_dispatch(
      secondary_index.tree,
      string_span(secondary),
      remaining_distance,
      algo,
      cost_map
    );

    for(size_t j = 0; j < secondary_matches.match.size(); ++j) {
      const size_t secondary_id = secondary_matches.match[j]->get_terminal_idx();
      if(secondary_id >= secondary_index.target_indices.size()) {
        throw std::runtime_error("Internal error: invalid split_search secondary index");
      }
      const int distance = primary_distance + secondary_matches.distance[j];
      for(size_t target_index : secondary_index.target_indices[secondary_id]) {
        out.push_back({target_index, distance});
      }
    }
  }
  return out;
}

} // namespace

// [[Rcpp::export(rng = false)]]
DataFrame c_split_search(CharacterVector query,
                         CharacterVector target,
                         IntegerVector query_split,
                         IntegerVector target_split,
                         const int edge_trim = 0,
                         IntegerVector max_distance = IntegerVector(),
                         Rcpp::Nullable<IntegerMatrix> cost_matrix = R_NilValue,
                         int gap_cost = NA_INTEGER,
                         int gap_open_cost = NA_INTEGER,
                         const int nthreads = 1,
                         const bool show_progress = false) {
  const size_t nquery = Rf_xlength(query);
  const size_t ntarget = Rf_xlength(target);

  if(nquery == 0 || ntarget == 0) {
    return DataFrame::create(_["query"] = CharacterVector(),
                             _["target"] = CharacterVector(),
                             _["distance"] = IntegerVector(),
                             _["stringsAsFactors"] = false);
  }

  AlignmentAlgo algo = decide_alignment_algo("anchored", cost_matrix, gap_cost, gap_open_cost);
  CostMap cost_map;
  if(algo == AlignmentAlgo::AnchoredLinear || algo == AlignmentAlgo::AnchoredAffine) {
    cost_map = convert_cost_matrix(cost_matrix.get(), gap_cost, gap_open_cost);
  }

  const SEXP * query_ptr = STRING_PTR_RO(query);
  const SEXP * target_ptr = STRING_PTR_RO(target);
  const int * query_split_ptr = INTEGER(query_split);
  const int * target_split_ptr = INTEGER(target_split);
  const int * max_distance_ptr = INTEGER(max_distance);

  std::unordered_map<std::string, size_t> unique_target_map;
  std::vector<size_t> unique_target_indices;
  std::vector<SplitParts> target_parts;
  size_t total_left_size = 0;
  size_t total_right_size = 0;

  unique_target_indices.reserve(ntarget);
  target_parts.reserve(ntarget);
  for(size_t i = 0; i < ntarget; ++i) {
    std::string target_string(CHAR(target_ptr[i]), static_cast<size_t>(Rf_xlength(target_ptr[i])));
    if(unique_target_map.find(target_string) != unique_target_map.end()) continue;

    const size_t target_index = unique_target_indices.size();
    unique_target_map.emplace(std::move(target_string), target_index);
    unique_target_indices.push_back(i);
    target_parts.push_back(split_sequence(target_ptr[i], target_split_ptr[i], edge_trim));
    total_left_size += target_parts.back().left.size();
    total_right_size += target_parts.back().right.size();
  }

  if(unique_target_indices.empty()) {
    return DataFrame::create(_["query"] = CharacterVector(),
                             _["target"] = CharacterVector(),
                             _["distance"] = IntegerVector(),
                             _["stringsAsFactors"] = false);
  }

  SplitIndex index;
  index.primary_is_left = total_left_size >= total_right_size;
  for(size_t i = 0; i < target_parts.size(); ++i) {
    insert_target(index, target_parts[i], i);
  }

  std::vector<SplitParts> query_parts;
  query_parts.reserve(nquery);
  for(size_t i = 0; i < nquery; ++i) {
    query_parts.push_back(split_sequence(query_ptr[i], query_split_ptr[i], edge_trim));
  }

  std::vector<std::vector<SplitMatch>> output(nquery);
  trqwe::simple_progress progress_bar(nquery, show_progress);
  do_parallel_for([&index, &query_parts, max_distance_ptr, algo, &cost_map, &output, &progress_bar](size_t begin, size_t end) {
    for(size_t i = begin; i < end; ++i) {
      output[i] = search_query(index, query_parts[i], max_distance_ptr[i], algo, cost_map);
      progress_bar.increment();
    }
  }, 0, nquery, 1, nthreads);

  size_t nresults = 0;
  for(const auto & matches : output) nresults += matches.size();

  CharacterVector query_results(nresults);
  CharacterVector target_results(nresults);
  IntegerVector distance_results(nresults);
  int * distance_ptr = INTEGER(distance_results);

  size_t k = 0;
  for(size_t i = 0; i < output.size(); ++i) {
    for(const SplitMatch & match : output[i]) {
      SET_STRING_ELT(query_results, k, STRING_ELT(query, i));
      SET_STRING_ELT(target_results, k, STRING_ELT(target, unique_target_indices[match.target_index]));
      distance_ptr[k] = match.distance;
      ++k;
    }
  }

  return DataFrame::create(_["query"] = query_results,
                           _["target"] = target_results,
                           _["distance"] = distance_results,
                           _["stringsAsFactors"] = false);
}
