#include "seqtrie_types.h"

#include <algorithm>
#include <memory>

////////////////////////////////////////////////////////////////////////////////
// StarTree R functions (global/Levenshtein, Hamming, and anchored modes)

namespace {

std::vector<std::string> character_vector_to_strings(CharacterVector x) {
  const SEXP * ptr = STRING_PTR_RO(x);
  const size_t n = checked_size_from_r_xlen(Rf_xlength(x));
  std::vector<std::string> out;
  out.reserve(n);
  for(size_t i = 0; i < n; ++i) {
    if(ptr[i] == NA_STRING) {
      throw std::runtime_error("StarTree does not accept missing sequences");
    }
    out.emplace_back(CHAR(ptr[i]), checked_size_from_r_xlen(Rf_xlength(ptr[i])));
  }
  return out;
}

void validate_star_sequence(const std::string& seq) {
  if(seq.empty()) {
    throw std::runtime_error("StarTree does not support empty sequences");
  }
  if(seq.size() > static_cast<size_t>(startree::kMaxSeqLen)) {
    throw std::runtime_error("StarTree sequence length cannot exceed " +
                             std::to_string(startree::kMaxSeqLen));
  }
  for(char ch : seq) {
    if(!startree::is_valid_dna(static_cast<unsigned char>(ch))) {
      throw std::runtime_error("StarTree only supports DNA characters A, C, G, T, and N");
    }
  }
}

void validate_star_sequences(const std::vector<std::string>& sequences,
                             const char* label = "sequences") {
  if(sequences.empty()) {
    throw std::runtime_error(std::string(label) + " must contain at least one sequence");
  }
  if(sequences.size() >= static_cast<size_t>(std::numeric_limits<uint32_t>::max())) {
    throw std::runtime_error("too many sequences for StarTree");
  }
  for(const std::string& seq : sequences) {
    validate_star_sequence(seq);
  }
}

void validate_star_params(const int max_distance,
                          const int mismatch_cost,
                          const int gap_cost) {
  if(max_distance < 0) {
    throw std::runtime_error("max_distance must be non-negative");
  }
  if(max_distance == std::numeric_limits<int>::max()) {
    throw std::runtime_error("max_distance must be less than the maximum R integer");
  }
  if(mismatch_cost <= 0) {
    throw std::runtime_error("mismatch_cost must be a positive integer");
  }
  if(gap_cost <= 0) {
    throw std::runtime_error("gap_cost must be a positive integer");
  }

  const startree::SearchParams params =
    startree::make_search_params(max_distance, mismatch_cost, gap_cost, false);
  if(params.lookup_distance > startree::kMaxTau ||
     params.band_radius > startree::kMaxTau) {
    throw std::runtime_error(
      "StarTree requires max_distance / min(mismatch_cost, gap_cost) <= " +
      std::to_string(startree::kMaxTau)
    );
  }
}

int max_string_length(const std::vector<std::string>& x) {
  size_t out = 0;
  for(const auto& value : x) {
    out = std::max(out, value.size());
  }
  if(out > static_cast<size_t>(startree::kMaxSeqLen)) {
    throw std::runtime_error("StarTree sequence length cannot exceed " +
                             std::to_string(startree::kMaxSeqLen));
  }
  if(out > static_cast<size_t>(std::numeric_limits<int>::max())) {
    throw std::runtime_error("sequence length exceeds integer range");
  }
  return static_cast<int>(out);
}

SEXP string_to_charsxp(const std::string& x) {
  return Rf_mkCharLen(x.data(), checked_r_len(x.size()));
}

DataFrame startree_empty_match_result() {
  return DataFrame::create(_["query"] = CharacterVector(),
                           _["target"] = CharacterVector(),
                           _["distance"] = IntegerVector(),
                           _["stringsAsFactors"] = false);
}

template <typename InputDataT>
DataFrame startree_pairs_to_dataframe(const std::vector<startree::PairRecord>& pairs,
                                      const InputDataT& target_data,
                                      const InputDataT& query_data) {
  if(pairs.empty()) {
    return startree_empty_match_result();
  }

  const size_t n = pairs.size();
  CharacterVector query(n);
  CharacterVector target(n);
  IntegerVector distance(n);

  for(size_t i = 0; i < n; ++i) {
    const startree::PairRecord& pair = pairs[i];
    SET_STRING_ELT(query, checked_r_xlen(i),
                   string_to_charsxp(query_data.seqs[pair.query_index].seq));
    SET_STRING_ELT(target, checked_r_xlen(i),
                   string_to_charsxp(target_data.seqs[pair.target_index].seq));
    distance[checked_r_xlen(i)] = pair.distance;
  }

  return DataFrame::create(_["query"] = query,
                           _["target"] = target,
                           _["distance"] = distance,
                           _["stringsAsFactors"] = false);
}

std::vector<startree::PairRecord> flatten_pair_chunks(
    std::vector<std::vector<startree::PairRecord>>& chunks) {
  size_t total = 0;
  for(const auto& chunk : chunks) {
    total += chunk.size();
  }

  std::vector<startree::PairRecord> out;
  out.reserve(total);
  for(auto& chunk : chunks) {
    out.insert(out.end(),
               std::make_move_iterator(chunk.begin()),
               std::make_move_iterator(chunk.end()));
  }
  return out;
}

std::vector<startree::PairRecord> run_self_similarity(const startree::InputData& data,
                                                      const startree::SearchParams& params,
                                                      const int nthreads,
                                                      const bool hamming = false) {
  const std::vector<startree::Block> blocks = startree::make_blocks(data.seqs, nthreads);
  std::vector<std::vector<startree::PairRecord>> per_block(blocks.size());

  // Parallelize over whole target blocks. A block builds its trie incrementally
  // and cannot be split across threads, so each block is one unit of work
  // (grainSize 1 never subdivides a block). Each block writes its own output
  // slot, so no synchronization is needed; the slots are concatenated after.
  do_parallel_for([&](std::size_t begin, std::size_t end) {
    for(std::size_t block_id = begin; block_id < end; ++block_id) {
      if(hamming) {
        startree::hamming::process_target_block_self(data.seqs,
                                                     blocks,
                                                     block_id,
                                                     data.height,
                                                     params,
                                                     &per_block[block_id]);
      } else {
        startree::process_target_block_self(data.seqs,
                                            blocks,
                                            block_id,
                                            data.height,
                                            data.median_len,
                                            params,
                                            &per_block[block_id]);
      }
    }
  }, 0, blocks.size(), 1, nthreads);

  return flatten_pair_chunks(per_block);
}

std::vector<startree::PairRecord> run_query_search(const startree::InputData& target_data,
                                                   const startree::InputData& query_data,
                                                   const startree::SearchParams& params,
                                                   const int nthreads,
                                                   const bool hamming = false) {
  const std::vector<startree::Block> target_blocks = startree::make_blocks(target_data.seqs, nthreads);
  const std::vector<startree::Block> query_blocks = startree::make_blocks(query_data.seqs, nthreads);
  std::vector<std::vector<startree::PairRecord>> per_block(target_blocks.size());

  // Parallelize over target blocks (see run_self_similarity): one block per
  // unit of work, each writing its own output slot, concatenated afterwards.
  do_parallel_for([&](std::size_t begin, std::size_t end) {
    for(std::size_t block_id = begin; block_id < end; ++block_id) {
      if(hamming) {
        startree::hamming::process_target_block_query(target_data.seqs,
                                                      query_data.seqs,
                                                      query_blocks,
                                                      target_blocks[block_id],
                                                      target_data.height,
                                                      params,
                                                      &per_block[block_id]);
      } else {
        startree::process_target_block_query(target_data.seqs,
                                             query_data.seqs,
                                             query_blocks,
                                             target_blocks[block_id],
                                             target_data.height,
                                             target_data.median_len,
                                             params,
                                             &per_block[block_id]);
      }
    }
  }, 0, target_blocks.size(), 1, nthreads);

  return flatten_pair_chunks(per_block);
}

// ---- anchored mode (semi-global) helpers ----

std::vector<std::pair<size_t, size_t>> make_lcp_blocks(const size_t n,
                                                       const int nthreads) {
  const size_t workers = static_cast<size_t>(std::max(1, nthreads));
  const size_t target_block =
    std::max<size_t>(256, n / (workers * 8 + 1));
  std::vector<std::pair<size_t, size_t>> blocks;
  for(size_t begin = 0; begin < n; begin += target_block) {
    blocks.emplace_back(begin, std::min(n, begin + target_block));
  }
  return blocks;
}

template <int MaxBand, typename TrieType, typename WorkspaceType>
std::vector<startree::PairRecord> run_anchored_lcp_search_impl(
    const startree::AnchoredInputData& target_data,
    const startree::AnchoredInputData& query_data,
    const startree::SearchParams& params,
    const int nthreads,
    const bool lower_triangle) {
  const TrieType trie =
    startree::anchored::build_trie<TrieType>(target_data.target_codes);
  const auto blocks = make_lcp_blocks(query_data.query_codes.size(), nthreads);
  std::vector<std::vector<startree::PairRecord>> per_block(blocks.size());

  do_parallel_for([&](std::size_t begin, std::size_t end) {
    for(std::size_t block_id = begin; block_id < end; ++block_id) {
      WorkspaceType ws;
      trie.search_batch_lcp(query_data.query_codes,
                            params,
                            blocks[block_id].first,
                            blocks[block_id].second,
                            &per_block[block_id],
                            &ws,
                            lower_triangle);
    }
  }, 0, blocks.size(), 1, nthreads);

  return flatten_pair_chunks(per_block);
}

template <int MaxBand>
std::vector<startree::PairRecord> run_anchored_lcp_search_scored(
    const startree::AnchoredInputData& target_data,
    const startree::AnchoredInputData& query_data,
    const startree::SearchParams& params,
    const int nthreads,
    const bool lower_triangle) {
  if(params.custom_cost) {
    return run_anchored_lcp_search_impl<
      MaxBand,
      startree::radix::AnchoredCustomTrie<MaxBand>,
      startree::radix::AnchoredCustomWorkspace<MaxBand>
    >(target_data, query_data, params, nthreads, lower_triangle);
  }

  return run_anchored_lcp_search_impl<
    MaxBand,
    startree::radix::AnchoredStandardTrie<MaxBand>,
    startree::radix::AnchoredStandardWorkspace<MaxBand>
  >(target_data, query_data, params, nthreads, lower_triangle);
}

std::vector<startree::PairRecord> run_anchored_lcp_search(
    const startree::AnchoredInputData& target_data,
    const startree::AnchoredInputData& query_data,
    const startree::SearchParams& params,
    const int nthreads,
    const bool lower_triangle) {
  return startree::anchored::dispatch_band(params, [&](auto band) {
    return run_anchored_lcp_search_scored<decltype(band)::value>(
      target_data, query_data, params, nthreads, lower_triangle
    );
  });
}

} // namespace

// [[Rcpp::export(rng = false)]]
StarTreeRXPtr StarTree_create(CharacterVector sequences,
                              const int max_distance,
                              const int mismatch_cost = 1,
                              const int gap_cost = 1,
                              const int nthreads = 1,
                              const bool show_progress = false,
                              const bool hamming = false) {
  if(show_progress) {
    Rcpp::warning("show_progress is not currently implemented for StarTree");
  }
  if(nthreads < 1) {
    throw std::runtime_error("nthreads must be a single positive integer");
  }

  auto raw_sequences = character_vector_to_strings(sequences);
  validate_star_sequences(raw_sequences);
  validate_star_params(max_distance, mismatch_cost, gap_cost);
  auto ptr = std::make_unique<StarTreeR>();
  ptr->data = startree::make_input_data(std::move(raw_sequences), true);
  ptr->sequences.reserve(ptr->data.seqs.size());
  for(const auto& seq : ptr->data.seqs) {
    ptr->sequences.push_back(seq.seq);
  }
  ptr->params = startree::make_search_params(max_distance, mismatch_cost, gap_cost, false);
  ptr->nthreads = nthreads;
  ptr->hamming = hamming;
  ptr->self_pairs = run_self_similarity(ptr->data, ptr->params, nthreads, hamming);
  return StarTreeRXPtr(ptr.release(), true);
}

// [[Rcpp::export(rng = false)]]
DataFrame StarTree_self_search(CharacterVector sequences,
                               const int max_distance,
                               const int mismatch_cost = 1,
                               const int gap_cost = 1,
                               const int nthreads = 1,
                               const bool show_progress = false,
                               const bool hamming = false) {
  if(show_progress) {
    Rcpp::warning("show_progress is not currently implemented for StarTree");
  }
  if(nthreads < 1) {
    throw std::runtime_error("nthreads must be a single positive integer");
  }

  auto raw_sequences = character_vector_to_strings(sequences);
  validate_star_sequences(raw_sequences);
  validate_star_params(max_distance, mismatch_cost, gap_cost);
  const startree::InputData data = startree::make_input_data(std::move(raw_sequences), true);
  const startree::SearchParams params =
    startree::make_search_params(max_distance, mismatch_cost, gap_cost, false);
  const auto pairs = run_self_similarity(data, params, nthreads, hamming);
  return startree_pairs_to_dataframe(pairs, data, data);
}

// [[Rcpp::export(rng = false)]]
double StarTree_size(StarTreeRXPtr xp) {
  return static_cast<double>(xp->data.seqs.size());
}

// [[Rcpp::export(rng = false)]]
CharacterVector StarTree_to_vector(StarTreeRXPtr xp) {
  CharacterVector out(xp->data.seqs.size());
  for(size_t i = 0; i < xp->data.seqs.size(); ++i) {
    SET_STRING_ELT(out, checked_r_xlen(i), string_to_charsxp(xp->data.seqs[i].seq));
  }
  return out;
}

// [[Rcpp::export(rng = false)]]
DataFrame StarTree_result(StarTreeRXPtr xp) {
  return startree_pairs_to_dataframe(xp->self_pairs, xp->data, xp->data);
}

// [[Rcpp::export(rng = false)]]
DataFrame StarTree_search(StarTreeRXPtr xp,
                          CharacterVector query,
                          const int nthreads = 1,
                          const bool show_progress = false) {
  if(show_progress) {
    Rcpp::warning("show_progress is not currently implemented for StarTree");
  }
  if(nthreads < 1) {
    throw std::runtime_error("nthreads must be a single positive integer");
  }
  if(Rf_xlength(query) == 0) {
    return startree_empty_match_result();
  }

  auto raw_query = character_vector_to_strings(query);
  if(raw_query.size() >= static_cast<size_t>(std::numeric_limits<uint32_t>::max())) {
    throw std::runtime_error("too many query sequences for StarTree");
  }
  for(const std::string& seq : raw_query) {
    validate_star_sequence(seq);
  }
  const int common_height = std::max(xp->data.height, max_string_length(raw_query));
  startree::InputData target_data =
    common_height == xp->data.height ?
      xp->data :
      startree::make_input_data(xp->sequences, true, common_height);
  startree::InputData query_data =
    startree::make_input_data(std::move(raw_query), false, common_height);

  startree::SearchParams params = xp->params;
  params.include_zero = true;
  const auto pairs = run_query_search(target_data, query_data, params, nthreads, xp->hamming);
  return startree_pairs_to_dataframe(pairs, target_data, query_data);
}

////////////////////////////////////////////////////////////////////////////////
// AnchoredStarTree R functions (star_tree mode == "anchored")

// [[Rcpp::export(rng = false)]]
AnchoredStarTreeRXPtr AnchoredStarTree_create(CharacterVector sequences,
                                              const int max_distance,
                                              const int mismatch_cost = 1,
                                              const int gap_cost = 1,
                                              const int nthreads = 1,
                                              const bool show_progress = false) {
  if(show_progress) {
    Rcpp::warning("show_progress is not currently implemented for StarTree anchored mode");
  }
  if(nthreads < 1) {
    throw std::runtime_error("nthreads must be a single positive integer");
  }

  auto raw_sequences = character_vector_to_strings(sequences);
  validate_star_sequences(raw_sequences, "sequences");
  validate_star_params(max_distance, mismatch_cost, gap_cost);
  auto ptr = std::make_unique<AnchoredStarTreeR>();
  ptr->data = startree::make_anchored_input_data(std::move(raw_sequences), true);
  ptr->sequences.reserve(ptr->data.seqs.size());
  for(const auto& seq : ptr->data.seqs) {
    ptr->sequences.push_back(seq.seq);
  }
  ptr->params = startree::make_search_params(max_distance, mismatch_cost, gap_cost, false);
  ptr->nthreads = nthreads;
  ptr->self_pairs =
    run_anchored_lcp_search(ptr->data, ptr->data, ptr->params, nthreads, true);
  return AnchoredStarTreeRXPtr(ptr.release(), true);
}

// [[Rcpp::export(rng = false)]]
DataFrame AnchoredStarTree_self_search(CharacterVector sequences,
                                       const int max_distance,
                                       const int mismatch_cost = 1,
                                       const int gap_cost = 1,
                                       const int nthreads = 1,
                                       const bool show_progress = false) {
  if(show_progress) {
    Rcpp::warning("show_progress is not currently implemented for StarTree anchored mode");
  }
  if(nthreads < 1) {
    throw std::runtime_error("nthreads must be a single positive integer");
  }

  auto raw_sequences = character_vector_to_strings(sequences);
  validate_star_sequences(raw_sequences, "sequences");
  validate_star_params(max_distance, mismatch_cost, gap_cost);
  const startree::AnchoredInputData data =
    startree::make_anchored_input_data(std::move(raw_sequences), true);
  const startree::SearchParams params =
    startree::make_search_params(max_distance, mismatch_cost, gap_cost, false);
  const auto pairs = run_anchored_lcp_search(data, data, params, nthreads, true);
  return startree_pairs_to_dataframe(pairs, data, data);
}

// [[Rcpp::export(rng = false)]]
double AnchoredStarTree_size(AnchoredStarTreeRXPtr xp) {
  return static_cast<double>(xp->data.seqs.size());
}

// [[Rcpp::export(rng = false)]]
CharacterVector AnchoredStarTree_to_vector(AnchoredStarTreeRXPtr xp) {
  CharacterVector out(xp->data.seqs.size());
  for(size_t i = 0; i < xp->data.seqs.size(); ++i) {
    SET_STRING_ELT(out, checked_r_xlen(i),
                   string_to_charsxp(xp->data.seqs[i].seq));
  }
  return out;
}

// [[Rcpp::export(rng = false)]]
DataFrame AnchoredStarTree_result(AnchoredStarTreeRXPtr xp) {
  return startree_pairs_to_dataframe(xp->self_pairs, xp->data, xp->data);
}

// [[Rcpp::export(rng = false)]]
DataFrame AnchoredStarTree_search(AnchoredStarTreeRXPtr xp,
                                  CharacterVector query,
                                  const int nthreads = 1,
                                  const bool show_progress = false) {
  if(show_progress) {
    Rcpp::warning("show_progress is not currently implemented for StarTree anchored mode");
  }
  if(nthreads < 1) {
    throw std::runtime_error("nthreads must be a single positive integer");
  }
  if(Rf_xlength(query) == 0) {
    return startree_empty_match_result();
  }

  auto raw_query = character_vector_to_strings(query);
  validate_star_sequences(raw_query, "query");
  const startree::AnchoredInputData query_data =
    startree::make_anchored_input_data(std::move(raw_query), false, true);
  const auto pairs =
    run_anchored_lcp_search(xp->data, query_data, xp->params, nthreads, false);
  return startree_pairs_to_dataframe(pairs, xp->data, query_data);
}
