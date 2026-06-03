#ifndef STARTREE_STARTREE_H
#define STARTREE_STARTREE_H

// StarTree is a modified version of the Starcode all-pairs search strategy
// described by Eduard Zorita, Pol Cusco, and Guillaume J. Filion (2015), doi:10.1093/bioinformatics/btv053, adapted to operate over a radix trie. The algorithm is credited to the Starcode authors; this package provides a separate, modified implementation.

#include "startree/common.h"
#include "startree/anchored_custom_trie.h"
#include "startree/anchored_standard_trie.h"
#include "startree/custom_trie.h"
#include "startree/hamming_trie.h"
#include "startree/prefilter.h"
#include "startree/standard_trie.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace startree {

inline void append_pair(const uint32_t target_index,
                        const uint32_t query_index,
                        const int distance,
                        std::vector<PairRecord>* out) {
  PairRecord rec;
  rec.target_index = target_index;
  rec.query_index = query_index;
  rec.distance = distance;
  out->push_back(std::move(rec));
}

inline void append_hit(const uint32_t target_index,
                       const uint32_t query_index,
                       const int distance,
                       const SearchParams& params,
                       std::vector<PairRecord>* out) {
  if(distance <= params.max_distance && (params.include_zero || distance > 0)) {
    append_pair(target_index, query_index, distance, out);
  }
}

inline void append_self_hit(const uint32_t target_index,
                            const uint32_t query_index,
                            const int distance,
                            const SearchParams& params,
                            std::vector<PairRecord>* out) {
  if(distance > 0 && distance <= params.max_distance) {
    append_pair(target_index, query_index, distance, out);
  }
}

inline size_t count_radix_nodes(const Block& block) {
  if(block.begin == block.end) {
    return 1;
  }
  return 1 + 2 * static_cast<size_t>(block.end - block.begin);
}

namespace radix {

template <int MaxBand, typename TrieType>
void process_target_block_self(const std::vector<Sequence>& seqs,
                               const std::vector<Block>& blocks,
                               const size_t target_block_id,
                               const int height,
                               const int median_len,
                               const SearchParams& params,
                               std::vector<PairRecord>* out) {
  const Block& target_block = blocks[target_block_id];
  const int target_pad_depth_max =
    height - static_cast<int>(target_block.min_len);

  TrieType trie;
  trie.init(height, target_pad_depth_max);
  trie.reserve_nodes(count_radix_nodes(target_block));
  UnitCostPrefilter unit_prefilter;
  WeightedCostPrefilter weighted_prefilter;
  if(params.custom_cost) {
    weighted_prefilter.reset(median_len, height, params.lookup_distance,
                             params.max_distance,
                             std::min(params.mismatch_cost, params.gap_cost),
                             params.gap_cost);
  } else {
    unit_prefilter.reset(median_len, height, params.lookup_distance);
  }

  std::vector<Hit> hits;
  hits.reserve(128);

  const size_t query_jobs = (blocks.size() + 1) / 2;
  const Sequence* last_query = nullptr;

  for(uint32_t qi = target_block.begin; qi < target_block.end; ++qi) {
    const Sequence& query = seqs[qi];
    if(params.custom_cost) {
      weighted_prefilter.insert(query.padded);
    } else {
      unit_prefilter.insert(query.padded);
    }
    const uint32_t leaf_node_idx = trie.insert_path(query.padded);

    int trail = 0;
    if(qi + 1 < target_block.end) {
      trail = common_prefix_length(query.padded, seqs[qi + 1].padded);
    }
    int start = 0;
    if(last_query != nullptr) {
      start = common_prefix_length(query.padded, last_query->padded);
    }

    hits.clear();
    trie.search(query.padded, params, start, trail, &hits);
    for(const Hit& hit : hits) {
      append_self_hit(hit.sequence_index, qi, hit.distance, params, out);
    }
    last_query = &query;
    trie.set_sequence(leaf_node_idx, qi);
  }

  for(size_t step = 1; step < query_jobs; ++step) {
    const size_t query_block_id = (target_block_id + step) % blocks.size();
    const Block& query_block = blocks[query_block_id];
    if(!length_ranges_can_match(target_block, query_block, params.max_distance,
                                params.gap_cost)) {
      continue;
    }

    const Sequence* last_xq = nullptr;
    for(uint32_t qi = query_block.begin; qi < query_block.end; ++qi) {
      const Sequence& query = seqs[qi];
      if(!sequence_length_can_match(target_block, query.seq.size(),
                                    params.max_distance, params.gap_cost)) {
        continue;
      }
      if(params.custom_cost) {
        if(!weighted_prefilter.maybe_contains(query.padded)) {
          continue;
        }
      } else if(!unit_prefilter.maybe_contains(query.padded)) {
        continue;
      }

      int trail = 0;
      if(qi + 1 < query_block.end) {
        trail = common_prefix_length(query.padded, seqs[qi + 1].padded);
      }
      int start = 0;
      if(last_xq != nullptr) {
        start = common_prefix_length(query.padded, last_xq->padded);
      }

      hits.clear();
      trie.search(query.padded, params, start, trail, &hits);
      for(const Hit& hit : hits) {
        append_self_hit(hit.sequence_index, qi, hit.distance, params, out);
      }
      last_xq = &query;
    }
  }
}

template <int MaxBand, typename TrieType>
void process_target_block_query(const std::vector<Sequence>& target_seqs,
                                const std::vector<Sequence>& query_seqs,
                                const std::vector<Block>& query_blocks,
                                const Block& target_block,
                                const int height,
                                const int median_len,
                                const SearchParams& params,
                                std::vector<PairRecord>* out) {
  TrieType trie;
  const int target_pad_depth_max =
    height - static_cast<int>(target_block.min_len);

  trie.init(height, target_pad_depth_max);
  trie.reserve_nodes(count_radix_nodes(target_block));
  UnitCostPrefilter unit_prefilter;
  WeightedCostPrefilter weighted_prefilter;
  if(params.custom_cost) {
    weighted_prefilter.reset(median_len, height, params.lookup_distance,
                             params.max_distance,
                             std::min(params.mismatch_cost, params.gap_cost),
                             params.gap_cost);
  } else {
    unit_prefilter.reset(median_len, height, params.lookup_distance);
  }

  for(uint32_t ti = target_block.begin; ti < target_block.end; ++ti) {
    const Sequence& target = target_seqs[ti];
    if(params.custom_cost) {
      weighted_prefilter.insert(target.padded);
    } else {
      unit_prefilter.insert(target.padded);
    }
    const uint32_t leaf_node_idx = trie.insert_path(target.padded);
    trie.set_sequence(leaf_node_idx, ti);
  }

  std::vector<Hit> hits;
  hits.reserve(128);

  for(const Block& query_block : query_blocks) {
    if(!length_ranges_can_match(target_block, query_block, params.max_distance,
                                params.gap_cost)) {
      continue;
    }

    const Sequence* last_query = nullptr;
    for(uint32_t qi = query_block.begin; qi < query_block.end; ++qi) {
      const Sequence& query = query_seqs[qi];
      if(!sequence_length_can_match(target_block, query.seq.size(),
                                    params.max_distance, params.gap_cost)) {
        continue;
      }
      if(params.custom_cost) {
        if(!weighted_prefilter.maybe_contains(query.padded)) {
          continue;
        }
      } else if(!unit_prefilter.maybe_contains(query.padded)) {
        continue;
      }

      int trail = 0;
      if(qi + 1 < query_block.end) {
        trail = common_prefix_length(query.padded, query_seqs[qi + 1].padded);
      }
      int start = 0;
      if(last_query != nullptr) {
        start = common_prefix_length(query.padded, last_query->padded);
      }

      hits.clear();
      trie.search(query.padded, params, start, trail, &hits);
      for(const Hit& hit : hits) {
        append_hit(hit.sequence_index, qi, hit.distance, params, out);
      }
      last_query = &query;
    }
  }
}

}  // namespace radix

template <int MaxBand>
inline void process_target_block_self_scored(const std::vector<Sequence>& seqs,
                                             const std::vector<Block>& blocks,
                                             const size_t target_block_id,
                                             const int height,
                                             const int median_len,
                                             const SearchParams& params,
                                             std::vector<PairRecord>* out) {
  if(params.custom_cost) {
    radix::process_target_block_self<MaxBand, radix::CustomTrie<MaxBand>>(
      seqs, blocks, target_block_id, height, median_len, params, out
    );
  } else {
    radix::process_target_block_self<MaxBand, radix::StandardTrie<MaxBand>>(
      seqs, blocks, target_block_id, height, median_len, params, out
    );
  }
}

inline void process_target_block_self(const std::vector<Sequence>& seqs,
                                      const std::vector<Block>& blocks,
                                      const size_t target_block_id,
                                      const int height,
                                      const int median_len,
                                      const SearchParams& params,
                                      std::vector<PairRecord>* out) {
  switch(std::max(1, params.band_radius)) {
    case 1:
      process_target_block_self_scored<1>(
        seqs, blocks, target_block_id, height, median_len, params, out
      );
      return;
    case 2:
      process_target_block_self_scored<2>(
        seqs, blocks, target_block_id, height, median_len, params, out
      );
      return;
    case 3:
      process_target_block_self_scored<3>(
        seqs, blocks, target_block_id, height, median_len, params, out
      );
      return;
    case 4:
      process_target_block_self_scored<4>(
        seqs, blocks, target_block_id, height, median_len, params, out
      );
      return;
    case 5:
      process_target_block_self_scored<5>(
        seqs, blocks, target_block_id, height, median_len, params, out
      );
      return;
    case 6:
      process_target_block_self_scored<6>(
        seqs, blocks, target_block_id, height, median_len, params, out
      );
      return;
    case 7:
      process_target_block_self_scored<7>(
        seqs, blocks, target_block_id, height, median_len, params, out
      );
      return;
    case 8:
      process_target_block_self_scored<8>(
        seqs, blocks, target_block_id, height, median_len, params, out
      );
      return;
  }
}

template <int MaxBand>
inline void process_target_block_query_scored(const std::vector<Sequence>& target_seqs,
                                              const std::vector<Sequence>& query_seqs,
                                              const std::vector<Block>& query_blocks,
                                              const Block& target_block,
                                              const int height,
                                              const int median_len,
                                              const SearchParams& params,
                                              std::vector<PairRecord>* out) {
  if(params.custom_cost) {
    radix::process_target_block_query<MaxBand, radix::CustomTrie<MaxBand>>(
      target_seqs, query_seqs, query_blocks, target_block, height, median_len,
      params, out
    );
  } else {
    radix::process_target_block_query<MaxBand, radix::StandardTrie<MaxBand>>(
      target_seqs, query_seqs, query_blocks, target_block, height, median_len,
      params, out
    );
  }
}

inline void process_target_block_query(const std::vector<Sequence>& target_seqs,
                                       const std::vector<Sequence>& query_seqs,
                                       const std::vector<Block>& query_blocks,
                                       const Block& target_block,
                                       const int height,
                                       const int median_len,
                                       const SearchParams& params,
                                       std::vector<PairRecord>* out) {
  switch(std::max(1, params.band_radius)) {
    case 1:
      process_target_block_query_scored<1>(
        target_seqs, query_seqs, query_blocks, target_block, height, median_len,
        params, out
      );
      return;
    case 2:
      process_target_block_query_scored<2>(
        target_seqs, query_seqs, query_blocks, target_block, height, median_len,
        params, out
      );
      return;
    case 3:
      process_target_block_query_scored<3>(
        target_seqs, query_seqs, query_blocks, target_block, height, median_len,
        params, out
      );
      return;
    case 4:
      process_target_block_query_scored<4>(
        target_seqs, query_seqs, query_blocks, target_block, height, median_len,
        params, out
      );
      return;
    case 5:
      process_target_block_query_scored<5>(
        target_seqs, query_seqs, query_blocks, target_block, height, median_len,
        params, out
      );
      return;
    case 6:
      process_target_block_query_scored<6>(
        target_seqs, query_seqs, query_blocks, target_block, height, median_len,
        params, out
      );
      return;
    case 7:
      process_target_block_query_scored<7>(
        target_seqs, query_seqs, query_blocks, target_block, height, median_len,
        params, out
      );
      return;
    case 8:
      process_target_block_query_scored<8>(
        target_seqs, query_seqs, query_blocks, target_block, height, median_len,
        params, out
      );
      return;
  }
}

namespace anchored {

template <typename TrieType>
inline TrieType build_trie(const std::vector<std::string>& target_codes) {
  TrieType trie;
  trie.reserve_nodes(1 + 2 * target_codes.size());
  for(size_t i = 0; i < target_codes.size(); ++i) {
    trie.insert_path(target_codes[i], static_cast<uint32_t>(i));
  }
  trie.compact_search_order();
  return trie;
}

template <typename Fn>
decltype(auto) dispatch_band(const SearchParams& params, Fn&& fn) {
  switch(std::max(1, params.band_radius)) {
    case 1: return fn(std::integral_constant<int, 1>{});
    case 2: return fn(std::integral_constant<int, 2>{});
    case 3: return fn(std::integral_constant<int, 3>{});
    case 4: return fn(std::integral_constant<int, 4>{});
    case 5: return fn(std::integral_constant<int, 5>{});
    case 6: return fn(std::integral_constant<int, 6>{});
    case 7: return fn(std::integral_constant<int, 7>{});
    default: return fn(std::integral_constant<int, kMaxTau>{});
  }
}

}  // namespace anchored

}  // namespace startree

#endif  // STARTREE_STARTREE_H
