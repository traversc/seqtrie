#ifndef STARTREE_HAMMING_TRIE_H
#define STARTREE_HAMMING_TRIE_H

// Hamming-distance specialization of the StarTree global self-join.
//
// Semantics (matching the StarTree DNA convention, see common.h):
//   - only equal-length sequences can match (substitutions only, no indels);
//   - A/C/G/T compare by identity; N ALWAYS mismatches (target N -> code 0,
//     query N -> kNoMatch, which never equals any target code).
//
// Reuses the global block self-join structure (radix::process_target_block_self):
// length-sorted blocks, a per-block radix trie built incrementally, deferred
// set_sequence for the within-block lower triangle, length blocking across
// blocks, and LCP (start/seed) pebble resume across sorted queries. The banded
// edit-distance DP collapses to a scalar mismatch count, so the per-node int8
// caches and the band disappear. (A k-mer prefilter was evaluated and removed:
// on immune-repertoire data the conserved germline framework makes it skip ~0%,
// and the trie already amortizes the conserved prefix via LCP resume.)
//
// The trie is built over the height-padded codes (leading kPad), so every
// root->leaf path has length == height and sequences are terminal only at
// leaves. Padding makes unequal-length pairs accrue >= |dLen| pad/base
// mismatches; the equal-length requirement is enforced exactly by the
// seq_len == query_real_len gate at emission.

#include "simple_array/small_array.h"
#include "startree/bits.h"
#include "startree/common.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace startree {

constexpr uint32_t kNoHammingChild = std::numeric_limits<uint32_t>::max();

template <typename T>
using HammingSmallArr =
  trqwe::small_array<T, std::allocator<T>, uint16_t,
                     std::integral_constant<uint16_t, 8>>;

class HammingTrie {
  struct Node {
    HammingSmallArr<uint8_t> edge;
    std::array<uint32_t, kAlphabet> child;
    uint32_t sequence_index = kNoSequence;
    uint32_t seq_len = 0;  // real (un-padded) length; valid at terminal leaves
    uint8_t child_mask = 0;

    Node() { child.fill(kNoHammingChild); }
  };

  // LCP-resume frontier entry: the trie position (node_idx, last consumed edge
  // index) and the scalar mismatch count of the shared query prefix at that
  // depth. Far simpler than the edit-distance pebble (no DP column to carry).
  struct Pebble {
    uint32_t node_idx;
    int32_t char_offset;
    int distance;
  };

 public:
  void init(const int height) {
    height_ = height;
    nodes_.clear();
    nodes_.emplace_back();
    pebbles_.assign(static_cast<size_t>(height + 2), {});
    pebbles_[0].push_back(Pebble{0, -1, 0});
  }

  void reserve_nodes(const size_t n) {
    nodes_.reserve(std::max<size_t>(n, 1));
  }

  // Insert a height-padded code string (codes 0..kPad). Returns the leaf index;
  // call set_sequence() afterwards to make it emit (deferred for the within-block
  // lower triangle).
  uint32_t insert_path(const std::string& padded) {
    uint32_t node_idx = 0;
    int pos = 0;
    while(pos < height_) {
      const int c = static_cast<unsigned char>(padded[static_cast<size_t>(pos)]);
      const uint32_t child_idx = nodes_[node_idx].child[static_cast<size_t>(c)];
      if(child_idx == kNoHammingChild) {
        const uint32_t new_idx = make_node();
        Node& nn = nodes_[new_idx];
        const uint16_t len = static_cast<uint16_t>(height_ - pos);
        nn.edge = HammingSmallArr<uint8_t>(
          reinterpret_cast<const uint8_t*>(padded.data() + pos), len);
        nodes_[node_idx].child[static_cast<size_t>(c)] = new_idx;
        nodes_[node_idx].child_mask |= static_cast<uint8_t>(1U << c);
        nodes_[new_idx].seq_len = padded_real_len(padded);
        return new_idx;
      }

      const Node& child = nodes_[child_idx];
      const int edge_len = static_cast<int>(child.edge.size());
      int matched = 0;
      while(matched < edge_len && pos + matched < height_ &&
            child.edge[static_cast<size_t>(matched)] ==
              static_cast<uint8_t>(padded[static_cast<size_t>(pos + matched)])) {
        ++matched;
      }
      if(matched == edge_len) {
        pos += matched;
        node_idx = child_idx;
        continue;
      }

      const uint32_t tail_idx = make_node();
      {
        Node& tail = nodes_[tail_idx];
        Node& ch = nodes_[child_idx];
        const uint16_t suffix_len = static_cast<uint16_t>(edge_len - matched);
        tail.edge = HammingSmallArr<uint8_t>(ch.edge.data() + matched, suffix_len);
        tail.child = ch.child;
        tail.child_mask = ch.child_mask;
        tail.sequence_index = ch.sequence_index;
        tail.seq_len = ch.seq_len;
      }
      {
        Node& ch = nodes_[child_idx];
        ch.edge.resize(static_cast<uint16_t>(matched));
        ch.child.fill(kNoHammingChild);
        ch.child_mask = 0;
        ch.sequence_index = kNoSequence;
        ch.seq_len = 0;
      }
      {
        Node& tail = nodes_[tail_idx];
        const int tc = tail.edge[0];
        nodes_[child_idx].child[static_cast<size_t>(tc)] = tail_idx;
        nodes_[child_idx].child_mask |= static_cast<uint8_t>(1U << tc);
      }

      pos += matched;
      if(pos == height_) {
        nodes_[child_idx].seq_len = padded_real_len(padded);
        return child_idx;
      }

      const uint32_t branch_idx = make_node();
      Node& branch = nodes_[branch_idx];
      const uint16_t blen = static_cast<uint16_t>(height_ - pos);
      branch.edge = HammingSmallArr<uint8_t>(
        reinterpret_cast<const uint8_t*>(padded.data() + pos), blen);
      branch.seq_len = padded_real_len(padded);
      const int bc = branch.edge[0];
      nodes_[child_idx].child[static_cast<size_t>(bc)] = branch_idx;
      nodes_[child_idx].child_mask |= static_cast<uint8_t>(1U << bc);
      return branch_idx;
    }
    return node_idx;
  }

  void set_sequence(const uint32_t leaf_idx, const uint32_t seq_idx) {
    nodes_[static_cast<size_t>(leaf_idx)].sequence_index = seq_idx;
  }

  // start_depth: resume the DFS from the frontier seeded at this depth (== LCP
  //   with the previous query). seed_depth: seed the frontier down to this depth
  //   for the NEXT query. Both 0 disables LCP reuse (fresh descent from root).
  void search(const std::string& padded_query, const SearchParams& params,
              int start_depth, int seed_depth, std::vector<Hit>* hits) {
    threshold_ = params.max_distance;
    trans_.resize(static_cast<size_t>(height_));
    static constexpr int kInsertToQuery[kAlphabet] = {kNoMatch, 1, 2, 3, 4, kPad};
    int pad = 0;
    bool counting_pad = true;
    for(int i = 0; i < height_; ++i) {
      const int code = kInsertToQuery[static_cast<unsigned char>(
        padded_query[static_cast<size_t>(i)])];
      trans_[static_cast<size_t>(i)] = code;
      if(counting_pad && code == kPad) {
        ++pad;
      } else {
        counting_pad = false;
      }
    }
    query_real_len_ = static_cast<uint32_t>(height_ - pad);

    start_depth = std::max(start_depth, 0);
    start_depth = std::min(start_depth, height_ - 1);
    seed_depth_ = std::min(std::max(seed_depth, 0), height_);
    for(int d = start_depth + 1; d <= seed_depth_; ++d) {
      pebbles_[static_cast<size_t>(d)].clear();
    }

    hits_ = hits;
    for(const Pebble& peb : pebbles_[static_cast<size_t>(start_depth)]) {
      descend(peb.node_idx, peb.char_offset, start_depth, peb.distance);
    }
    hits_ = nullptr;
  }

  size_t node_count() const noexcept { return nodes_.size(); }

 private:
  uint32_t make_node() {
    nodes_.emplace_back();
    return static_cast<uint32_t>(nodes_.size() - 1);
  }

  uint32_t padded_real_len(const std::string& padded) const noexcept {
    int pad = 0;
    while(pad < height_ &&
          static_cast<unsigned char>(padded[static_cast<size_t>(pad)]) == kPad) {
      ++pad;
    }
    return static_cast<uint32_t>(height_ - pad);
  }

  void emit(const Node& node, const int distance) {
    if(node.sequence_index != kNoSequence && distance <= threshold_ &&
       node.seq_len == query_real_len_) {
      hits_->push_back(Hit{node.sequence_index, distance});
    }
  }

  // Continue from (node_idx, char_offset) at target_depth == depth with running
  // mismatch count dist; recurse into children at the node end.
  void descend(const uint32_t node_idx, const int char_offset, int depth,
               int dist) {
    const Node& node = nodes_[static_cast<size_t>(node_idx)];
    const int edge_len = static_cast<int>(node.edge.size());
    for(int i = char_offset + 1; i < edge_len; ++i) {
      const int qc = trans_[static_cast<size_t>(depth)];
      dist += (static_cast<int>(node.edge[static_cast<size_t>(i)]) == qc) ? 0 : 1;
      ++depth;
      if(dist > threshold_) {
        return;
      }
      if(depth <= seed_depth_) {
        pebbles_[static_cast<size_t>(depth)].push_back(Pebble{node_idx, i, dist});
      }
    }
    emit(node, dist);
    unsigned mask = node.child_mask;
    while(mask) {
      const int code = ctz_u(mask);
      mask &= mask - 1;
      descend(node.child[static_cast<size_t>(code)], -1, depth, dist);
    }
  }

  int height_ = 0;
  int threshold_ = 0;
  int seed_depth_ = 0;
  uint32_t query_real_len_ = 0;
  std::vector<Node> nodes_;
  std::vector<std::vector<Pebble>> pebbles_;
  std::vector<int> trans_;
  std::vector<Hit>* hits_ = nullptr;
};

namespace hamming {

// Hamming length blocking: only equal-length sequences can match, so two length
// ranges interact iff they overlap, and a single query length is relevant to a
// block iff it lies within the block's [min_len, max_len].
inline bool blocks_can_match(const Block& a, const Block& b) noexcept {
  return a.min_len <= b.max_len && b.min_len <= a.max_len;
}
inline bool seq_len_can_match(const Block& block, const size_t qlen) noexcept {
  return block.min_len <= qlen && qlen <= block.max_len;
}

// Self-join over one target block (Hamming analogue of
// radix::process_target_block_self). Emits each unordered pair once via the
// deferred-set_sequence lower triangle within the block, and once per block
// pair across blocks. lcp toggles the LCP pebble resume (kept for benchmarking;
// production uses the default true).
inline void process_target_block_self(const std::vector<Sequence>& seqs,
                                       const std::vector<Block>& blocks,
                                       const size_t target_block_id,
                                       const int height,
                                       const SearchParams& params,
                                       std::vector<PairRecord>* out,
                                       const bool lcp = true) {
  const Block& target_block = blocks[target_block_id];

  HammingTrie trie;
  trie.init(height);
  trie.reserve_nodes(1 + 2 * static_cast<size_t>(target_block.end - target_block.begin));

  std::vector<Hit> hits;
  hits.reserve(128);

  const auto append = [&](const Hit& hit, const uint32_t qi) {
    if(hit.distance > 0 && hit.distance <= params.max_distance) {
      out->push_back(PairRecord{hit.sequence_index, qi, hit.distance});
    }
  };

  const Sequence* last_query = nullptr;
  for(uint32_t qi = target_block.begin; qi < target_block.end; ++qi) {
    const Sequence& query = seqs[qi];
    const uint32_t leaf_node_idx = trie.insert_path(query.padded);

    int trail = 0;
    int start = 0;
    if(lcp) {
      if(qi + 1 < target_block.end) {
        trail = common_prefix_length(query.padded, seqs[qi + 1].padded);
      }
      if(last_query != nullptr) {
        start = common_prefix_length(query.padded, last_query->padded);
      }
    }

    hits.clear();
    trie.search(query.padded, params, start, trail, &hits);
    for(const Hit& hit : hits) {
      append(hit, qi);
    }
    last_query = &query;
    trie.set_sequence(leaf_node_idx, qi);  // deferred: qi never matches itself
  }

  const size_t query_jobs = (blocks.size() + 1) / 2;
  for(size_t step = 1; step < query_jobs; ++step) {
    const size_t query_block_id = (target_block_id + step) % blocks.size();
    const Block& query_block = blocks[query_block_id];
    if(!blocks_can_match(target_block, query_block)) {
      continue;
    }

    const Sequence* last_xq = nullptr;
    for(uint32_t qi = query_block.begin; qi < query_block.end; ++qi) {
      const Sequence& query = seqs[qi];
      if(!seq_len_can_match(target_block, query.seq.size())) {
        continue;
      }

      int trail = 0;
      int start = 0;
      if(lcp) {
        if(qi + 1 < query_block.end) {
          trail = common_prefix_length(query.padded, seqs[qi + 1].padded);
        }
        if(last_xq != nullptr) {
          start = common_prefix_length(query.padded, last_xq->padded);
        }
      }

      hits.clear();
      trie.search(query.padded, params, start, trail, &hits);
      for(const Hit& hit : hits) {
        append(hit, qi);
      }
      last_xq = &query;
    }
  }
}

// Query-vs-database over one target block (Hamming analogue of
// radix::process_target_block_query). The trie is built once from the target
// block, then every length-compatible query is searched against it.
inline void process_target_block_query(const std::vector<Sequence>& target_seqs,
                                        const std::vector<Sequence>& query_seqs,
                                        const std::vector<Block>& query_blocks,
                                        const Block& target_block,
                                        const int height,
                                        const SearchParams& params,
                                        std::vector<PairRecord>* out,
                                        const bool lcp = true) {
  HammingTrie trie;
  trie.init(height);
  trie.reserve_nodes(1 + 2 * static_cast<size_t>(target_block.end - target_block.begin));
  for(uint32_t ti = target_block.begin; ti < target_block.end; ++ti) {
    const uint32_t leaf = trie.insert_path(target_seqs[ti].padded);
    trie.set_sequence(leaf, ti);
  }

  std::vector<Hit> hits;
  hits.reserve(128);

  const auto append = [&](const Hit& hit, const uint32_t qi) {
    if(hit.distance <= params.max_distance &&
       (params.include_zero || hit.distance > 0)) {
      out->push_back(PairRecord{hit.sequence_index, qi, hit.distance});
    }
  };

  for(const Block& query_block : query_blocks) {
    if(!blocks_can_match(target_block, query_block)) {
      continue;
    }
    const Sequence* last_query = nullptr;
    for(uint32_t qi = query_block.begin; qi < query_block.end; ++qi) {
      const Sequence& query = query_seqs[qi];
      if(!seq_len_can_match(target_block, query.seq.size())) {
        continue;
      }

      int trail = 0;
      int start = 0;
      if(lcp) {
        if(qi + 1 < query_block.end) {
          trail = common_prefix_length(query.padded, query_seqs[qi + 1].padded);
        }
        if(last_query != nullptr) {
          start = common_prefix_length(query.padded, last_query->padded);
        }
      }

      hits.clear();
      trie.search(query.padded, params, start, trail, &hits);
      for(const Hit& hit : hits) {
        append(hit, qi);
      }
      last_query = &query;
    }
  }
}

}  // namespace hamming
}  // namespace startree

#endif  // STARTREE_HAMMING_TRIE_H
