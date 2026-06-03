#ifndef STARTREE_ANCHORED_CUSTOM_TRIE_H
#define STARTREE_ANCHORED_CUSTOM_TRIE_H

// Anchored (semi-global) search with the custom_trie.h weighted-cost model:
//  - tight banded DP over integer caches
//  - copy-free LCP pebble-frontier resume across sorted queries
//  - dash() exact-tail shortcut once the edit budget is exhausted
//
// DP orientation (no padding; both sequences anchored at position 0):
//   target_depth d = number of target chars consumed.
//   cache offset a in [-radius, radius]; cell(a) = editdist(query[0:d+a], target[0:d]).
//   col_min = min over band = best alignment of full target-so-far to a query prefix
//             (Ukkonen: non-decreasing in d, so col_min > tau prunes the subtree).
//   row candidate at depth d = cell(query_len - d) = editdist(full query, target[0:d])
//             (case a: query fully aligned to a target prefix). row_min tracks its min.
// Emission:
//   case (b): a target node (sequence_index set) with col_min <= tau emits at col_min.
//   case (a): when row_min <= tau and row_min <= col_min, every sequence in the subtree
//             matches at row_min -> add_all, stop descending.

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
#include <type_traits>
#include <vector>

#ifndef STARTREE_ALWAYS_INLINE
#  if defined(__GNUC__) || defined(__clang__)
#    define STARTREE_ALWAYS_INLINE inline __attribute__((always_inline))
#  elif defined(_MSC_VER)
#    define STARTREE_ALWAYS_INLINE __forceinline
#  else
#    define STARTREE_ALWAYS_INLINE inline
#  endif
#endif

namespace startree {
namespace radix {

#ifdef ANCHORED_CUSTOM_TRIE_COUNT
inline thread_local size_t g_anchored_custom_columns = 0;     // update_cache calls (non-empty)
inline thread_local size_t g_anchored_custom_cells = 0;       // band cells computed
inline thread_local size_t g_anchored_custom_live_cells = 0;  // cells with value <= threshold
#endif

constexpr uint32_t kNoAnchoredCustomChild = std::numeric_limits<uint32_t>::max();
constexpr int kAnchoredCustomInf = std::numeric_limits<int>::max() / 4;

// One sentinel slot beyond the band on each side keeps band-edge reads/writes
// (offset +/- (radius+1)) in bounds even when radius == MaxBand. Center index is
// kAnchoredCustomCenter<MaxBand> = MaxBand + 1.
template <int MaxBand>
using AnchoredCustomCache = std::array<int, 2 * MaxBand + 3>;

template <int MaxBand>
constexpr int kAnchoredCustomCenter = MaxBand + 1;

template <typename T>
using AnchoredCustomSmallArr =
  trqwe::small_array<T, std::allocator<T>, uint16_t,
                     std::integral_constant<uint16_t, 8>>;

template <int MaxBand>
struct AnchoredCustomNode {
  using Cache = AnchoredCustomCache<MaxBand>;

  AnchoredCustomSmallArr<uint8_t> edge;
  std::array<uint32_t, kAlphabet> child;
  uint32_t sequence_index = kNoSequence;
  // Min sequence_index over this node's subtree (self-join lower-triangle prune).
  // Valid only after compact_search_order(); kNoSequence for an empty subtree.
  uint32_t subtree_min_seq = kNoSequence;
  uint8_t child_mask = 0;

  AnchoredCustomNode() { child.fill(kNoAnchoredCustomChild); }
};

// Frontier entry for LCP resume, recorded at target_depth == seed_depth.
// edge_offset = index of the last consumed edge char of node_idx; the snapshot is
// the DP column right after consuming it. To resume: if edge_offset+1 < edge_len,
// continue the edge from edge_offset+1; otherwise iterate the node's children.
// At resume depth (<= query_len - radius - 1) the running row_min is always inf,
// so it is not stored.
template <int MaxBand>
struct AnchoredCustomPebble {
  uint32_t node_idx;
  int32_t edge_offset;
  AnchoredCustomCache<MaxBand> cache;
};

// Per-thread scratch for LCP-reuse search. Kept outside the trie so a single
// const trie can be shared across worker threads (only this workspace is
// duplicated per thread). pebbles[d] holds the DFS frontier seeded at
// target_depth d; hit_log[d] replays shallow case-b hits on resume.
template <int MaxBand>
struct AnchoredCustomWorkspace {
  using Pebble = AnchoredCustomPebble<MaxBand>;
  std::vector<std::vector<Pebble>> pebbles;
  std::vector<std::vector<std::pair<uint32_t, int>>> hit_log;  // (seq_idx, distance)
};

template <int MaxBand>
class AnchoredCustomTrie {
  using NodeType = AnchoredCustomNode<MaxBand>;
  using Cache = AnchoredCustomCache<MaxBand>;
  using Workspace = AnchoredCustomWorkspace<MaxBand>;

  struct SearchState {
    const unsigned char* query = nullptr;
    int query_len = 0;
    int threshold = 0;
    int radius = 0;
    int mismatch_cost = 1;
    int gap_cost = 1;
    int min_edit_cost = 1;
    int inf = kAnchoredCustomInf;
    uint32_t query_index = 0;
    std::vector<Hit>* hits = nullptr;
    std::vector<PairRecord>* pairs = nullptr;
    Workspace* ws = nullptr;  // only used when seeding/resuming (LCP path)
    int start_depth = 0;
    int seed_depth = 0;
    // Self-join lower-triangle: emit only targets with sequence_index < cutoff,
    // and skip subtrees whose min sequence_index >= cutoff. kNoSequence disables
    // it (keep everything) for plain query-vs-database search.
    uint32_t triangle_cutoff = kNoSequence;
  };

 public:
  AnchoredCustomTrie() { nodes_.emplace_back(); }

  void clear();
  void reserve_nodes(size_t n);
  uint32_t insert_path(const std::string& target_codes, uint32_t seq_idx);
  // Reorder nodes into DFS pre-order so a parent's subtree is contiguous in
  // memory (prefetch-friendly traversal). Call once after all inserts.
  void compact_search_order();

  // Single-query entry points (no LCP reuse).
  void search(const std::string& query_codes, const SearchParams& params,
              std::vector<Hit>* hits) const;
  void search_pairs(const std::string& query_codes, const SearchParams& params,
                    uint32_t query_index, std::vector<PairRecord>* pairs) const;

  // Batched LCP-reuse entry point: query_codes[begin:end) must be sorted.
  // The workspace is per-thread scratch; the trie itself stays const/shared.
  // When lower_triangle is true the query set must equal the (lex-sorted) target
  // set with query position qi == target sequence_index; only pairs (j, qi) with
  // j < qi are emitted (each unordered self pair once, diagonal excluded).
  void search_batch_lcp(const std::vector<std::string>& query_codes,
                        const SearchParams& params, size_t begin, size_t end,
                        std::vector<PairRecord>* pairs, Workspace* ws,
                        bool lower_triangle = false) const;

  size_t node_count() const noexcept { return nodes_.size(); }

 private:
  uint32_t make_node();
  uint32_t compact_subtree(uint32_t old_idx, std::vector<NodeType>* out) const;
  static Cache inf_cache() noexcept;
  static int cell(const Cache& c, int offset) noexcept;
  static int substitution_cost(int t, int q, const SearchState& st) noexcept;
  static int add_cost(int value, int cost, const SearchState& st) noexcept;
  static int gap_prefix_cost(int len, const SearchState& st) noexcept;
  static void validate_params(const SearchParams& params);

  Cache make_root_cache(const SearchState& st) const noexcept;
  // Fast banded update: parent column -> child column for target char c at depth d.
  // Returns col_min. Writes child cache (only band cells + boundary sentinels).
  int update_cache(const Cache& parent, int target_code, int depth,
                   const SearchState& st, Cache* out) const noexcept;
  int row_candidate(const Cache& c, int depth, const SearchState& st) const noexcept;

  void emit_hit(uint32_t seq_idx, int distance, const SearchState& st) const;
  void add_all(uint32_t node_idx, int distance, const SearchState& st) const;

  void descend_edge(uint32_t node_idx, int edge_start, const Cache& cache_in,
                    int depth_in, int row_min_in, const SearchState& st) const;
  bool can_dash(const Cache& c, int depth, const SearchState& st) const noexcept;
  void dash(uint32_t node_idx, int edge_i, int depth, int distance,
            const SearchState& st) const;
  void iterate_children(uint32_t node_idx, const Cache& parent_cache,
                        int parent_depth, int parent_row_min,
                        const SearchState& st) const;
  void run_query(const SearchState& st, const Cache& root_cache,
                 int root_row_min) const;

  using Pebble = AnchoredCustomPebble<MaxBand>;
  std::vector<NodeType> nodes_;
};

template <int MaxBand>
void AnchoredCustomTrie<MaxBand>::clear() {
  nodes_.clear();
  nodes_.emplace_back();
}

template <int MaxBand>
void AnchoredCustomTrie<MaxBand>::reserve_nodes(const size_t n) {
  nodes_.reserve(std::max<size_t>(n, 1));
}

template <int MaxBand>
uint32_t AnchoredCustomTrie<MaxBand>::insert_path(
    const std::string& target_codes, const uint32_t seq_idx) {
  uint32_t node_idx = 0;
  int pos = 0;
  const int n = static_cast<int>(target_codes.size());

  if(n == 0) {
    nodes_[0].sequence_index = seq_idx;
    return 0;
  }

  while(pos < n) {
    const int c = static_cast<unsigned char>(target_codes[static_cast<size_t>(pos)]);
    if(c < 0 || c >= kAlphabet || c == kPad) {
      throw std::runtime_error("AnchoredCustomTrie target code is out of range");
    }
    const uint32_t child_idx = nodes_[node_idx].child[static_cast<size_t>(c)];
    if(child_idx == kNoAnchoredCustomChild) {
      const uint32_t new_idx = make_node();
      NodeType& nn = nodes_[new_idx];
      const uint16_t len = static_cast<uint16_t>(n - pos);
      nn.edge = AnchoredCustomSmallArr<uint8_t>(
        reinterpret_cast<const uint8_t*>(target_codes.data() + pos), len);
      nn.sequence_index = seq_idx;
      nodes_[node_idx].child[static_cast<size_t>(c)] = new_idx;
      nodes_[node_idx].child_mask |= static_cast<uint8_t>(1U << c);
      return new_idx;
    }

    const NodeType& child = nodes_[child_idx];
    const int edge_len = static_cast<int>(child.edge.size());
    int matched = 0;
    while(matched < edge_len && pos + matched < n &&
          child.edge[static_cast<size_t>(matched)] ==
            static_cast<uint8_t>(target_codes[static_cast<size_t>(pos + matched)])) {
      ++matched;
    }
    if(matched == edge_len) {
      pos += matched;
      node_idx = child_idx;
      continue;
    }

    const uint32_t tail_idx = make_node();
    {
      NodeType& tail = nodes_[tail_idx];
      NodeType& ch = nodes_[child_idx];
      const uint16_t suffix_len = static_cast<uint16_t>(edge_len - matched);
      tail.edge = AnchoredCustomSmallArr<uint8_t>(ch.edge.data() + matched, suffix_len);
      tail.child = ch.child;
      tail.child_mask = ch.child_mask;
      tail.sequence_index = ch.sequence_index;
    }
    {
      NodeType& ch = nodes_[child_idx];
      ch.edge.resize(static_cast<uint16_t>(matched));
      ch.child.fill(kNoAnchoredCustomChild);
      ch.child_mask = 0;
      ch.sequence_index = kNoSequence;
    }
    {
      NodeType& tail = nodes_[tail_idx];
      const int tc = tail.edge[0];
      nodes_[child_idx].child[static_cast<size_t>(tc)] = tail_idx;
      nodes_[child_idx].child_mask |= static_cast<uint8_t>(1U << tc);
    }

    pos += matched;
    if(pos == n) {
      nodes_[child_idx].sequence_index = seq_idx;
      return child_idx;
    }

    const uint32_t branch_idx = make_node();
    NodeType& branch = nodes_[branch_idx];
    const uint16_t blen = static_cast<uint16_t>(n - pos);
    branch.edge = AnchoredCustomSmallArr<uint8_t>(
      reinterpret_cast<const uint8_t*>(target_codes.data() + pos), blen);
    branch.sequence_index = seq_idx;
    const int bc = branch.edge[0];
    nodes_[child_idx].child[static_cast<size_t>(bc)] = branch_idx;
    nodes_[child_idx].child_mask |= static_cast<uint8_t>(1U << bc);
    return branch_idx;
  }

  nodes_[node_idx].sequence_index = seq_idx;
  return node_idx;
}

template <int MaxBand>
uint32_t AnchoredCustomTrie<MaxBand>::make_node() {
  nodes_.emplace_back();
  return static_cast<uint32_t>(nodes_.size() - 1);
}

template <int MaxBand>
uint32_t AnchoredCustomTrie<MaxBand>::compact_subtree(
    const uint32_t old_idx, std::vector<NodeType>* out) const {
  const uint32_t new_idx = static_cast<uint32_t>(out->size());
  out->push_back(nodes_[static_cast<size_t>(old_idx)]);
  // Re-link children to their new (DFS pre-order) indices. Recursion appends
  // each child's subtree immediately after this node.
  uint32_t min_seq = (*out)[new_idx].sequence_index;  // own; kNoSequence if none
  unsigned mask = (*out)[new_idx].child_mask;
  while(mask) {
    const int code = ctz_u(mask);
    mask &= mask - 1;
    const uint32_t old_child = nodes_[static_cast<size_t>(old_idx)]
                                 .child[static_cast<size_t>(code)];
    const uint32_t new_child = compact_subtree(old_child, out);
    (*out)[new_idx].child[static_cast<size_t>(code)] = new_child;
    const uint32_t cm = (*out)[new_child].subtree_min_seq;
    if(cm < min_seq) min_seq = cm;
  }
  (*out)[new_idx].subtree_min_seq = min_seq;
  return new_idx;
}

template <int MaxBand>
void AnchoredCustomTrie<MaxBand>::compact_search_order() {
  std::vector<NodeType> compacted;
  compacted.reserve(nodes_.size());
  compact_subtree(0, &compacted);
  nodes_ = std::move(compacted);
}

template <int MaxBand>
auto AnchoredCustomTrie<MaxBand>::inf_cache() noexcept -> Cache {
  Cache out{};
  out.fill(kAnchoredCustomInf);
  return out;
}

template <int MaxBand>
STARTREE_ALWAYS_INLINE
int AnchoredCustomTrie<MaxBand>::cell(const Cache& c, const int offset) noexcept {
  return c[static_cast<size_t>(offset + MaxBand + 1)];
}

template <int MaxBand>
STARTREE_ALWAYS_INLINE
int AnchoredCustomTrie<MaxBand>::substitution_cost(
    const int t, const int q, const SearchState& st) noexcept {
  return t == q ? 0 : st.mismatch_cost;
}

template <int MaxBand>
STARTREE_ALWAYS_INLINE
int AnchoredCustomTrie<MaxBand>::add_cost(
    const int value, const int cost, const SearchState& st) noexcept {
  if(value >= st.inf || cost >= st.inf - value) {
    return st.inf;
  }
  return value + cost;
}

template <int MaxBand>
STARTREE_ALWAYS_INLINE
int AnchoredCustomTrie<MaxBand>::gap_prefix_cost(
    const int len, const SearchState& st) noexcept {
  if(len <= 0) {
    return 0;
  }
  const int64_t cost =
    static_cast<int64_t>(len) * static_cast<int64_t>(st.gap_cost);
  return cost >= st.inf ? st.inf : static_cast<int>(cost);
}

template <int MaxBand>
void AnchoredCustomTrie<MaxBand>::validate_params(const SearchParams& params) {
  if(params.max_distance < 0) {
    throw std::runtime_error("AnchoredCustomTrie max_distance must be non-negative");
  }
  if(params.band_radius < 0 || params.band_radius > MaxBand) {
    throw std::runtime_error("AnchoredCustomTrie band radius exceeds MaxBand");
  }
  if(params.mismatch_cost <= 0 || params.gap_cost <= 0) {
    throw std::runtime_error("AnchoredCustomTrie costs must be positive");
  }
  if(params.max_distance >= kAnchoredCustomInf ||
     params.mismatch_cost >= kAnchoredCustomInf ||
     params.gap_cost >= kAnchoredCustomInf) {
    throw std::runtime_error("AnchoredCustomTrie costs are too large");
  }
}

template <int MaxBand>
auto AnchoredCustomTrie<MaxBand>::make_root_cache(
    const SearchState& st) const noexcept -> Cache {
  Cache out = inf_cache();
  // depth 0: cell(a) = editdist(query[0:a], target[0:0]).
  for(int a = 0; a <= st.radius && a <= st.query_len; ++a) {
    out[static_cast<size_t>(a + MaxBand + 1)] = gap_prefix_cost(a, st);
  }
  return out;
}

template <int MaxBand>
STARTREE_ALWAYS_INLINE
int AnchoredCustomTrie<MaxBand>::update_cache(
    const Cache& parent, const int target_code, const int depth,
    const SearchState& st, Cache* out) const noexcept {
  const int radius = st.radius;
  const int query_len = st.query_len;
  int lo = depth >= radius ? -radius : -depth;     // i = depth+a >= 0
  int hi = depth + radius <= query_len ? radius : query_len - depth;  // i <= query_len

  constexpr int center = MaxBand + 1;
  if(lo > hi) {
    // band empty: clear a sentinel so parent reads stay inf next step
    (*out)[static_cast<size_t>(center)] = st.inf;
    return kAnchoredCustomInf;
  }

  const int* p = parent.data();
  int* o = out->data();
  int col_min = st.inf;
  int left = st.inf;  // out(a-1); inf at a==lo
#ifdef ANCHORED_CUSTOM_TRIE_COUNT
  ++g_anchored_custom_columns;
  g_anchored_custom_cells += static_cast<size_t>(hi - lo + 1);
#endif
  for(int a = lo; a <= hi; ++a) {
    const int idx = a + center;
    const int qi = depth + a;  // query index i
    int best;
    if(qi >= 1) {
      best = add_cost(
        p[idx],
        substitution_cost(
          target_code,
          static_cast<int>(st.query[static_cast<size_t>(qi - 1)]),
          st
        ),
        st
      );
    } else {
      best = st.inf;
    }
    const int del = add_cost(p[idx + 1], st.gap_cost, st);  // parent(a+1) + gap
    if(del < best) best = del;
    const int ins = add_cost(left, st.gap_cost, st);         // out(a-1) + gap
    if(ins < best) best = ins;
    o[idx] = best;
    left = best;
    if(best < col_min) col_min = best;
#ifdef ANCHORED_CUSTOM_TRIE_COUNT
    if(best <= st.threshold) ++g_anchored_custom_live_cells;
#endif
  }
  // boundary sentinels (so reads at lo-1 / hi+1 by the next depth see inf);
  // the +1 sentinel slots make lo-1 == -(radius+1) and hi+1 == radius+1 in bounds.
  o[lo - 1 + center] = st.inf;
  o[hi + 1 + center] = st.inf;
  return col_min;
}

template <int MaxBand>
STARTREE_ALWAYS_INLINE
int AnchoredCustomTrie<MaxBand>::row_candidate(
    const Cache& c, const int depth, const SearchState& st) const noexcept {
  const int offset = st.query_len - depth;
  if(offset < -st.radius || offset > st.radius) {
    return kAnchoredCustomInf;
  }
  return cell(c, offset);
}

template <int MaxBand>
void AnchoredCustomTrie<MaxBand>::emit_hit(const uint32_t seq_idx,
                                             const int distance,
                                             const SearchState& st) const {
  if(seq_idx >= st.triangle_cutoff) {
    return;  // upper triangle (or self pair) in a self-join; skip
  }
  if(st.hits != nullptr) {
    st.hits->push_back(Hit{seq_idx, distance});
  } else {
    st.pairs->push_back(PairRecord{seq_idx, st.query_index, distance});
  }
}

template <int MaxBand>
void AnchoredCustomTrie<MaxBand>::add_all(const uint32_t node_idx,
                                            const int distance,
                                            const SearchState& st) const {
  const NodeType& node = nodes_[static_cast<size_t>(node_idx)];
  if(node.subtree_min_seq >= st.triangle_cutoff) {
    return;  // whole subtree is upper triangle; nothing to emit
  }
  if(node.sequence_index != kNoSequence) {
    emit_hit(node.sequence_index, distance, st);
  }
  unsigned mask = node.child_mask;
  while(mask) {
    const int code = ctz_u(mask);
    mask &= mask - 1;
    add_all(node.child[static_cast<size_t>(code)], distance, st);
  }
}

// True when the only band cell that is still <= threshold is the diagonal and
// the remaining budget cannot pay for the cheapest possible edit. The remaining
// search then reduces to an exact match of the query tail.
template <int MaxBand>
STARTREE_ALWAYS_INLINE
bool AnchoredCustomTrie<MaxBand>::can_dash(const Cache& c, const int depth,
                                            const SearchState& st) const noexcept {
  const int diagonal = cell(c, 0);
  if(diagonal > st.threshold ||
     st.threshold - diagonal >= st.min_edit_cost) {
    return false;
  }
  const int lo = depth >= st.radius ? -st.radius : -depth;
  const int hi = depth + st.radius <= st.query_len ? st.radius : st.query_len - depth;
  for(int a = lo; a <= hi; ++a) {
    if(a != 0 && cell(c, a) <= st.threshold) {
      return false;
    }
  }
  return true;
}

// Walk the query exactly through the trie from (node_idx, edge_i) at target_depth
// == depth. Every surviving alignment keeps the current diagonal distance: emit case-b at each
// terminal node passed (target consumed, query suffix free), and case-a (add_all)
// when the query is consumed (query aligned to a target prefix, target suffix free).
template <int MaxBand>
void AnchoredCustomTrie<MaxBand>::dash(const uint32_t node_idx, const int edge_i,
                                       const int depth, const int distance,
                                       const SearchState& st) const {
  uint32_t node = node_idx;
  int off = edge_i;
  int d = depth;
  while(d < st.query_len) {
    const NodeType& nd = nodes_[static_cast<size_t>(node)];
    const int elen = static_cast<int>(nd.edge.size());
    if(off == elen - 1 && nd.sequence_index != kNoSequence) {
      emit_hit(nd.sequence_index, distance, st);  // case (b)
    }
    const int qc = static_cast<int>(st.query[static_cast<size_t>(d)]);
    if(qc == kNoMatch) {
      return;  // query N can never match exactly
    }
    if(off + 1 < elen) {
      if(static_cast<int>(nd.edge[static_cast<size_t>(off + 1)]) != qc) {
        return;
      }
      ++off;
    } else {
      const uint32_t ch = nd.child[static_cast<size_t>(qc)];
      if(ch == kNoAnchoredCustomChild) {
        return;
      }
      node = ch;
      off = 0;
    }
    ++d;
  }
  add_all(node, distance, st);  // case (a): query fully aligned to target prefix
}

// Run the edge of node_idx starting at edge index edge_start, given the DP column
// (cache_in) and running row_min valid at depth_in (= just before consuming
// edge[edge_start]). Emits hits, seeds pebbles / logs case-b hits for the next
// query's LCP resume, then recurses into the node's children.
template <int MaxBand>
void AnchoredCustomTrie<MaxBand>::descend_edge(const uint32_t node_idx,
                                                 const int edge_start,
                                                 const Cache& cache_in,
                                                 const int depth_in,
                                                 const int row_min_in,
                                                 const SearchState& st) const {
  const NodeType& node = nodes_[static_cast<size_t>(node_idx)];
  // Lower-triangle prune: an entirely-upper-triangle subtree emits nothing.
  // Gate on depth_in >= seed_depth so shallow LCP seeding (depths <= seed_depth,
  // index-independent, needed by the higher-indexed next query) is never skipped;
  // for non-LCP search seed_depth == -1 so the prune applies throughout.
  if(depth_in >= st.seed_depth && node.subtree_min_seq >= st.triangle_cutoff) {
    return;
  }
  const int edge_len = static_cast<int>(node.edge.size());
  int depth = depth_in;
  int row_min = row_min_in;
  int col_min = kAnchoredCustomInf;

  Cache cache_a = cache_in;
  Cache cache_b;
  Cache* cur = &cache_a;
  Cache* nxt = &cache_b;

  for(int i = edge_start; i < edge_len; ++i) {
    ++depth;
    col_min = update_cache(*cur, node.edge[static_cast<size_t>(i)], depth, st, nxt);
    row_min = std::min(row_min, row_candidate(*nxt, depth, st));
    std::swap(cur, nxt);  // *cur is the column at this depth

    if(col_min > st.threshold && row_min > st.threshold) {
      return;
    }
    if(row_min <= st.threshold && row_min <= col_min) {
      add_all(node_idx, row_min, st);
      return;
    }
    // Exact-tail shortcut: once the edit budget is exhausted everywhere except the
    // diagonal, the only surviving alignment matches the query exactly from here.
    if(depth > st.seed_depth && depth < st.query_len && can_dash(*cur, depth, st)) {
      dash(node_idx, i, depth, cell(*cur, 0), st);
      return;
    }
    if(i < edge_len - 1 && depth > st.start_depth && depth <= st.seed_depth) {
      st.ws->pebbles[static_cast<size_t>(depth)].push_back(Pebble{node_idx, i, *cur});
    }
  }

  // node-end at target_depth == depth
  if(node.sequence_index != kNoSequence && col_min <= st.threshold) {
    emit_hit(node.sequence_index, col_min, st);
    if(depth > st.start_depth && depth <= st.seed_depth) {
      st.ws->hit_log[static_cast<size_t>(depth)].emplace_back(node.sequence_index, col_min);
    }
  }
  if(depth > st.start_depth && depth <= st.seed_depth) {
    st.ws->pebbles[static_cast<size_t>(depth)].push_back(Pebble{node_idx, edge_len - 1, *cur});
  }
  iterate_children(node_idx, *cur, depth, row_min, st);
}

template <int MaxBand>
void AnchoredCustomTrie<MaxBand>::iterate_children(const uint32_t node_idx,
                                                     const Cache& parent_cache,
                                                     const int parent_depth,
                                                     const int parent_row_min,
                                                     const SearchState& st) const {
  const NodeType& node = nodes_[static_cast<size_t>(node_idx)];
  unsigned mask = node.child_mask;
  while(mask) {
    const int code = ctz_u(mask);
    mask &= mask - 1;
    descend_edge(node.child[static_cast<size_t>(code)], 0, parent_cache,
                 parent_depth, parent_row_min, st);
  }
}

template <int MaxBand>
void AnchoredCustomTrie<MaxBand>::run_query(const SearchState& st,
                                              const Cache& root_cache,
                                              const int root_row_min) const {
  if(st.start_depth == 0) {
    // Fresh descent from the root.
    if(root_row_min <= st.threshold && root_row_min <= 0) {
      add_all(0, root_row_min, st);
      return;
    }
    if(nodes_[0].sequence_index != kNoSequence) {  // empty target, col_min(root)=0
      emit_hit(nodes_[0].sequence_index, 0, st);
    }
    iterate_children(0, root_cache, 0, root_row_min, st);
    return;
  }

  // Resume from the seeded frontier at start_depth.
  // 1) replay shallow case-b hits (depths 1..start) and the empty-target hit.
  if(nodes_[0].sequence_index != kNoSequence) {
    emit_hit(nodes_[0].sequence_index, 0, st);  // start_depth>0 => query long => no case-a
  }
  for(int d = 1; d <= st.start_depth; ++d) {
    for(const std::pair<uint32_t, int>& h : st.ws->hit_log[static_cast<size_t>(d)]) {
      emit_hit(h.first, h.second, st);
    }
  }
  // 2) continue the DFS from each frontier pebble (row_min == inf at this depth).
  for(const Pebble& peb : st.ws->pebbles[static_cast<size_t>(st.start_depth)]) {
    const NodeType& node = nodes_[static_cast<size_t>(peb.node_idx)];
    const int edge_len = static_cast<int>(node.edge.size());
    if(peb.edge_offset + 1 < edge_len) {
      descend_edge(peb.node_idx, peb.edge_offset + 1, peb.cache,
                   st.start_depth, kAnchoredCustomInf, st);
    } else {
      iterate_children(peb.node_idx, peb.cache, st.start_depth,
                       kAnchoredCustomInf, st);
    }
  }
}

template <int MaxBand>
void AnchoredCustomTrie<MaxBand>::search(const std::string& query_codes,
                                           const SearchParams& params,
                                           std::vector<Hit>* hits) const {
  validate_params(params);
  for(char ch : query_codes) {
    const int code = static_cast<unsigned char>(ch);
    if(code < 0 || code > kNoMatch) {
      throw std::runtime_error("AnchoredCustomTrie query code is out of range");
    }
  }
  SearchState st;
  st.query = reinterpret_cast<const unsigned char*>(query_codes.data());
  st.query_len = static_cast<int>(query_codes.size());
  st.threshold = params.max_distance;
  st.radius = params.band_radius;
  st.mismatch_cost = params.mismatch_cost;
  st.gap_cost = params.gap_cost;
  st.min_edit_cost = std::min(params.mismatch_cost, params.gap_cost);
  st.hits = hits;
  st.start_depth = 0;
  st.seed_depth = -1;  // no seeding

  const Cache root_cache = make_root_cache(st);
  const int root_row_min = st.query_len <= st.radius ? gap_prefix_cost(st.query_len, st)
                                                     : kAnchoredCustomInf;
  run_query(st, root_cache, root_row_min);
}

template <int MaxBand>
void AnchoredCustomTrie<MaxBand>::search_pairs(const std::string& query_codes,
                                                 const SearchParams& params,
                                                 const uint32_t query_index,
                                                 std::vector<PairRecord>* pairs) const {
  validate_params(params);
  for(char ch : query_codes) {
    const int code = static_cast<unsigned char>(ch);
    if(code < 0 || code > kNoMatch) {
      throw std::runtime_error("AnchoredCustomTrie query code is out of range");
    }
  }
  SearchState st;
  st.query = reinterpret_cast<const unsigned char*>(query_codes.data());
  st.query_len = static_cast<int>(query_codes.size());
  st.threshold = params.max_distance;
  st.radius = params.band_radius;
  st.mismatch_cost = params.mismatch_cost;
  st.gap_cost = params.gap_cost;
  st.min_edit_cost = std::min(params.mismatch_cost, params.gap_cost);
  st.query_index = query_index;
  st.pairs = pairs;
  st.start_depth = 0;
  st.seed_depth = -1;

  const Cache root_cache = make_root_cache(st);
  const int root_row_min = st.query_len <= st.radius ? gap_prefix_cost(st.query_len, st)
                                                     : kAnchoredCustomInf;
  run_query(st, root_cache, root_row_min);
}

template <int MaxBand>
void AnchoredCustomTrie<MaxBand>::search_batch_lcp(
    const std::vector<std::string>& query_codes, const SearchParams& params,
    const size_t begin, const size_t end, std::vector<PairRecord>* pairs,
    Workspace* ws, const bool lower_triangle) const {
  validate_params(params);
  if(begin >= end) {
    return;
  }
  const int radius = params.band_radius;

  int max_query_len = 0;
  for(size_t qi = begin; qi < end; ++qi) {
    max_query_len = std::max(max_query_len, static_cast<int>(query_codes[qi].size()));
  }
  // Reuse the workspace's per-level vectors across calls (blocks): grow the
  // outer vectors if needed, then clear() each in-range level so the inner
  // vectors keep their allocated capacity instead of reallocating per block.
  const size_t levels = static_cast<size_t>(max_query_len + 2);
  if(ws->pebbles.size() < levels) ws->pebbles.resize(levels);
  if(ws->hit_log.size() < levels) ws->hit_log.resize(levels);
  for(size_t d = 0; d < levels; ++d) {
    ws->pebbles[d].clear();
    ws->hit_log[d].clear();
  }

  int previous_seed = 0;
  for(size_t qi = begin; qi < end; ++qi) {
    const std::string& query = query_codes[qi];
    for(char ch : query) {
      const int code = static_cast<unsigned char>(ch);
      if(code < 0 || code > kNoMatch) {
        throw std::runtime_error("AnchoredCustomTrie query code is out of range");
      }
    }
    const int qlen = static_cast<int>(query.size());

    int start_depth = 0;
    if(qi > begin) {
      const int lcp = common_prefix_length(query_codes[qi - 1], query);
      start_depth = lcp - radius;
      start_depth = std::min(start_depth, qlen - radius - 1);
      start_depth = std::min(start_depth,
                             static_cast<int>(query_codes[qi - 1].size()) - radius - 1);
      start_depth = std::min(start_depth, previous_seed);
      start_depth = std::max(0, start_depth);
      while(start_depth > 0 && ws->pebbles[static_cast<size_t>(start_depth)].empty()) {
        --start_depth;
      }
    }

    int seed_depth = 0;
    if(qi + 1 < end) {
      const int lcp = common_prefix_length(query, query_codes[qi + 1]);
      seed_depth = lcp - radius;
      seed_depth = std::min(seed_depth, qlen - radius - 1);
      seed_depth = std::min(seed_depth,
                            static_cast<int>(query_codes[qi + 1].size()) - radius - 1);
      seed_depth = std::max(0, seed_depth);
    }
    previous_seed = seed_depth;

    // Clear the levels we are about to re-seed; reuse levels <= start_depth.
    for(int d = start_depth + 1; d <= seed_depth; ++d) {
      ws->pebbles[static_cast<size_t>(d)].clear();
      ws->hit_log[static_cast<size_t>(d)].clear();
    }

    SearchState st;
    st.query = reinterpret_cast<const unsigned char*>(query.data());
    st.query_len = qlen;
    st.threshold = params.max_distance;
    st.radius = radius;
    st.mismatch_cost = params.mismatch_cost;
    st.gap_cost = params.gap_cost;
    st.min_edit_cost = std::min(params.mismatch_cost, params.gap_cost);
    st.query_index = static_cast<uint32_t>(qi);
    st.pairs = pairs;
    st.ws = ws;
    st.start_depth = start_depth;
    st.seed_depth = seed_depth;
    if(lower_triangle) {
      st.triangle_cutoff = static_cast<uint32_t>(qi);
    }

    const Cache root_cache = make_root_cache(st);
    const int root_row_min = qlen <= radius ? gap_prefix_cost(qlen, st)
                                            : kAnchoredCustomInf;
    run_query(st, root_cache, root_row_min);
  }
}

}  // namespace radix
}  // namespace startree

#endif  // STARTREE_ANCHORED_CUSTOM_TRIE_H
