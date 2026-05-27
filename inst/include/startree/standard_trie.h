#ifndef STARTREE_STANDARD_TRIE_H
#define STARTREE_STANDARD_TRIE_H

#include "simple_array/small_array.h"
#include "startree/bits.h"
#include "startree/common.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
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

constexpr uint32_t kNoStandardChild = std::numeric_limits<uint32_t>::max();

template <int MaxBand>
using StandardCache = std::array<int8_t, 2 * MaxBand + 1>;

template <typename T>
using StandardSmallArr = trqwe::small_array<T, std::allocator<T>, uint16_t,
                                            std::integral_constant<uint16_t, 8>>;

template <int MaxBand>
constexpr StandardCache<MaxBand> standard_initial_cache() {
  StandardCache<MaxBand> c{};
  for(int i = -MaxBand; i <= MaxBand; ++i) {
    c[static_cast<size_t>(i + MaxBand)] =
      static_cast<int8_t>(i < 0 ? -i : i);
  }
  return c;
}

template <int MaxBand>
struct StandardNode {
  using Cache = StandardCache<MaxBand>;

  StandardSmallArr<Cache> caches;
  StandardSmallArr<uint8_t> edge;
  std::array<uint32_t, kAlphabet> child;
  uint32_t sequence_index = kNoSequence;
  uint8_t child_mask = 0;

  StandardNode() { child.fill(kNoStandardChild); }
};

struct StandardPebble {
  uint32_t node_idx;
  int32_t char_offset;
  uint32_t path;
};

enum class ForcedEdgeAdvance {
  Advanced,
  SwitchToNoPad,
  AtBranch,
  Done
};

enum class StandardPadMode {
  NoPad,
  PaddingAware
};

template <int MaxBand>
class StandardTrie {
  using NodeType = StandardNode<MaxBand>;
  using Cache = StandardCache<MaxBand>;

 public:
  void init(const int height, const int max_target_pad_depth = -1);
  void reserve_nodes(const size_t n);
  uint32_t insert_path(const std::string& seq);
  void set_sequence(const uint32_t leaf_idx, const uint32_t seq_idx);
  void search(const std::string& padded_query,
              const SearchParams& params,
              int start_depth,
              const int seed_depth,
              std::vector<Hit>* hits);

 private:
  uint32_t make_node();
  static const Cache& initial_cache_ref();
  static int cell(const Cache& c, int offset) noexcept;
  static void set_cell(Cache* c, int offset, int v) noexcept;
  static int substitution_cost(const int target_code,
                               const int query_code) noexcept;

  template <StandardPadMode Mode>
  static int target_gap_cost(const int target_code) noexcept;
  template <StandardPadMode Mode>
  static int query_gap_cost(const int query_code) noexcept;
  static int active_band(const int depth, const int threshold) noexcept;
  static bool path_band_touches_padding(const uint32_t path,
                                        const int maxa) noexcept;
  bool common_band_touches_padding(const uint32_t path,
                                   const int depth,
                                   const int maxa) const noexcept;

  template <StandardPadMode Mode>
  std::array<int, MaxBand + 1> compute_child_independent_dp(
    const Cache& fcache,
    const uint32_t fpath,
    const int depth
  ) const;

  template <StandardPadMode Mode>
  void compute_dp_update_and_cache(
    const Cache& fcache,
    const std::array<int, MaxBand + 1>& common,
    const int i,
    const int depth,
    Cache* ccache
  ) const;

  void push_hit(const NodeType& node, const int distance);
  void dash(uint32_t node_idx, int offset, int from_depth);
  template <StandardPadMode Mode>
  bool post_step(uint32_t node_idx,
                 int offset,
                 uint32_t path,
                 int depth,
                 const Cache& cache);
  template <StandardPadMode Mode>
  ForcedEdgeAdvance advance_forced_edge_character(uint32_t* fnode,
                                                  int* foffset,
                                                  uint32_t* fpath,
                                                  int* focus_depth);
  template <StandardPadMode Mode>
  void iterate_children(const uint32_t fnode,
                        const int foffset,
                        const uint32_t fpath,
                        const int focus_depth);
  template <StandardPadMode Mode>
  void recurse(uint32_t fnode, int foffset, uint32_t fpath, int focus_depth);

  int height_ = 0;
  int threshold_ = 0;
  int seed_depth_ = 0;
  int query_pad_depth_ = 0;
  std::vector<NodeType> nodes_;
  std::vector<std::vector<StandardPebble>> pebbles_;
  std::vector<int> translated_;
  std::vector<uint8_t> prefix_has_no_match_;
  std::vector<Hit>* hits_ = nullptr;
};

template <int MaxBand>
void StandardTrie<MaxBand>::init(const int height,
                                 const int /*max_target_pad_depth*/) {
  height_ = height;
  nodes_.clear();
  nodes_.emplace_back();
  pebbles_.assign(static_cast<size_t>(height + 2), {});
  pebbles_[0].push_back(StandardPebble{0, -1, 0});
}

template <int MaxBand>
void StandardTrie<MaxBand>::reserve_nodes(const size_t n) {
  nodes_.reserve(std::max<size_t>(n, 1));
}

template <int MaxBand>
uint32_t StandardTrie<MaxBand>::insert_path(const std::string& seq) {
  uint32_t node_idx = 0;
  int pos = 0;

  while(pos < height_) {
    const int c = static_cast<unsigned char>(seq[static_cast<size_t>(pos)]);
    const uint32_t child_idx = nodes_[node_idx].child[static_cast<size_t>(c)];

    if(child_idx == kNoStandardChild) {
      const uint32_t new_idx = make_node();
      NodeType& nn = nodes_[new_idx];
      const uint16_t len = static_cast<uint16_t>(height_ - pos);
      nn.edge = StandardSmallArr<uint8_t>(
        reinterpret_cast<const uint8_t*>(seq.data() + pos), len
      );
      nn.caches = StandardSmallArr<Cache>(len, standard_initial_cache<MaxBand>());
      nodes_[node_idx].child[static_cast<size_t>(c)] = new_idx;
      nodes_[node_idx].child_mask |= static_cast<uint8_t>(1U << c);
      return new_idx;
    }

    const NodeType& child = nodes_[child_idx];
    const int edge_len = static_cast<int>(child.edge.size());
    int matched = 0;
    while(matched < edge_len && pos + matched < height_ &&
          child.edge[static_cast<size_t>(matched)] ==
            static_cast<uint8_t>(seq[static_cast<size_t>(pos + matched)])) {
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
      tail.edge = StandardSmallArr<uint8_t>(ch.edge.data() + matched, suffix_len);
      tail.caches = StandardSmallArr<Cache>(ch.caches.data() + matched, suffix_len);
      tail.child = ch.child;
      tail.child_mask = ch.child_mask;
      tail.sequence_index = ch.sequence_index;
    }
    {
      NodeType& ch = nodes_[child_idx];
      ch.edge.resize(static_cast<uint16_t>(matched));
      ch.caches.resize(static_cast<uint16_t>(matched));
      ch.child.fill(kNoStandardChild);
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
    if(pos == height_) {
      return child_idx;
    }

    const uint32_t branch_idx = make_node();
    NodeType& branch = nodes_[branch_idx];
    const uint16_t blen = static_cast<uint16_t>(height_ - pos);
    branch.edge = StandardSmallArr<uint8_t>(
      reinterpret_cast<const uint8_t*>(seq.data() + pos), blen
    );
    branch.caches = StandardSmallArr<Cache>(blen, standard_initial_cache<MaxBand>());
    const int bc = branch.edge[0];
    nodes_[child_idx].child[static_cast<size_t>(bc)] = branch_idx;
    nodes_[child_idx].child_mask |= static_cast<uint8_t>(1U << bc);
    return branch_idx;
  }

  return node_idx;
}

template <int MaxBand>
void StandardTrie<MaxBand>::set_sequence(const uint32_t leaf_idx,
                                         const uint32_t seq_idx) {
  nodes_[static_cast<size_t>(leaf_idx)].sequence_index = seq_idx;
}

template <int MaxBand>
void StandardTrie<MaxBand>::search(const std::string& padded_query,
                                   const SearchParams& params,
                                   int start_depth,
                                   const int seed_depth,
                                   std::vector<Hit>* hits) {
  // Every entry below is written explicitly or in the loop, so resize (a no-op
  // that keeps stale values when height is unchanged across queries) is enough;
  // the prior full kEos/0 fill was redundant. translated_ holds [0]=height
  // sentinel, [1..height] query codes, and [height+1]=kEos terminator for dash().
  translated_.resize(static_cast<size_t>(height_ + 2));
  prefix_has_no_match_.resize(static_cast<size_t>(height_ + 1));
  translated_[0] = height_;
  translated_[static_cast<size_t>(height_ + 1)] = kEos;
  prefix_has_no_match_[0] = 0;
  static constexpr int kInsertToQuery[6] = {kNoMatch, 1, 2, 3, 4, kPad};
  query_pad_depth_ = 0;
  for(int i = 0; i < height_; ++i) {
    const int code = kInsertToQuery[static_cast<unsigned char>(
      padded_query[static_cast<size_t>(i)]
    )];
    if(code == kPad && query_pad_depth_ == i) {
      ++query_pad_depth_;
    }
    translated_[static_cast<size_t>(i + 1)] = code;
    prefix_has_no_match_[static_cast<size_t>(i + 1)] =
      static_cast<uint8_t>(prefix_has_no_match_[static_cast<size_t>(i)] ||
                           code == kNoMatch);
  }

  start_depth = std::max(start_depth, 0);
  start_depth = std::min(start_depth, height_ - 1);
  seed_depth_ = std::min(std::max(seed_depth, 0), height_);
  for(int d = start_depth + 1; d <= seed_depth_; ++d) {
    pebbles_[static_cast<size_t>(d)].clear();
  }

  hits_ = hits;
  threshold_ = params.max_distance;
  for(const StandardPebble& peb : pebbles_[static_cast<size_t>(start_depth)]) {
    recurse<StandardPadMode::PaddingAware>(
      peb.node_idx, peb.char_offset, peb.path, start_depth
    );
  }
  hits_ = nullptr;
}

template <int MaxBand>
uint32_t StandardTrie<MaxBand>::make_node() {
  nodes_.emplace_back();
  return static_cast<uint32_t>(nodes_.size() - 1);
}

template <int MaxBand>
STARTREE_ALWAYS_INLINE
auto StandardTrie<MaxBand>::initial_cache_ref() -> const Cache& {
  static constexpr Cache kInit = standard_initial_cache<MaxBand>();
  return kInit;
}

template <int MaxBand>
STARTREE_ALWAYS_INLINE
int StandardTrie<MaxBand>::cell(const Cache& c, int offset) noexcept {
  return c[static_cast<size_t>(offset + MaxBand)];
}

template <int MaxBand>
STARTREE_ALWAYS_INLINE
void StandardTrie<MaxBand>::set_cell(Cache* c, int offset, int v) noexcept {
  (*c)[static_cast<size_t>(offset + MaxBand)] = static_cast<int8_t>(v);
}

template <int MaxBand>
STARTREE_ALWAYS_INLINE
int StandardTrie<MaxBand>::substitution_cost(const int target_code,
                                             const int query_code) noexcept {
  return target_code == query_code ? 0 : 1;
}

template <int MaxBand>
template <StandardPadMode Mode>
STARTREE_ALWAYS_INLINE
int StandardTrie<MaxBand>::target_gap_cost(const int target_code) noexcept {
  if constexpr(Mode == StandardPadMode::PaddingAware) {
    return target_code == kPad ? 0 : 1;
  }
  (void)target_code;
  return 1;
}

template <int MaxBand>
template <StandardPadMode Mode>
STARTREE_ALWAYS_INLINE
int StandardTrie<MaxBand>::query_gap_cost(const int query_code) noexcept {
  if constexpr(Mode == StandardPadMode::PaddingAware) {
    return query_code == kPad ? 0 : 1;
  }
  (void)query_code;
  return 1;
}

template <int MaxBand>
STARTREE_ALWAYS_INLINE
int StandardTrie<MaxBand>::active_band(const int depth,
                                       const int threshold) noexcept {
  return std::min(depth - 1, std::min(threshold, MaxBand));
}

template <int MaxBand>
STARTREE_ALWAYS_INLINE
bool StandardTrie<MaxBand>::path_band_touches_padding(
    const uint32_t path,
    const int maxa) noexcept {
  if(maxa <= 0) {
    return false;
  }
  // Padding is a leading run, so if the oldest character in the active
  // target band is not PAD, no newer character in that band can be PAD.
  return static_cast<int>((path >> (4 * (maxa - 1))) & 15U) == kPad;
}

template <int MaxBand>
STARTREE_ALWAYS_INLINE
bool StandardTrie<MaxBand>::common_band_touches_padding(
    const uint32_t path,
    const int depth,
    const int maxa) const noexcept {
  return maxa > 0 &&
    (depth <= query_pad_depth_ + 1 ||
     path_band_touches_padding(path, maxa));
}

template <int MaxBand>
template <StandardPadMode Mode>
STARTREE_ALWAYS_INLINE
std::array<int, MaxBand + 1>
StandardTrie<MaxBand>::compute_child_independent_dp(
    const Cache& fcache,
    const uint32_t fpath,
    const int depth) const {
  std::array<int, MaxBand + 1> common{};
  for(int i = 0; i <= MaxBand; ++i) {
    common[static_cast<size_t>(i)] = i + 1;
  }

  const int maxa = active_band(depth, threshold_);
  if(maxa <= 0) {
    return common;
  }

  const int query_code = translated_[static_cast<size_t>(depth)];
  const int query_gap = query_gap_cost<Mode>(query_code);
  if constexpr(Mode == StandardPadMode::PaddingAware) {
    if(maxa == depth - 1) {
      common[static_cast<size_t>(maxa)] = cell(fcache, maxa) + query_gap;
    }
  }

  int target_code = static_cast<int>((fpath >> (4 * (maxa - 1))) & 15U);
  int mmatch = cell(fcache, maxa) + substitution_cost(target_code, query_code);
  int shift;
  if constexpr(Mode == StandardPadMode::PaddingAware) {
    shift = std::min(
      cell(fcache, maxa - 1) + query_gap,
      common[static_cast<size_t>(maxa)] + target_gap_cost<Mode>(target_code)
    );
  } else {
    shift =
      std::min(cell(fcache, maxa - 1), common[static_cast<size_t>(maxa)]) +
      query_gap;
  }
  common[static_cast<size_t>(maxa - 1)] = std::min(mmatch, shift);
  for(int a = maxa - 1; a > 0; --a) {
    target_code = static_cast<int>((fpath >> (4 * (a - 1))) & 15U);
    mmatch = cell(fcache, a) + substitution_cost(target_code, query_code);
    if constexpr(Mode == StandardPadMode::PaddingAware) {
      shift = std::min(
        cell(fcache, a - 1) + query_gap,
        common[static_cast<size_t>(a)] + target_gap_cost<Mode>(target_code)
      );
    } else {
      shift = std::min(cell(fcache, a - 1), common[static_cast<size_t>(a)]) +
        query_gap;
    }
    common[static_cast<size_t>(a - 1)] = std::min(mmatch, shift);
  }
  return common;
}

template <int MaxBand>
template <StandardPadMode Mode>
STARTREE_ALWAYS_INLINE
void StandardTrie<MaxBand>::compute_dp_update_and_cache(
    const Cache& fcache,
    const std::array<int, MaxBand + 1>& common,
    const int i,
    const int depth,
    Cache* ccache) const {
  *ccache = initial_cache_ref();
  for(int a = 1; a <= MaxBand; ++a) {
    set_cell(ccache, a, common[static_cast<size_t>(a - 1)]);
  }

  if constexpr(Mode == StandardPadMode::PaddingAware) {
    const int target_gap = target_gap_cost<Mode>(i);
    if(depth <= MaxBand) {
      set_cell(ccache, -depth, cell(fcache, 1 - depth) + target_gap);
      set_cell(ccache, depth, cell(fcache, depth - 1) +
        query_gap_cost<Mode>(translated_[static_cast<size_t>(depth)]));
    }
  }

  const int maxa = active_band(depth, threshold_);
  if(maxa > 0) {
    int query_code = translated_[static_cast<size_t>(depth - maxa)];
    int mmatch = cell(fcache, -maxa) + substitution_cost(i, query_code);
    int shift;
    if constexpr(Mode == StandardPadMode::PaddingAware) {
      const int target_gap = target_gap_cost<Mode>(i);
      const int insert_from =
        maxa < MaxBand ? cell(*ccache, -maxa - 1) : maxa + 1;
      shift = std::min(
        cell(fcache, 1 - maxa) + target_gap,
        insert_from + query_gap_cost<Mode>(query_code)
      );
    } else {
      shift = std::min(cell(fcache, 1 - maxa), maxa + 1) + 1;
    }
    set_cell(ccache, -maxa, std::min(mmatch, shift));

    if constexpr(MaxBand > 1) {
      for(int a = maxa - 1; a > 0; --a) {
        query_code = translated_[static_cast<size_t>(depth - a)];
        mmatch = cell(fcache, -a) + substitution_cost(i, query_code);
        if constexpr(Mode == StandardPadMode::PaddingAware) {
          const int target_gap = target_gap_cost<Mode>(i);
          shift = std::min(
            cell(fcache, 1 - a) + target_gap,
            cell(*ccache, -a - 1) + query_gap_cost<Mode>(query_code)
          );
        } else {
          shift = std::min(cell(fcache, 1 - a), cell(*ccache, -a - 1)) + 1;
        }
        set_cell(ccache, -a, std::min(mmatch, shift));
      }
    }
  }

  const int query_code = translated_[static_cast<size_t>(depth)];
  const int mmatch = cell(fcache, 0) + substitution_cost(i, query_code);
  int shift;
  if constexpr(Mode == StandardPadMode::PaddingAware) {
    shift = std::min(
      cell(*ccache, -1) + query_gap_cost<Mode>(query_code),
      cell(*ccache, 1) + target_gap_cost<Mode>(i)
    );
  } else {
    shift = std::min(cell(*ccache, -1), cell(*ccache, 1)) + 1;
  }
  set_cell(ccache, 0, std::min(mmatch, shift));
}

template <int MaxBand>
STARTREE_ALWAYS_INLINE
void StandardTrie<MaxBand>::push_hit(const NodeType& node,
                                     const int distance) {
  if(node.sequence_index != kNoSequence) {
    hits_->push_back(Hit{node.sequence_index, distance});
  }
}

template <int MaxBand>
void StandardTrie<MaxBand>::dash(uint32_t node_idx,
                                 int offset,
                                 int from_depth) {
  size_t sfx = static_cast<size_t>(from_depth + 1);
  while(sfx < translated_.size() && translated_[sfx] != kEos) {
    const int qc = translated_[sfx++];
    if(qc == kNoMatch) {
      return;
    }
    const NodeType* node = &nodes_[static_cast<size_t>(node_idx)];
    if(offset + 1 < static_cast<int>(node->edge.size())) {
      const int ec = node->edge[static_cast<size_t>(offset + 1)];
      if(ec != qc) {
        return;
      }
      ++offset;
    } else {
      const uint32_t ch = node->child[static_cast<size_t>(qc)];
      if(ch == kNoStandardChild) {
        return;
      }
      node_idx = ch;
      offset = 0;
    }
  }
  push_hit(nodes_[static_cast<size_t>(node_idx)], threshold_);
}

template <int MaxBand>
template <StandardPadMode Mode>
STARTREE_ALWAYS_INLINE
bool StandardTrie<MaxBand>::post_step(uint32_t node_idx,
                                      int offset,
                                      uint32_t path,
                                      int depth,
                                      const Cache& cache) {
  if(depth == height_) {
    push_hit(nodes_[static_cast<size_t>(node_idx)], cell(cache, 0));
    return false;
  }
  if(depth <= seed_depth_) {
    pebbles_[static_cast<size_t>(depth)].push_back(
      StandardPebble{node_idx, offset, path}
    );
  }
  if constexpr(Mode == StandardPadMode::NoPad) {
    if(threshold_ <= 3 && !prefix_has_no_match_[static_cast<size_t>(depth)] &&
       depth > seed_depth_) {
      bool can_dash = true;
      const int maxa = std::min(depth - 1, std::min(threshold_, MaxBand));
      for(int a = -maxa; a < maxa + 1; ++a) {
        if(cell(cache, a) < threshold_) {
          can_dash = false;
          break;
        }
      }
      if(can_dash) {
        dash(node_idx, offset, depth);
        return false;
      }
    }
  }
  return true;
}

template <int MaxBand>
template <StandardPadMode Mode>
STARTREE_ALWAYS_INLINE
ForcedEdgeAdvance StandardTrie<MaxBand>::advance_forced_edge_character(
    uint32_t* fnode,
    int* foffset,
    uint32_t* fpath,
    int* focus_depth) {
  const NodeType& fnode_ref = nodes_[static_cast<size_t>(*fnode)];
  const int edge_len = static_cast<int>(fnode_ref.edge.size());
  if(*foffset + 1 >= edge_len) {
    return ForcedEdgeAdvance::AtBranch;
  }

  const Cache& fcache =
    (*foffset >= 0) ? fnode_ref.caches[static_cast<size_t>(*foffset)]
                    : initial_cache_ref();
  const int next_depth = *focus_depth + 1;
  const int next_off = *foffset + 1;
  const int code = fnode_ref.edge[static_cast<size_t>(next_off)];
  Cache& ncache = nodes_[static_cast<size_t>(*fnode)]
    .caches[static_cast<size_t>(next_off)];

  ForcedEdgeAdvance next_advance = ForcedEdgeAdvance::Advanced;
  bool child_touches_padding = false;
  if constexpr(Mode == StandardPadMode::NoPad) {
    const auto common =
      compute_child_independent_dp<StandardPadMode::NoPad>(
        fcache, *fpath, next_depth
      );
    compute_dp_update_and_cache<StandardPadMode::NoPad>(
      fcache, common, code, next_depth, &ncache
    );
  } else {
    const int maxa = active_band(next_depth, threshold_);
    const bool common_touches_padding =
      common_band_touches_padding(*fpath, next_depth, maxa);
    const auto common = common_touches_padding
      ? compute_child_independent_dp<StandardPadMode::PaddingAware>(
          fcache, *fpath, next_depth
        )
      : compute_child_independent_dp<StandardPadMode::NoPad>(
          fcache, *fpath, next_depth
        );

    child_touches_padding =
      common_touches_padding || code == kPad ||
      next_depth - maxa <= query_pad_depth_;
    if(child_touches_padding) {
      compute_dp_update_and_cache<StandardPadMode::PaddingAware>(
        fcache, common, code, next_depth, &ncache
      );
    } else {
      compute_dp_update_and_cache<StandardPadMode::NoPad>(
        fcache, common, code, next_depth, &ncache
      );
      next_advance = ForcedEdgeAdvance::SwitchToNoPad;
    }
  }
  if(cell(ncache, 0) > threshold_) {
    return ForcedEdgeAdvance::Done;
  }

  const uint32_t npath = (*fpath << 4) | static_cast<uint32_t>(code & 0xF);
  bool keep_going;
  if constexpr(Mode == StandardPadMode::NoPad) {
    keep_going =
      post_step<StandardPadMode::NoPad>(
        *fnode, next_off, npath, next_depth, ncache
      );
  } else {
    keep_going = child_touches_padding
      ? post_step<StandardPadMode::PaddingAware>(
          *fnode, next_off, npath, next_depth, ncache
        )
      : post_step<StandardPadMode::NoPad>(
          *fnode, next_off, npath, next_depth, ncache
        );
  }
  if(!keep_going) {
    return ForcedEdgeAdvance::Done;
  }

  *foffset = next_off;
  *fpath = npath;
  *focus_depth = next_depth;
  return next_advance;
}

template <int MaxBand>
template <StandardPadMode Mode>
void StandardTrie<MaxBand>::iterate_children(const uint32_t fnode,
                                             const int foffset,
                                             const uint32_t fpath,
                                             const int focus_depth) {
  const NodeType& fnode_ref = nodes_[static_cast<size_t>(fnode)];
  const Cache& fcache =
    (foffset >= 0) ? fnode_ref.caches[static_cast<size_t>(foffset)]
                   : initial_cache_ref();
  const int next_depth = focus_depth + 1;

  unsigned mask = fnode_ref.child_mask;
  if constexpr(Mode == StandardPadMode::NoPad) {
    const auto common =
      compute_child_independent_dp<StandardPadMode::NoPad>(
        fcache, fpath, next_depth
      );
    while(mask) {
      const int code = ctz_u(mask);
      mask &= mask - 1;
      const uint32_t child_idx = fnode_ref.child[static_cast<size_t>(code)];
      Cache& ncache = nodes_[static_cast<size_t>(child_idx)].caches[0];
      compute_dp_update_and_cache<StandardPadMode::NoPad>(
        fcache, common, code, next_depth, &ncache
      );
      if(cell(ncache, 0) > threshold_) {
        continue;
      }
      const uint32_t npath = (fpath << 4) | static_cast<uint32_t>(code & 0xF);
      if(!post_step<StandardPadMode::NoPad>(
          child_idx, 0, npath, next_depth, ncache)) {
        continue;
      }
      recurse<StandardPadMode::NoPad>(child_idx, 0, npath, next_depth);
    }
  } else {
    const int maxa = active_band(next_depth, threshold_);
    const bool common_touches_padding =
      common_band_touches_padding(fpath, next_depth, maxa);
    const auto common = common_touches_padding
      ? compute_child_independent_dp<StandardPadMode::PaddingAware>(
          fcache, fpath, next_depth
        )
      : compute_child_independent_dp<StandardPadMode::NoPad>(
          fcache, fpath, next_depth
        );
    const bool child_independent_touches_padding =
      common_touches_padding || next_depth - maxa <= query_pad_depth_;

    while(mask) {
      const int code = ctz_u(mask);
      mask &= mask - 1;
      const uint32_t child_idx = fnode_ref.child[static_cast<size_t>(code)];
      Cache& ncache = nodes_[static_cast<size_t>(child_idx)].caches[0];
      const bool child_touches_padding =
        child_independent_touches_padding || code == kPad;
      if(child_touches_padding) {
        compute_dp_update_and_cache<StandardPadMode::PaddingAware>(
          fcache, common, code, next_depth, &ncache
        );
      } else {
        compute_dp_update_and_cache<StandardPadMode::NoPad>(
          fcache, common, code, next_depth, &ncache
        );
      }
      if(cell(ncache, 0) > threshold_) {
        continue;
      }
      const uint32_t npath = (fpath << 4) | static_cast<uint32_t>(code & 0xF);
      const bool keep_going = child_touches_padding
        ? post_step<StandardPadMode::PaddingAware>(
            child_idx, 0, npath, next_depth, ncache
          )
        : post_step<StandardPadMode::NoPad>(
            child_idx, 0, npath, next_depth, ncache
          );
      if(!keep_going) {
        continue;
      }
      if(child_touches_padding) {
        recurse<StandardPadMode::PaddingAware>(child_idx, 0, npath, next_depth);
      } else {
        recurse<StandardPadMode::NoPad>(child_idx, 0, npath, next_depth);
      }
    }
  }
}

template <int MaxBand>
template <StandardPadMode Mode>
void StandardTrie<MaxBand>::recurse(uint32_t fnode,
                                    int foffset,
                                    uint32_t fpath,
                                    int focus_depth) {
  while(true) {
    const ForcedEdgeAdvance advance =
      advance_forced_edge_character<Mode>(
        &fnode, &foffset, &fpath, &focus_depth
      );
    if(advance == ForcedEdgeAdvance::Advanced) {
      continue;
    }
    if constexpr(Mode == StandardPadMode::PaddingAware) {
      if(advance == ForcedEdgeAdvance::SwitchToNoPad) {
        recurse<StandardPadMode::NoPad>(fnode, foffset, fpath, focus_depth);
        return;
      }
    }
    if(advance == ForcedEdgeAdvance::AtBranch) {
      iterate_children<Mode>(fnode, foffset, fpath, focus_depth);
    }
    return;
  }
}

}  // namespace radix
}  // namespace startree

#endif  // STARTREE_STANDARD_TRIE_H
