#ifndef STARTREE_CUSTOM_TRIE_H
#define STARTREE_CUSTOM_TRIE_H

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

constexpr uint32_t kNoCustomChild = std::numeric_limits<uint32_t>::max();

template <int MaxBand>
using CustomCache = std::array<int, 2 * MaxBand + 1>;

template <typename T, uint16_t StackSize>
using CustomSmallArrBase =
  trqwe::small_array<T, std::allocator<T>, uint16_t,
                     std::integral_constant<uint16_t, StackSize>>;

template <typename T>
using CustomSmallArr = CustomSmallArrBase<T, 8>;

// Weighted DP caches are much larger than edge bytes. Keeping only two caches
// inline keeps radix nodes cache-friendly; longer compressed edges spill.
template <int MaxBand>
using CustomCacheArr = CustomSmallArrBase<CustomCache<MaxBand>, 2>;

template <int MaxBand>
constexpr CustomCache<MaxBand> custom_initial_cache(const int gap_cost = 1) {
  CustomCache<MaxBand> c{};
  for(int i = -MaxBand; i <= MaxBand; ++i) {
    c[static_cast<size_t>(i + MaxBand)] = (i < 0 ? -i : i) * gap_cost;
  }
  return c;
}

template <int MaxBand>
struct CustomNode {
  using Cache = CustomCache<MaxBand>;

  CustomCacheArr<MaxBand> caches;
  CustomSmallArr<uint8_t> edge;
  std::array<uint32_t, kAlphabet> child;
  uint32_t sequence_index = kNoSequence;
  uint8_t child_mask = 0;

  CustomNode() { child.fill(kNoCustomChild); }
};

struct CustomPebble {
  uint32_t node_idx;
  int32_t char_offset;
  uint32_t path;
};

enum class CustomEdgeAdvance {
  Advanced,
  SwitchToNoPad,
  AtBranch,
  Done
};

enum class CustomPadMode {
  NoPad,
  PaddingAware
};

template <int MaxBand>
class CustomTrie {
  using NodeType = CustomNode<MaxBand>;
  using Cache = CustomCache<MaxBand>;

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
  static int cell(const Cache& c, int offset) noexcept;
  static void set_cell(Cache* c, int offset, int v) noexcept;

  template <CustomPadMode Mode>
  int substitution_cost(const int target_code,
                        const int query_code) const noexcept;
  template <CustomPadMode Mode>
  int target_gap_cost(const int target_code) const noexcept;
  template <CustomPadMode Mode>
  int query_gap_cost(const int query_code) const noexcept;
  int active_band(const int depth) const noexcept;
  static bool path_band_touches_padding(const uint32_t path,
                                        const int maxa) noexcept;
  bool common_band_touches_padding(const uint32_t path,
                                   const int depth,
                                   const int maxa) const noexcept;

  template <CustomPadMode Mode>
  std::array<int, MaxBand + 1> compute_child_independent_dp(
    const Cache& fcache,
    const uint32_t fpath,
    const int depth
  ) const;

  template <CustomPadMode Mode>
  void compute_dp_update_and_cache(
    const Cache& fcache,
    const std::array<int, MaxBand + 1>& common,
    const int i,
    const int depth,
    Cache* ccache
  ) const;

  void push_hit(const NodeType& node, const int distance);
  void dash(uint32_t node_idx, int offset, int from_depth);
  template <CustomPadMode Mode>
  bool post_step(uint32_t node_idx,
                 int offset,
                 uint32_t path,
                 int depth,
                 const Cache& cache);
  template <CustomPadMode Mode>
  CustomEdgeAdvance advance_forced_edge_character(uint32_t* fnode,
                                                  int* foffset,
                                                  uint32_t* fpath,
                                                  int* focus_depth);
  template <CustomPadMode Mode>
  void iterate_children(const uint32_t fnode,
                        const int foffset,
                        const uint32_t fpath,
                        const int focus_depth);
  template <CustomPadMode Mode>
  void recurse(uint32_t fnode, int foffset, uint32_t fpath, int focus_depth);

  int height_ = 0;
  int threshold_ = 0;
  int band_radius_ = 0;
  int seed_depth_ = 0;
  int query_pad_depth_ = 0;
  int mismatch_cost_ = 1;
  int gap_cost_ = 1;
  Cache initial_cache_ = custom_initial_cache<MaxBand>();
  std::vector<NodeType> nodes_;
  std::vector<std::vector<CustomPebble>> pebbles_;
  std::vector<int> translated_;
  std::vector<uint8_t> prefix_has_no_match_;
  std::vector<Hit>* hits_ = nullptr;
};

template <int MaxBand>
void CustomTrie<MaxBand>::init(const int height,
                               const int /*max_target_pad_depth*/) {
  height_ = height;
  nodes_.clear();
  nodes_.emplace_back();
  pebbles_.assign(static_cast<size_t>(height + 2), {});
  pebbles_[0].push_back(CustomPebble{0, -1, 0});
}

template <int MaxBand>
void CustomTrie<MaxBand>::reserve_nodes(const size_t n) {
  nodes_.reserve(std::max<size_t>(n, 1));
}

template <int MaxBand>
uint32_t CustomTrie<MaxBand>::insert_path(const std::string& seq) {
  uint32_t node_idx = 0;
  int pos = 0;

  while(pos < height_) {
    const int c = static_cast<unsigned char>(seq[static_cast<size_t>(pos)]);
    const uint32_t child_idx = nodes_[node_idx].child[static_cast<size_t>(c)];

    if(child_idx == kNoCustomChild) {
      const uint32_t new_idx = make_node();
      NodeType& nn = nodes_[new_idx];
      const uint16_t len = static_cast<uint16_t>(height_ - pos);
      nn.edge = CustomSmallArr<uint8_t>(
        reinterpret_cast<const uint8_t*>(seq.data() + pos), len
      );
      nn.caches = CustomCacheArr<MaxBand>(len, custom_initial_cache<MaxBand>());
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
      tail.edge = CustomSmallArr<uint8_t>(ch.edge.data() + matched, suffix_len);
      tail.caches = CustomCacheArr<MaxBand>(ch.caches.data() + matched, suffix_len);
      tail.child = ch.child;
      tail.child_mask = ch.child_mask;
      tail.sequence_index = ch.sequence_index;
    }
    {
      NodeType& ch = nodes_[child_idx];
      ch.edge.resize(static_cast<uint16_t>(matched));
      ch.caches.resize(static_cast<uint16_t>(matched));
      ch.child.fill(kNoCustomChild);
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
    branch.edge = CustomSmallArr<uint8_t>(
      reinterpret_cast<const uint8_t*>(seq.data() + pos), blen
    );
    branch.caches = CustomCacheArr<MaxBand>(blen, custom_initial_cache<MaxBand>());
    const int bc = branch.edge[0];
    nodes_[child_idx].child[static_cast<size_t>(bc)] = branch_idx;
    nodes_[child_idx].child_mask |= static_cast<uint8_t>(1U << bc);
    return branch_idx;
  }

  return node_idx;
}

template <int MaxBand>
void CustomTrie<MaxBand>::set_sequence(const uint32_t leaf_idx,
                                       const uint32_t seq_idx) {
  nodes_[static_cast<size_t>(leaf_idx)].sequence_index = seq_idx;
}

template <int MaxBand>
void CustomTrie<MaxBand>::search(const std::string& padded_query,
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
  band_radius_ = params.band_radius;
  mismatch_cost_ = params.mismatch_cost;
  gap_cost_ = params.gap_cost;
  initial_cache_ = custom_initial_cache<MaxBand>(gap_cost_);
  for(const CustomPebble& peb : pebbles_[static_cast<size_t>(start_depth)]) {
    recurse<CustomPadMode::PaddingAware>(
      peb.node_idx, peb.char_offset, peb.path, start_depth
    );
  }
  hits_ = nullptr;
}

template <int MaxBand>
uint32_t CustomTrie<MaxBand>::make_node() {
  nodes_.emplace_back();
  return static_cast<uint32_t>(nodes_.size() - 1);
}

template <int MaxBand>
STARTREE_ALWAYS_INLINE
int CustomTrie<MaxBand>::cell(const Cache& c, int offset) noexcept {
  return c[static_cast<size_t>(offset + MaxBand)];
}

template <int MaxBand>
STARTREE_ALWAYS_INLINE
void CustomTrie<MaxBand>::set_cell(Cache* c, int offset, int v) noexcept {
  (*c)[static_cast<size_t>(offset + MaxBand)] = v;
}

template <int MaxBand>
template <CustomPadMode Mode>
STARTREE_ALWAYS_INLINE
int CustomTrie<MaxBand>::substitution_cost(const int target_code,
                                           const int query_code) const noexcept {
  if(target_code == query_code) {
    return 0;
  }
  if constexpr(Mode == CustomPadMode::PaddingAware) {
    if(target_code == kPad || query_code == kPad) {
      return gap_cost_;
    }
  }
  return mismatch_cost_;
}

template <int MaxBand>
template <CustomPadMode Mode>
STARTREE_ALWAYS_INLINE
int CustomTrie<MaxBand>::target_gap_cost(const int target_code) const noexcept {
  if constexpr(Mode == CustomPadMode::PaddingAware) {
    return target_code == kPad ? 0 : gap_cost_;
  }
  (void)target_code;
  return gap_cost_;
}

template <int MaxBand>
template <CustomPadMode Mode>
STARTREE_ALWAYS_INLINE
int CustomTrie<MaxBand>::query_gap_cost(const int query_code) const noexcept {
  if constexpr(Mode == CustomPadMode::PaddingAware) {
    return query_code == kPad ? 0 : gap_cost_;
  }
  (void)query_code;
  return gap_cost_;
}

template <int MaxBand>
STARTREE_ALWAYS_INLINE
int CustomTrie<MaxBand>::active_band(const int depth) const noexcept {
  return std::min(depth - 1, std::min(band_radius_, MaxBand));
}

template <int MaxBand>
STARTREE_ALWAYS_INLINE
bool CustomTrie<MaxBand>::path_band_touches_padding(
    const uint32_t path,
    const int maxa) noexcept {
  if(maxa <= 0) {
    return false;
  }
  return static_cast<int>((path >> (4 * (maxa - 1))) & 15U) == kPad;
}

template <int MaxBand>
STARTREE_ALWAYS_INLINE
bool CustomTrie<MaxBand>::common_band_touches_padding(
    const uint32_t path,
    const int depth,
    const int maxa) const noexcept {
  return maxa > 0 &&
    (depth <= query_pad_depth_ + 1 ||
     path_band_touches_padding(path, maxa));
}

template <int MaxBand>
template <CustomPadMode Mode>
STARTREE_ALWAYS_INLINE
std::array<int, MaxBand + 1>
CustomTrie<MaxBand>::compute_child_independent_dp(
    const Cache& fcache,
    const uint32_t fpath,
    const int depth) const {
  std::array<int, MaxBand + 1> common{};
  for(int i = 0; i <= MaxBand; ++i) {
    common[static_cast<size_t>(i)] = (i + 1) * gap_cost_;
  }

  const int maxa = active_band(depth);
  if(maxa <= 0) {
    return common;
  }

  const int query_code = translated_[static_cast<size_t>(depth)];
  const int query_gap = query_gap_cost<Mode>(query_code);
  if constexpr(Mode == CustomPadMode::PaddingAware) {
    if(maxa == depth - 1) {
      common[static_cast<size_t>(maxa)] = cell(fcache, maxa) + query_gap;
    }
  }

  int target_code = static_cast<int>((fpath >> (4 * (maxa - 1))) & 15U);
  int mmatch =
    cell(fcache, maxa) + substitution_cost<Mode>(target_code, query_code);
  int shift;
  if constexpr(Mode == CustomPadMode::PaddingAware) {
    shift = std::min(
      cell(fcache, maxa - 1) + query_gap,
      common[static_cast<size_t>(maxa)] + target_gap_cost<Mode>(target_code)
    );
  } else {
    shift =
      std::min(cell(fcache, maxa - 1), common[static_cast<size_t>(maxa)]) +
      gap_cost_;
  }
  common[static_cast<size_t>(maxa - 1)] = std::min(mmatch, shift);
  for(int a = maxa - 1; a > 0; --a) {
    target_code = static_cast<int>((fpath >> (4 * (a - 1))) & 15U);
    mmatch =
      cell(fcache, a) + substitution_cost<Mode>(target_code, query_code);
    if constexpr(Mode == CustomPadMode::PaddingAware) {
      shift = std::min(
        cell(fcache, a - 1) + query_gap,
        common[static_cast<size_t>(a)] + target_gap_cost<Mode>(target_code)
      );
    } else {
      shift = std::min(cell(fcache, a - 1), common[static_cast<size_t>(a)]) +
        gap_cost_;
    }
    common[static_cast<size_t>(a - 1)] = std::min(mmatch, shift);
  }
  return common;
}

template <int MaxBand>
template <CustomPadMode Mode>
STARTREE_ALWAYS_INLINE
void CustomTrie<MaxBand>::compute_dp_update_and_cache(
    const Cache& fcache,
    const std::array<int, MaxBand + 1>& common,
    const int i,
    const int depth,
    Cache* ccache) const {
  for(int a = 1; a <= MaxBand; ++a) {
    set_cell(ccache, a, common[static_cast<size_t>(a - 1)]);
  }

  const int maxa = active_band(depth);
  if(maxa < MaxBand) {
    set_cell(ccache, -maxa - 1, (maxa + 1) * gap_cost_);
  }

  if constexpr(Mode == CustomPadMode::PaddingAware) {
    const int target_gap = target_gap_cost<Mode>(i);
    if(depth <= MaxBand) {
      set_cell(ccache, -depth, cell(fcache, 1 - depth) + target_gap);
      set_cell(ccache, depth, cell(fcache, depth - 1) +
        query_gap_cost<Mode>(translated_[static_cast<size_t>(depth)]));
    }
  }

  if(maxa > 0) {
    int query_code = translated_[static_cast<size_t>(depth - maxa)];
    int mmatch = cell(fcache, -maxa) + substitution_cost<Mode>(i, query_code);
    int shift;
    if constexpr(Mode == CustomPadMode::PaddingAware) {
      const int target_gap = target_gap_cost<Mode>(i);
      const int insert_from =
        maxa < MaxBand ? cell(*ccache, -maxa - 1) : (maxa + 1) * gap_cost_;
      shift = std::min(
        cell(fcache, 1 - maxa) + target_gap,
        insert_from + query_gap_cost<Mode>(query_code)
      );
    } else {
      shift =
        std::min(cell(fcache, 1 - maxa), (maxa + 1) * gap_cost_) + gap_cost_;
    }
    set_cell(ccache, -maxa, std::min(mmatch, shift));

    if constexpr(MaxBand > 1) {
      for(int a = maxa - 1; a > 0; --a) {
        query_code = translated_[static_cast<size_t>(depth - a)];
        mmatch = cell(fcache, -a) + substitution_cost<Mode>(i, query_code);
        if constexpr(Mode == CustomPadMode::PaddingAware) {
          const int target_gap = target_gap_cost<Mode>(i);
          shift = std::min(
            cell(fcache, 1 - a) + target_gap,
            cell(*ccache, -a - 1) + query_gap_cost<Mode>(query_code)
          );
        } else {
          shift = std::min(cell(fcache, 1 - a), cell(*ccache, -a - 1)) +
            gap_cost_;
        }
        set_cell(ccache, -a, std::min(mmatch, shift));
      }
    }
  }

  const int query_code = translated_[static_cast<size_t>(depth)];
  const int mmatch = cell(fcache, 0) + substitution_cost<Mode>(i, query_code);
  int shift;
  if constexpr(Mode == CustomPadMode::PaddingAware) {
    shift = std::min(
      cell(*ccache, -1) + query_gap_cost<Mode>(query_code),
      cell(*ccache, 1) + target_gap_cost<Mode>(i)
    );
  } else {
    shift = std::min(cell(*ccache, -1), cell(*ccache, 1)) + gap_cost_;
  }
  set_cell(ccache, 0, std::min(mmatch, shift));
}

template <int MaxBand>
STARTREE_ALWAYS_INLINE
void CustomTrie<MaxBand>::push_hit(const NodeType& node,
                                   const int distance) {
  if(node.sequence_index != kNoSequence) {
    hits_->push_back(Hit{node.sequence_index, distance});
  }
}

template <int MaxBand>
void CustomTrie<MaxBand>::dash(uint32_t node_idx,
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
      if(ch == kNoCustomChild) {
        return;
      }
      node_idx = ch;
      offset = 0;
    }
  }
  push_hit(nodes_[static_cast<size_t>(node_idx)], threshold_);
}

template <int MaxBand>
template <CustomPadMode Mode>
STARTREE_ALWAYS_INLINE
bool CustomTrie<MaxBand>::post_step(uint32_t node_idx,
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
      CustomPebble{node_idx, offset, path}
    );
  }
  if constexpr(Mode == CustomPadMode::NoPad) {
    if(!prefix_has_no_match_[static_cast<size_t>(depth)] &&
       depth > seed_depth_) {
      bool can_dash = true;
      const int maxa = active_band(depth);
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
template <CustomPadMode Mode>
STARTREE_ALWAYS_INLINE
CustomEdgeAdvance CustomTrie<MaxBand>::advance_forced_edge_character(
    uint32_t* fnode,
    int* foffset,
    uint32_t* fpath,
    int* focus_depth) {
  const NodeType& fnode_ref = nodes_[static_cast<size_t>(*fnode)];
  const int edge_len = static_cast<int>(fnode_ref.edge.size());
  if(*foffset + 1 >= edge_len) {
    return CustomEdgeAdvance::AtBranch;
  }

  const Cache& fcache =
    (*foffset >= 0) ? fnode_ref.caches[static_cast<size_t>(*foffset)]
                    : initial_cache_;
  const int next_depth = *focus_depth + 1;
  const int next_off = *foffset + 1;
  const int code = fnode_ref.edge[static_cast<size_t>(next_off)];
  Cache& ncache = nodes_[static_cast<size_t>(*fnode)]
    .caches[static_cast<size_t>(next_off)];

  CustomEdgeAdvance next_advance = CustomEdgeAdvance::Advanced;
  bool child_touches_padding = false;
  if constexpr(Mode == CustomPadMode::NoPad) {
    const auto common = compute_child_independent_dp<CustomPadMode::NoPad>(
      fcache, *fpath, next_depth
    );
    compute_dp_update_and_cache<CustomPadMode::NoPad>(
      fcache, common, code, next_depth, &ncache
    );
  } else {
    const int maxa = active_band(next_depth);
    const bool common_touches_padding =
      common_band_touches_padding(*fpath, next_depth, maxa);
    const auto common = common_touches_padding
      ? compute_child_independent_dp<CustomPadMode::PaddingAware>(
          fcache, *fpath, next_depth
        )
      : compute_child_independent_dp<CustomPadMode::NoPad>(
          fcache, *fpath, next_depth
        );

    child_touches_padding =
      common_touches_padding || code == kPad ||
      next_depth - maxa <= query_pad_depth_;
    if(child_touches_padding) {
      compute_dp_update_and_cache<CustomPadMode::PaddingAware>(
        fcache, common, code, next_depth, &ncache
      );
    } else {
      compute_dp_update_and_cache<CustomPadMode::NoPad>(
        fcache, common, code, next_depth, &ncache
      );
      next_advance = CustomEdgeAdvance::SwitchToNoPad;
    }
  }
  if(cell(ncache, 0) > threshold_) {
    return CustomEdgeAdvance::Done;
  }

  const uint32_t npath = (*fpath << 4) | static_cast<uint32_t>(code & 0xF);
  bool keep_going;
  if constexpr(Mode == CustomPadMode::NoPad) {
    keep_going =
      post_step<CustomPadMode::NoPad>(
        *fnode, next_off, npath, next_depth, ncache
      );
  } else {
    keep_going = child_touches_padding
      ? post_step<CustomPadMode::PaddingAware>(
          *fnode, next_off, npath, next_depth, ncache
        )
      : post_step<CustomPadMode::NoPad>(
          *fnode, next_off, npath, next_depth, ncache
        );
  }
  if(!keep_going) {
    return CustomEdgeAdvance::Done;
  }

  *foffset = next_off;
  *fpath = npath;
  *focus_depth = next_depth;
  return next_advance;
}

template <int MaxBand>
template <CustomPadMode Mode>
void CustomTrie<MaxBand>::iterate_children(const uint32_t fnode,
                                           const int foffset,
                                           const uint32_t fpath,
                                           const int focus_depth) {
  const NodeType& fnode_ref = nodes_[static_cast<size_t>(fnode)];
  const Cache& fcache =
    (foffset >= 0) ? fnode_ref.caches[static_cast<size_t>(foffset)]
                   : initial_cache_;
  const int next_depth = focus_depth + 1;

  unsigned mask = fnode_ref.child_mask;
  if constexpr(Mode == CustomPadMode::NoPad) {
    const auto common = compute_child_independent_dp<CustomPadMode::NoPad>(
      fcache, fpath, next_depth
    );
    while(mask) {
      const int code = ctz_u(mask);
      mask &= mask - 1;
      const uint32_t child_idx = fnode_ref.child[static_cast<size_t>(code)];
      Cache& ncache = nodes_[static_cast<size_t>(child_idx)].caches[0];
      compute_dp_update_and_cache<CustomPadMode::NoPad>(
        fcache, common, code, next_depth, &ncache
      );
      if(cell(ncache, 0) > threshold_) {
        continue;
      }
      const uint32_t npath = (fpath << 4) | static_cast<uint32_t>(code & 0xF);
      if(!post_step<CustomPadMode::NoPad>(
          child_idx, 0, npath, next_depth, ncache)) {
        continue;
      }
      recurse<CustomPadMode::NoPad>(child_idx, 0, npath, next_depth);
    }
  } else {
    const int maxa = active_band(next_depth);
    const bool common_touches_padding =
      common_band_touches_padding(fpath, next_depth, maxa);
    const auto common = common_touches_padding
      ? compute_child_independent_dp<CustomPadMode::PaddingAware>(
          fcache, fpath, next_depth
        )
      : compute_child_independent_dp<CustomPadMode::NoPad>(
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
        compute_dp_update_and_cache<CustomPadMode::PaddingAware>(
          fcache, common, code, next_depth, &ncache
        );
      } else {
        compute_dp_update_and_cache<CustomPadMode::NoPad>(
          fcache, common, code, next_depth, &ncache
        );
      }
      if(cell(ncache, 0) > threshold_) {
        continue;
      }
      const uint32_t npath = (fpath << 4) | static_cast<uint32_t>(code & 0xF);
      const bool keep_going = child_touches_padding
        ? post_step<CustomPadMode::PaddingAware>(
            child_idx, 0, npath, next_depth, ncache
          )
        : post_step<CustomPadMode::NoPad>(
            child_idx, 0, npath, next_depth, ncache
          );
      if(!keep_going) {
        continue;
      }
      if(child_touches_padding) {
        recurse<CustomPadMode::PaddingAware>(child_idx, 0, npath, next_depth);
      } else {
        recurse<CustomPadMode::NoPad>(child_idx, 0, npath, next_depth);
      }
    }
  }
}

template <int MaxBand>
template <CustomPadMode Mode>
void CustomTrie<MaxBand>::recurse(uint32_t fnode,
                                  int foffset,
                                  uint32_t fpath,
                                  int focus_depth) {
  while(true) {
    const CustomEdgeAdvance advance =
      advance_forced_edge_character<Mode>(
        &fnode, &foffset, &fpath, &focus_depth
      );
    if(advance == CustomEdgeAdvance::Advanced) {
      continue;
    }
    if constexpr(Mode == CustomPadMode::PaddingAware) {
      if(advance == CustomEdgeAdvance::SwitchToNoPad) {
        recurse<CustomPadMode::NoPad>(fnode, foffset, fpath, focus_depth);
        return;
      }
    }
    if(advance == CustomEdgeAdvance::AtBranch) {
      iterate_children<Mode>(fnode, foffset, fpath, focus_depth);
    }
    return;
  }
}

}  // namespace radix
}  // namespace startree

#endif  // STARTREE_CUSTOM_TRIE_H
