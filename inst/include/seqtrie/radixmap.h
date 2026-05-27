#ifndef seqtrie_RADIXMAP_H
#define seqtrie_RADIXMAP_H

#include "seqtrie/utility.h"
#include "seqtrie/radix_child_map.h"
#include "seqtrie/node_pool.h"
#include "simple_array/simple_array.h"
#include "simple_array/small_array.h"
#include <algorithm>
#include <functional>
#include <type_traits>
#include <limits>
#include <stdexcept>
#include <utility>
#include <tuple>
#include <vector>

// Inline stack capacity for branch_type's small_array. The default lives here
// so the header is self-contained; override with -DSEQTRIE_SMALL_ARRAY_SIZE=<n>
// (e.g. via PKG_CPPFLAGS) before including to change it.
#ifndef SEQTRIE_SMALL_ARRAY_SIZE
#define SEQTRIE_SMALL_ARRAY_SIZE 32
#endif

namespace seqtrie {

template <typename T>
class NodePool;
class RadixForest;

class RadixMap {
public:
  // fundamental types
  using atomic_type = char;
  using branch_type = trqwe::small_array<atomic_type,
                                         std::allocator<atomic_type>,
                                         size_t,
                                         std::integral_constant<size_t, SEQTRIE_SMALL_ARRAY_SIZE>>;
  using index_type = size_t;
  using size_type = size_t;
  using self_type = RadixMap;
  using pointer_type = self_type*;
  using const_weak_pointer_type = const self_type*;
  using weak_pointer_type = self_type*;
  using map_type = seqtrie::radix_child_map<pointer_type>;
  using pairchar_type = std::pair<atomic_type, atomic_type>;
  using span_type = nonstd::span<const atomic_type>;
  using affine_lane_type = trqwe::simple_array<int>;
  using affine_col_type = std::tuple<affine_lane_type, affine_lane_type, affine_lane_type>;

  // constants
  static constexpr index_type nullidx  = std::numeric_limits<index_type>::max();
  static constexpr atomic_type GAP_CHAR      = atomic_type(0);                                 // '\0' gap for non-affine
  static constexpr atomic_type GAP_OPEN_CHAR = std::numeric_limits<atomic_type>::min();        // '\255' gap open for affine
  static constexpr atomic_type GAP_EXTN_CHAR = atomic_type(0);                                 // '\0' gap extension for affine
  static constexpr int         NO_ALIGN     = std::numeric_limits<int>::max() / 2;             // impossible positions

private:
  branch_type             branch;
  index_type              terminal_idx;
  map_type                child_nodes;
  const_weak_pointer_type parent_node;

public:
  RadixMap() : branch(), terminal_idx(nullidx), child_nodes(), parent_node(nullptr) {}

  struct path {
    const_weak_pointer_type m;
    path() : m(nullptr) {}
    path(const_weak_pointer_type x) : m(x) {}
    const_weak_pointer_type operator->() const { return m; }
    path parent() const { return m ? path(m->get_parent_node()) : path(); }
    path root() const {
      const_weak_pointer_type cur = m;
      const_weak_pointer_type parent_ptr = cur ? cur->get_parent_node() : nullptr;
      if(parent_ptr == nullptr) return path(cur);
      path p(parent_ptr);
      return p.root();
    }
    bool is_root() const { return m && (m->get_parent_node() == nullptr); }
  };

  struct BandLimits {
    size_t lower;
    size_t upper;

    BandLimits() noexcept : lower(0), upper(0) {}

    BandLimits(size_t char_depth, size_t band_radius, size_t query_len) noexcept
      : lower((char_depth > band_radius) ? char_depth - band_radius : 0),
        upper(std::min(query_len, char_depth + band_radius)) {
      if(lower > query_len) lower = query_len;
    }

    static inline size_t linear_radius(int max_distance, int gap_cost, size_t /*query_len*/) noexcept {
      return static_cast<size_t>(max_distance / gap_cost);
    }

    static inline size_t affine_radius(int max_distance, int gap_cost, int gap_open_including_first_extension, size_t /*query_len*/) noexcept {
      if(max_distance < gap_open_including_first_extension) return 0;
      const int remaining = max_distance - gap_open_including_first_extension;
      return 1 + static_cast<size_t>(remaining / gap_cost);
    }

    inline void increment(size_t previous_char_depth, size_t band_radius, size_t query_len) noexcept {
      if(previous_char_depth >= band_radius && lower < query_len) ++lower;
      if(upper < query_len) ++upper;
    }
  };

  struct search_result {
    std::vector<path> match;
    std::vector<int> distance;

    void append(const search_result & other) {
      match.insert(match.end(), other.match.begin(), other.match.end());
      distance.insert(distance.end(), other.distance.begin(), other.distance.end());
    }
  };

  struct search_context : search_result {
    span_type query;
    int max_distance;
    size_t band_radius;

    search_context() : search_result(), query(), max_distance(0), band_radius(0) {}
    search_context(span_type query_value, int max_distance_value)
      : search_result(), query(query_value), max_distance(max_distance_value),
        band_radius(static_cast<size_t>(max_distance_value)) {}

    void update_max_distance(int value) noexcept {
      max_distance = value;
      band_radius = static_cast<size_t>(value);
    }

    void update_band_radius(size_t value) noexcept {
      band_radius = value;
    }

    BandLimits band(size_t char_depth) const noexcept {
      return BandLimits(char_depth, band_radius, query.size());
    }

    search_result release_result() && {
      search_result out;
      out.match = std::move(match);
      out.distance = std::move(distance);
      return out;
    }
  };

  struct search_context_cached : search_result {
    const QueryCostCache & query_costs;
    int max_distance;
    size_t band_radius;

    search_context_cached(const QueryCostCache & query_costs_value, int max_distance_value, size_t band_radius_value)
      : search_result(), query_costs(query_costs_value), max_distance(max_distance_value),
        band_radius(band_radius_value) {}

    void update_max_distance(int value, size_t band_radius_value) noexcept {
      max_distance = value;
      band_radius = band_radius_value;
    }

    void update_band_radius(size_t value) noexcept {
      band_radius = value;
    }

    BandLimits band(size_t char_depth) const noexcept {
      return BandLimits(char_depth, band_radius, query_costs.query_len);
    }

    search_result release_result() && {
      search_result out;
      out.match = std::move(match);
      out.distance = std::move(distance);
      return out;
    }
  };

  enum class SearchHookAction : unsigned char {
    add_continue,
    add_stop,
    skip_continue,
    skip_stop
  };

  struct NullSearchHook {
    template <typename Context>
    SearchHookAction operator()(Context &, path, int) const noexcept {
      return SearchHookAction::add_continue;
    }
  };

  template <typename Hook>
  using hook_ret = std::conditional_t<std::is_same_v<Hook, NullSearchHook>, void, bool>;

  // Each per-depth column is a full DP column of
  // size query_len+1, banded inner loop, NO_ALIGN sentinels outside the band. Storage
  // is vector<simple_array<int>>; simple_array's stored size acts as capacity (grows but
  // never shrinks across queries); column_len tracks the used prefix = query_len + 1.
  struct UnitFullSimpleArrayWorkspace {
    using column_type = trqwe::simple_array<int>;
    std::vector<column_type> columns;
    size_t column_len = 0; // used prefix = query_len + 1

    void reset_global(size_t query_len, size_t depth_hint = 0) {
      column_len = query_len + 1;
      if(depth_hint > 0 && columns.capacity() < depth_hint) columns.reserve(depth_hint);
      reset_global_root(query_len);
    }

    inline void ensure_capacity(column_type & slot) const {
      if(slot.size() < column_len) {
        size_t new_cap = slot.size() * 2;
        if(new_cap < column_len) new_cap = column_len;
        slot.reset(new_cap); // grow-only realloc, no copy needed (overwritten below)
      }
    }

    void reset_global_root(size_t /*query_len*/) {
      if(columns.empty()) columns.emplace_back();
      ensure_capacity(columns[0]);
      int * root_col = columns[0].data();
      for(size_t i = 0; i < column_len; ++i) {
        root_col[i] = static_cast<int>(i);
      }
    }

    void ensure_child_slot(size_t depth) {
      const size_t required = depth + 2;
      if(columns.size() < required) columns.resize(required);
    }

    int *       at(size_t depth)       { return columns[depth].data(); }
    const int * at(size_t depth) const { return columns[depth].data(); }

    int * clone_from_parent(size_t depth) {
      ensure_child_slot(depth);
      int * parent = columns[depth].data();
      column_type & child_slot = columns[depth + 1];
      ensure_capacity(child_slot);
      int * child = child_slot.data();
      std::copy(parent, parent + column_len, child);
      return child;
    }
  };

  struct AffineWorkspace {
    std::vector<affine_col_type> columns;
    size_t column_len = 0; // used prefix = query_len + 1

    void initialize(size_t query_len, size_t depth_hint = 0) {
      column_len = query_len + 1;
      columns.clear();
      if(depth_hint > 0 && columns.capacity() < depth_hint) columns.reserve(depth_hint);
      columns.emplace_back();
      ensure_capacity(columns[0]);
    }

    inline void ensure_capacity(affine_lane_type & slot) const {
      if(slot.size() < column_len) {
        size_t new_cap = slot.size() * 2;
        if(new_cap < column_len) new_cap = column_len;
        slot.reset(new_cap); // grow-only realloc, no copy needed (overwritten below)
      }
    }

    inline void ensure_capacity(affine_col_type & col) const {
      ensure_capacity(std::get<0>(col));
      ensure_capacity(std::get<1>(col));
      ensure_capacity(std::get<2>(col));
    }

    void ensure_child_slot(size_t depth) {
      const size_t required_size = depth + 2;
      if(columns.size() < required_size) columns.resize(required_size);
    }

    affine_col_type & at(const size_t depth) { return columns[depth]; }
    const affine_col_type & at(const size_t depth) const { return columns[depth]; }

    affine_col_type & clone_from_parent(size_t depth) {
      ensure_child_slot(depth);
      const affine_col_type & parent_col = columns[depth];
      affine_col_type & child_col = columns[depth + 1];
      ensure_capacity(child_col);
      std::copy(std::get<0>(parent_col).data(), std::get<0>(parent_col).data() + column_len, std::get<0>(child_col).data());
      std::copy(std::get<1>(parent_col).data(), std::get<1>(parent_col).data() + column_len, std::get<1>(child_col).data());
      std::copy(std::get<2>(parent_col).data(), std::get<2>(parent_col).data() + column_len, std::get<2>(child_col).data());
      return child_col;
    }
  };

  struct SingleGapCol {
    int lower;
    int diag;
    int upper;
  };
  struct SingleGapWorkspace {
    std::vector<SingleGapCol> columns;

    void initialize(SingleGapCol initial, size_t depth_hint = 0) {
      columns.clear();
      if(depth_hint > 0) {
        if(columns.capacity() < depth_hint) columns.reserve(depth_hint);
      }
      columns.emplace_back(initial);
    }

    void ensure_child_slot(size_t depth) {
      size_t required_size = depth + 2;
      if(columns.size() < required_size) {
        columns.resize(required_size);
      }
    }

    SingleGapCol & at(const size_t depth) { return columns[depth]; }
    const SingleGapCol & at(const size_t depth) const { return columns[depth]; }

    SingleGapCol & clone_from_parent(size_t depth) {
      ensure_child_slot(depth);
      const auto & parent_col = columns[depth];
      SingleGapCol & child_col = columns[depth + 1];
      child_col = parent_col;
      return child_col;
    }
  };

  const map_type & get_child_nodes() const { return child_nodes; }
  const branch_type & get_branch() const { return branch; }
  const_weak_pointer_type get_parent_node() const { return parent_node; }
  index_type get_terminal_idx() const { return terminal_idx; }

  template <typename ST>
  ST sequence() const;

  size_type size() const;
  bool validate(const bool is_root = true) const;
  std::vector<path> all(size_t max_depth = std::numeric_limits<size_t>::max()) const;
  std::string print() const { return print_impl(0); }
  std::pair<std::vector<path>, std::vector<path>> graph(size_t max_depth = std::numeric_limits<size_t>::max()) const;
  path find(const span_type query) const;
  path insert(const span_type sequence, index_type idx);
  path insert(const span_type sequence, index_type idx, NodePool<self_type> & pool);
  // Ensure sequence exists and return its terminal path
  path insert_get_path(const span_type sequence, index_type idx);
  path insert_get_path(const span_type sequence, index_type idx, NodePool<self_type> & pool);
  index_type erase(const span_type sequence);
  std::vector<path> prefix_search(const span_type query) const;

  // Rebuild pooled descendants into a fresh pool in the same child iteration
  // order used by the search implementations. The root object itself is kept
  // in place so external ownership of the root stays stable.
  void compact_search_order(NodePool<self_type> & new_pool);

  // search entry points
  search_context hamming_search(const span_type query, int max_distance) const;
  template <typename Hook>
  search_context hamming_search(const span_type query, int max_distance, Hook & hook) const;
  search_context global_search(const span_type query, int max_distance) const;
  template <typename Hook>
  search_context global_search(const span_type query, int max_distance, Hook & hook) const;
  search_context anchored_search(const span_type query, int max_distance) const;
  template <typename Hook>
  search_context anchored_search(const span_type query, int max_distance, Hook & hook) const;

  search_result global_search_linear(const span_type query, int max_distance, const CostMap & cost_map) const;
  template <typename Hook>
  search_result global_search_linear(const span_type query, int max_distance, const CostMap & cost_map, Hook & hook) const;
  search_result anchored_search_linear(const span_type query, int max_distance, const CostMap & cost_map) const;
  template <typename Hook>
  search_result anchored_search_linear(const span_type query, int max_distance, const CostMap & cost_map, Hook & hook) const;

  search_result global_search_affine(const span_type query, int max_distance, const CostMap & cost_map) const;
  template <typename Hook>
  search_result global_search_affine(const span_type query, int max_distance, const CostMap & cost_map, Hook & hook) const;
  search_result anchored_search_affine(const span_type query, int max_distance, const CostMap & cost_map) const;
  template <typename Hook>
  search_result anchored_search_affine(const span_type query, int max_distance, const CostMap & cost_map, Hook & hook) const;
  search_context single_gap_search(const span_type query, int max_distance, const int gap_cost) const;
  template <typename Hook>
  search_context single_gap_search(const span_type query, int max_distance, const int gap_cost, Hook & hook) const;

private:
  // Implementation functions
  template <bool ReturnPathAlways>
  path insert_impl(const span_type sequence, index_type idx);
  template <bool ReturnPathAlways>
  path insert_impl(const span_type sequence, index_type idx, NodePool<self_type> & pool);

  friend class RadixForest;

  // search helpers
  template <typename>
  struct dependent_false : std::false_type {};

  static SearchHookAction normalize_hook_result(SearchHookAction action) noexcept { return action; }
  static bool hook_action_adds(SearchHookAction action) noexcept {
    return action == SearchHookAction::add_continue || action == SearchHookAction::add_stop;
  }
  static bool hook_action_continues(SearchHookAction action) noexcept {
    return action == SearchHookAction::add_continue || action == SearchHookAction::skip_continue;
  }

  template <typename Context, typename Hook>
  static SearchHookAction call_search_hook(Context & ctx, path node_path, int distance, Hook & hook);

  template <typename Context, typename Hook>
  static hook_ret<Hook> search_add(Context & ctx, path node_path, int distance, Hook & hook);

  template <typename Context, typename Hook>
  static hook_ret<Hook> search_add_all(Context & ctx, path node_path, int distance, Hook & hook);

  static int min3(const int a, const int b, const int c) noexcept {
    int out = a < b ? a : b;
    return c < out ? c : out;
  }

  static pointer_type compact_subtree_search_order(weak_pointer_type source,
                                                   const_weak_pointer_type new_parent,
                                                   NodePool<self_type> & new_pool);

  static void mark_removed_subtree(weak_pointer_type node);

  std::string print_impl(size_t depth) const;
  enum class erase_action { erase, merge, keep };
  static erase_action erase_impl(weak_pointer_type node, const span_type sequence, index_type & result);


  // search implementations
  template <typename Hook>
  static hook_ret<Hook> hamming_search_impl(const_weak_pointer_type node, size_t position, int distance, search_context & ctx, Hook & hook);
  static int update_col_banded(atomic_type branchval, const span_type query, std::vector<int> & col, size_t lower, size_t upper);
  // Same logic as update_col_banded but operates on a raw int* of length col_len.
  static int update_col_banded_ptr(atomic_type branchval, const span_type query, int * __restrict__ col, size_t col_len, size_t lower, size_t upper);
  template <typename Hook>
  static hook_ret<Hook> global_search_impl(const_weak_pointer_type node, size_t node_depth, size_t char_depth, search_context & ctx, UnitFullSimpleArrayWorkspace & workspace, Hook & hook, int col_min); // unitary cost
  template <typename Hook>
  static hook_ret<Hook> anchored_search_impl(const_weak_pointer_type node, size_t node_depth, size_t char_depth, int row_min, search_context & ctx, UnitFullSimpleArrayWorkspace & workspace, Hook & hook, int col_min); // unitary cost

  static int update_col_linear_banded(atomic_type branchval, std::vector<int> & col, size_t lower, size_t upper, const QueryCostCache & query_costs);
  static int update_col_linear_banded_ptr(atomic_type branchval, int * __restrict__ col, size_t col_len, size_t lower, size_t upper, const QueryCostCache & query_costs);
  template <typename Hook>
  static hook_ret<Hook> global_search_linear_impl(const_weak_pointer_type node, size_t node_depth, size_t char_depth, search_context_cached & ctx, UnitFullSimpleArrayWorkspace & workspace, Hook & hook, int col_min);
  template <typename Hook>
  static hook_ret<Hook> anchored_search_linear_impl(const_weak_pointer_type node, size_t node_depth, size_t char_depth, int row_min, search_context_cached & ctx, UnitFullSimpleArrayWorkspace & workspace, Hook & hook, int col_min);

  static int update_col_affine_banded(atomic_type branchval, affine_col_type & col, size_t lower, size_t upper, const QueryCostCache & query_costs);
  template <typename Hook>
  static hook_ret<Hook> global_search_affine_impl(const_weak_pointer_type node, size_t node_depth, size_t char_depth, search_context_cached & ctx, AffineWorkspace & workspace, Hook & hook, int col_min);
  template <typename Hook>
  static hook_ret<Hook> anchored_search_affine_impl(const_weak_pointer_type node, size_t node_depth, size_t char_depth, int row_min, search_context_cached & ctx, AffineWorkspace & workspace, Hook & hook, int col_min);

  static int update_col_single_gap(atomic_type branchval, const span_type query, size_t ref_len, SingleGapCol & col, const int gap_cost);
  template <typename Hook>
  static hook_ret<Hook> single_gap_search_impl(const_weak_pointer_type node, size_t node_depth, size_t char_depth, int row_min, search_context & ctx, SingleGapWorkspace & workspace, const int gap_cost, Hook & hook, int col_min);
};

// implementations


template <typename ST>
inline ST RadixMap::sequence() const {
  static_assert(std::is_same<typename ST::value_type, atomic_type>::value, "Output sequence value_type must be the same as RadixMap::atomic_type.");
  const_weak_pointer_type current = this;
  std::vector<const_weak_pointer_type> h;
  size_t size = 0;
  while(current != nullptr) {
    h.push_back(current);
    size += current->branch.size();
    current = current->parent_node;
  }
  ST result = array_allocate<ST>(size);
  atomic_type * result_ptr = array_data(result);
  for(auto it = h.rbegin(); it != h.rend(); ++it) {
    current = *it;
    result_ptr = std::copy(current->branch.begin(), current->branch.end(), result_ptr);
  }
  return result;
}

inline std::vector<RadixMap::path> RadixMap::all(size_t max_depth) const {
  std::vector<path> result;
  if(terminal_idx != nullidx) {
    result.push_back(this);
  }
  if(max_depth == 0) return result;
  for(auto ch : child_nodes) {
    size_t next_depth = (max_depth == std::numeric_limits<size_t>::max()) ? max_depth : (max_depth - 1);
    std::vector<path> x = ch.second->all(next_depth);
    appendspan(result, x);
  }
  return result;
}

inline RadixMap::size_type RadixMap::size() const {
  size_type result = terminal_idx == nullidx ? 0 : 1;
  for(auto ch : child_nodes) { result += ch.second->size(); }
  return result;
}

inline bool RadixMap::validate(const bool is_root) const {
  // checks:
  // 1) all parent_node values are correct
  // root can be distinguished by parent_node == nullptr
  // 2) map key matches value branch
  // 3) A nullidx node cannot have only <=1 child otherwise they should've been combined, except root
  if(is_root) {
    if(parent_node != nullptr) return false;
  } else {
    if(parent_node == nullptr) return false;
    if(terminal_idx == nullidx) {
      if(child_nodes.size() <= 1) return false;
    }
  }
  for(auto ch : child_nodes) {
    if(ch.first != ch.second->branch[0]) return false;
    if(ch.second->parent_node != this) return false;
    if(!ch.second->validate(false)) return false;
  }
  return true;
}

inline std::string RadixMap::print_impl(size_t depth) const {
  std::string result;
  if(depth == 0) {
    result += "(root)";
  } else {
    for(size_t i=0; i<depth; ++i) result += " ";
  }
  std::string x(branch.size(), 0);
  for(size_t j=0; j<branch.size(); ++j) {
    x[j] = static_cast<char>(branch[j]);
  }
  result += x;
  
  if(terminal_idx != nullidx) {
    result += ":";
    result += std::to_string(terminal_idx);
  }
  // result += ",";
  // result += ptr_tostring(this);
  result += "\n";
  std::vector<atomic_type> child_node_keys;
  for(auto ch : child_nodes) {
    child_node_keys.push_back(ch.first);
  }
  std::sort(child_node_keys.begin(), child_node_keys.end()); // sort for reproducible printing; unnecessary for std::map, but doesn't matter much
  for(atomic_type k : child_node_keys) {
    result += child_nodes.at(k)->print_impl(depth + 1);
  }
  if(depth == 0) {
    result += "\n";
  }
  return result;
}

inline std::pair<std::vector<RadixMap::path>, std::vector<RadixMap::path>> RadixMap::graph(size_t max_depth) const {
  std::pair<std::vector<path>, std::vector<path>> result;
  if(parent_node != nullptr) {
    result.first.push_back(path(parent_node));
    result.second.push_back(path(this));
  }
  if(max_depth == 0) return result;
  for(auto ch : child_nodes) {
    size_t next_depth = (max_depth == std::numeric_limits<size_t>::max()) ? max_depth : (max_depth - 1);
    auto x = ch.second->graph(next_depth);
    appendspan(result.first, x.first);
    appendspan(result.second, x.second);
  }
  return result;
}

inline RadixMap::path RadixMap::find(const RadixMap::span_type query) const {
  const_weak_pointer_type node = this;
  size_t position=0;
  while(position < query.size()) {
    auto it = node->child_nodes.find(query[position]);
    if(it == node->child_nodes.end()) return path();
    node = it->second;
    if(position + node->branch.size() > query.size()) return path();
    for(size_t j=0; j<node->branch.size(); ++j) {
      if(node->branch[j] != query[position+j]) return path();
    }
    position += node->branch.size();
  }
  return path(node);
}

inline RadixMap::path RadixMap::insert(const RadixMap::span_type, index_type) {
  throw std::runtime_error("RadixMap::insert requires an external NodePool");
}

inline RadixMap::path RadixMap::insert(const RadixMap::span_type sequence, index_type idx, NodePool<self_type> & pool) {
  return insert_impl<false>(sequence, idx, pool);
}

inline RadixMap::path RadixMap::insert_get_path(const RadixMap::span_type, index_type) {
  throw std::runtime_error("RadixMap::insert_get_path requires an external NodePool");
}

inline RadixMap::path RadixMap::insert_get_path(const RadixMap::span_type sequence, index_type idx, NodePool<self_type> & pool) {
  return insert_impl<true>(sequence, idx, pool);
}

template <bool ReturnPathAlways>
inline RadixMap::path RadixMap::insert_impl(const RadixMap::span_type, index_type) {
  throw std::runtime_error("RadixMap::insert_impl requires an external NodePool");
}

template <bool ReturnPathAlways>
inline RadixMap::path RadixMap::insert_impl(const RadixMap::span_type sequence, index_type idx, NodePool<self_type> & pool) {
  if(sequence.size() == 0) {
    if(terminal_idx == nullidx) {
      terminal_idx = idx;
      if constexpr (ReturnPathAlways) {
        return path(this);
      } else {
        return path();
      }
    } else {
      return path(this);
    }
  }
  atomic_type s = sequence[0];
  auto it = child_nodes.find(s);
  if(it == child_nodes.end()) {
    auto * new_child = pool.create();
    child_nodes.emplace(s, new_child);
    new_child->parent_node = this;
    new_child->branch = subvector<branch_type>(sequence, 0);
    new_child->terminal_idx = idx;
    if constexpr (ReturnPathAlways) {
      return path(new_child);
    } else {
      return path();
    }
  }
  pointer_type & child_ptr = it->second;
  const size_t child_branch_size = child_ptr->branch.size();
  size_t i = 0;
  while(i < child_branch_size && i < sequence.size() && sequence[i] == child_ptr->branch[i]) { ++i; }

  if(i == sequence.size() && i == child_branch_size) {
    if(child_ptr->terminal_idx == nullidx) {
      child_ptr->terminal_idx = idx;
      if constexpr (ReturnPathAlways) {
        return path(child_ptr);
      } else {
        return path();
      }
    } else {
      return path(child_ptr);
    }
  } else if(i == sequence.size()) {
    branch_type branch_prefix = subvector<branch_type>(child_ptr->branch, 0, i);
    branch_type branch_suffix = subvector<branch_type>(child_ptr->branch, i);
    atomic_type s_insert = branch_suffix[0];
    auto * inserted_node = pool.create();
    inserted_node->parent_node = this;
    inserted_node->child_nodes[s_insert] = child_ptr;
    inserted_node->child_nodes[s_insert]->parent_node = inserted_node;
    inserted_node->branch = std::move(branch_prefix);
    inserted_node->terminal_idx = idx;
    inserted_node->child_nodes[s_insert]->branch = std::move(branch_suffix);
    child_ptr = inserted_node;
    if constexpr (ReturnPathAlways) {
      return path(child_ptr);
    } else {
      return path();
    }
  } else if(i == child_branch_size) {
    span_type seq_suffix = sequence.subspan(i);
    return child_ptr->insert_impl<ReturnPathAlways>(seq_suffix, idx, pool);
  } else {
    branch_type branch_prefix = subvector<branch_type>(child_ptr->branch, 0, i);
    branch_type branch_suffix = subvector<branch_type>(child_ptr->branch, i);
    branch_type seq_suffix = subvector<branch_type>(sequence, i);
    atomic_type s_insert_branch = branch_suffix[0];
    atomic_type s_insert_seq = seq_suffix[0];

    auto * inserted_node = pool.create();
    inserted_node->parent_node = this;
    inserted_node->child_nodes[s_insert_branch] = child_ptr;
    inserted_node->child_nodes[s_insert_branch]->parent_node = inserted_node;
    inserted_node->child_nodes[s_insert_branch]->branch = std::move(branch_suffix);
    auto * new_seq_child = pool.create();
    inserted_node->child_nodes[s_insert_seq] = new_seq_child;
    new_seq_child->parent_node = inserted_node;
    new_seq_child->branch = std::move(seq_suffix);
    new_seq_child->terminal_idx = idx;
    child_ptr = inserted_node;
    child_ptr->branch = std::move(branch_prefix);
    if constexpr (ReturnPathAlways) {
      return path(new_seq_child);
    } else {
      return path();
    }
  }
}


inline RadixMap::index_type RadixMap::erase(const RadixMap::span_type sequence) {
  index_type result = nullidx;
  erase_impl(this, sequence, result);
  return result;
}

inline std::vector<RadixMap::path> RadixMap::prefix_search(const RadixMap::span_type query) const {
  const_weak_pointer_type node = this;
  size_t query_position = 0;
  size_t branch_position = 0;
  while(query_position < query.size()) {
    if(branch_position >= node->branch.size()) {
      auto it = node->child_nodes.find(query[query_position]);
      if(it == node->child_nodes.end()) return std::vector<path>();
      node = it->second;
      branch_position = 0;
    }
    if(node->branch[branch_position] == query[query_position]) {
      branch_position++;
      query_position++;
    } else {
      return std::vector<path>();
    }
  }
  return node->all();
}

inline RadixMap::pointer_type
RadixMap::compact_subtree_search_order(RadixMap::weak_pointer_type source,
                                       const RadixMap::const_weak_pointer_type new_parent,
                                       NodePool<self_type> & new_pool) {
  pointer_type copy = new_pool.create();
  copy->terminal_idx = source->terminal_idx;
  copy->parent_node = new_parent;

  copy->child_nodes.reserve_exact(source->child_nodes.size());
  for(auto ch : source->child_nodes) {
    pointer_type child_copy = compact_subtree_search_order(ch.second, copy, new_pool);
    copy->child_nodes.append_unchecked(ch.first, child_copy);
  }

  copy->branch = std::move(source->branch);

  return copy;
}

inline void RadixMap::compact_search_order(NodePool<self_type> & new_pool) {
  map_type new_children;
  new_children.reserve_exact(child_nodes.size());

  for(auto ch : child_nodes) {
    pointer_type child_copy = compact_subtree_search_order(ch.second, this, new_pool);
    new_children.append_unchecked(ch.first, child_copy);
  }

  child_nodes = std::move(new_children);
}

inline RadixMap::search_context RadixMap::hamming_search(const span_type query, const int max_distance) const {
  NullSearchHook hook;
  return hamming_search(query, max_distance, hook);
}

template <typename Hook>
inline RadixMap::search_context RadixMap::hamming_search(const span_type query, const int max_distance, Hook & hook) const {
  search_context ctx(query, max_distance);
  hamming_search_impl<Hook>(this, 0, 0, ctx, hook);
  return ctx;
}

inline RadixMap::search_context RadixMap::global_search(const RadixMap::span_type query, const int max_distance) const {
  NullSearchHook hook;
  return global_search(query, max_distance, hook);
}

template <typename Hook>
inline RadixMap::search_context RadixMap::global_search(const RadixMap::span_type query, const int max_distance, Hook & hook) const {
  search_context ctx(query, max_distance);
  UnitFullSimpleArrayWorkspace workspace;
  workspace.reset_global(query.size(), query.size() + ctx.band_radius + 1);
  global_search_impl<Hook>(this, 0, 0, ctx, workspace, hook, 0);
  return ctx;
}

// an "anchored" search can end on the last column or col of the dynamic programming array
// unlike global which must end on the bottom right corner
// we need to keep track of the minimum value in the last row
inline RadixMap::search_context RadixMap::anchored_search(const RadixMap::span_type query, const int max_distance) const {
  NullSearchHook hook;
  return anchored_search(query, max_distance, hook);
}

template <typename Hook>
inline RadixMap::search_context RadixMap::anchored_search(const RadixMap::span_type query, const int max_distance, Hook & hook) const {
  search_context ctx(query, max_distance);
  UnitFullSimpleArrayWorkspace workspace;
  workspace.reset_global(query.size(), query.size() + 1);
  anchored_search_impl<Hook>(this, 0, 0, static_cast<int>(query.size()), ctx, workspace, hook, 0);
  return ctx;
}

inline RadixMap::search_result RadixMap::global_search_linear(const RadixMap::span_type query, const int max_distance, const CostMap & cost_map) const {
  NullSearchHook hook;
  return global_search_linear(query, max_distance, cost_map, hook);
}

template <typename Hook>
inline RadixMap::search_result RadixMap::global_search_linear(const RadixMap::span_type query, const int max_distance, const CostMap & cost_map, Hook & hook) const {
  QueryCostCache query_costs(cost_map, query);
  const size_t band_radius = BandLimits::linear_radius(max_distance, query_costs.gap_cost, query_costs.query_len);
  search_context_cached ctx(query_costs, max_distance, band_radius);
  UnitFullSimpleArrayWorkspace workspace;
  workspace.reset_global(query_costs.query_len, query_costs.query_len + 1);
  int * root_col = workspace.at(0);
  for(size_t i = 1; i < workspace.column_len; ++i) {
    root_col[i] = root_col[i - 1] + query_costs.gap_cost; // gap in target
  }
  global_search_linear_impl<Hook>(this, 0, 0, ctx, workspace, hook, 0);
  return std::move(ctx).release_result();
}

inline RadixMap::search_result RadixMap::anchored_search_linear(const RadixMap::span_type query, const int max_distance, const CostMap & cost_map) const {
  NullSearchHook hook;
  return anchored_search_linear(query, max_distance, cost_map, hook);
}

template <typename Hook>
inline RadixMap::search_result RadixMap::anchored_search_linear(const RadixMap::span_type query, const int max_distance, const CostMap & cost_map, Hook & hook) const {
  QueryCostCache query_costs(cost_map, query);
  const size_t band_radius = BandLimits::linear_radius(max_distance, query_costs.gap_cost, query_costs.query_len);
  search_context_cached ctx(query_costs, max_distance, band_radius);
  UnitFullSimpleArrayWorkspace workspace;
  workspace.reset_global(query_costs.query_len, query_costs.query_len + 1);
  int * root_col = workspace.at(0);
  for(size_t i = 1; i < workspace.column_len; ++i) {
    root_col[i] = root_col[i - 1] + query_costs.gap_cost; // gap in target
  }
  int row_min = root_col[workspace.column_len - 1];
  anchored_search_linear_impl<Hook>(this, 0, 0, row_min, ctx, workspace, hook, 0);
  return std::move(ctx).release_result();
}

inline RadixMap::search_result RadixMap::global_search_affine(const RadixMap::span_type query, const int max_distance, const CostMap & cost_map) const {
  NullSearchHook hook;
  return global_search_affine(query, max_distance, cost_map, hook);
}

template <typename Hook>
inline RadixMap::search_result RadixMap::global_search_affine(const RadixMap::span_type query, const int max_distance, const CostMap & cost_map, Hook & hook) const {
  QueryCostCache query_costs(cost_map, query);
  const size_t band_radius = BandLimits::affine_radius(max_distance, query_costs.gap_cost, query_costs.gap_open_including_first_extension, query_costs.query_len);
  search_context_cached ctx(query_costs, max_distance, band_radius);
  const size_t col_size = query_costs.query_len + 1;
  AffineWorkspace workspace;
  workspace.initialize(query_costs.query_len, query_costs.query_len + 1);
  affine_col_type & col = workspace.at(0);
  auto & M_col =  std::get<0>(col);
  auto & X_col = std::get<1>(col);
  auto & Y_col = std::get<2>(col);
  M_col[0] = 0;
  X_col[0] = NO_ALIGN;
  Y_col[0] = NO_ALIGN;
  for(size_t i = 1; i < col_size; ++i) {
    M_col[i] = NO_ALIGN;
    X_col[i] = NO_ALIGN;
    if(i == 1) Y_col[i] = query_costs.gap_open_including_first_extension;
    else       Y_col[i] = Y_col[i - 1] + query_costs.gap_cost;
  }
  global_search_affine_impl<Hook>(this, 0, 0, ctx, workspace, hook, 0);
  return std::move(ctx).release_result();
}

inline RadixMap::search_result RadixMap::anchored_search_affine(const RadixMap::span_type query, const int max_distance, const CostMap & cost_map) const {
  NullSearchHook hook;
  return anchored_search_affine(query, max_distance, cost_map, hook);
}

template <typename Hook>
inline RadixMap::search_result RadixMap::anchored_search_affine(const RadixMap::span_type query, const int max_distance, const CostMap & cost_map, Hook & hook) const {
  QueryCostCache query_costs(cost_map, query);
  const size_t band_radius = BandLimits::affine_radius(max_distance, query_costs.gap_cost, query_costs.gap_open_including_first_extension, query_costs.query_len);
  search_context_cached ctx(query_costs, max_distance, band_radius);
  const size_t col_size = query_costs.query_len + 1;
  AffineWorkspace workspace;
  workspace.initialize(query_costs.query_len, query_costs.query_len + 1);
  affine_col_type & col = workspace.at(0);
  auto & M_col =  std::get<0>(col);
  auto & X_col = std::get<1>(col);
  auto & Y_col = std::get<2>(col);
  M_col[0] = 0;
  X_col[0] = NO_ALIGN;
  Y_col[0] = NO_ALIGN;
  for(size_t i = 1; i < col_size; ++i) {
    M_col[i] = NO_ALIGN;
    X_col[i] = NO_ALIGN;
    if(i == 1) Y_col[i] = query_costs.gap_open_including_first_extension;
    else       Y_col[i] = Y_col[i - 1] + query_costs.gap_cost;
  }
  // X_col[col_size - 1] is NO_ALIGN here (initial column), so it is omitted from the min.
  const int row_min = std::min(M_col[col_size - 1], Y_col[col_size - 1]);
  anchored_search_affine_impl<Hook>(this, 0, 0, row_min, ctx, workspace, hook, 0);
  return std::move(ctx).release_result();
}

inline RadixMap::search_context RadixMap::single_gap_search(const RadixMap::span_type query,
                                                            const int max_distance,
                                                            const int gap_cost) const {
  NullSearchHook hook;
  return single_gap_search(query, max_distance, gap_cost, hook);
}

template <typename Hook>
inline RadixMap::search_context RadixMap::single_gap_search(const RadixMap::span_type query,
                                                            const int max_distance,
                                                            const int gap_cost,
                                                            Hook & hook) const {
  search_context ctx(query, max_distance);
  SingleGapWorkspace workspace;
  SingleGapCol initial{NO_ALIGN, 0, (query.size() > 0 ? gap_cost : 0)}; // if query is empty, lower stores the row min which is 0
  workspace.initialize(initial, query.size() + 1);
  int row_min = query.size() == 0 ? 0 :
                query.size() == 1 ? gap_cost :
                NO_ALIGN;
  single_gap_search_impl<Hook>(this, 0, 0, row_min, ctx, workspace, gap_cost, hook, 0);
  return ctx;
}

template <typename Context, typename Hook>
inline RadixMap::SearchHookAction RadixMap::call_search_hook(Context & ctx,
                                                             RadixMap::path node_path,
                                                             int distance,
                                                             Hook & hook) {
  if constexpr (std::is_invocable_r_v<SearchHookAction, Hook &, Context &, path, int>) {
    return normalize_hook_result(hook(ctx, node_path, distance));
  } else {
    static_assert(dependent_false<Hook>::value,
                  "Search hook must be callable as SearchHookAction hook(ctx, path, distance)");
  }
}

template <typename Context, typename Hook>
inline RadixMap::hook_ret<Hook> RadixMap::search_add(Context & ctx,
                                                     RadixMap::path node_path,
                                                     int distance,
                                                     Hook & hook) {
  if constexpr (!std::is_same_v<Hook, NullSearchHook>) {
    const SearchHookAction action = call_search_hook(ctx, node_path, distance, hook);
    if(hook_action_adds(action)) {
      ctx.match.push_back(node_path);
      ctx.distance.push_back(distance);
    }
    return hook_action_continues(action);
  } else {
    ctx.match.push_back(node_path);
    ctx.distance.push_back(distance);
  }
}

template <typename Context, typename Hook>
inline RadixMap::hook_ret<Hook> RadixMap::search_add_all(Context & ctx,
                                                        RadixMap::path node_path,
                                                        int distance,
                                                        Hook & hook) {
  constexpr bool use_hook = !std::is_same_v<Hook, NullSearchHook>;
  const_weak_pointer_type node = node_path.m;
  if(node->terminal_idx != nullidx) {
    if constexpr (use_hook) {
      const SearchHookAction action = call_search_hook(ctx, node_path, distance, hook);
      if(hook_action_adds(action)) {
        ctx.match.push_back(node_path);
        ctx.distance.push_back(distance);
      }
      if(!hook_action_continues(action)) return false;
    } else {
      ctx.match.push_back(node_path);
      ctx.distance.push_back(distance);
    }
  }
  for(const auto & ch : node->child_nodes) {
    path child_path(ch.second);
    if constexpr (use_hook) {
      if(!search_add_all(ctx, child_path, distance, hook)) return false;
    } else {
      search_add_all(ctx, child_path, distance, hook);
    }
  }
  if constexpr (use_hook) return true;
}

inline int RadixMap::update_col_single_gap(const RadixMap::atomic_type branchval,
                                           const RadixMap::span_type query,
                                           const size_t char_depth,
                                           RadixMap::SingleGapCol & col,
                                           const int gap_cost) {
  const RadixMap::SingleGapCol prev = col;
  auto match_score = [branchval, &query](size_t pos) -> int {
    return branchval == query[pos] ? 0 : 1;
  };
  if(char_depth <= 1) {
    col.lower = gap_cost;
  } else if(char_depth <= query.size() + 1) {
    col.lower = std::min(prev.lower + match_score(char_depth - 2), prev.diag + gap_cost);
  } else {
    col.lower = NO_ALIGN;
  }
  // if diag exceeds query size, it is NO_ALIGN
  if(char_depth <= query.size()) {
    col.diag = prev.diag + match_score(char_depth - 1);
  } else {
    col.diag = NO_ALIGN;
  }
  // if upper exceeds query size, it is NO_ALIGN
  if(char_depth < query.size()) {
    col.upper = std::min(prev.upper + match_score(char_depth), col.diag + gap_cost); // depends on current col.diag, must come after col.diag update
  } else {
    col.upper = NO_ALIGN;
  }
  return min3(col.diag, col.upper, col.lower);
}

template <typename Hook>
inline RadixMap::hook_ret<Hook> RadixMap::single_gap_search_impl(RadixMap::const_weak_pointer_type node,
                                                                const size_t node_depth,
                                                                const size_t char_depth,
                                                                int row_min,
                                                                RadixMap::search_context & ctx,
                                                                SingleGapWorkspace & workspace,
                                                                const int gap_cost,
                                                                Hook & hook,
                                                                int current_col_min) {
  constexpr bool use_hook = !std::is_same_v<Hook, NullSearchHook>;
  const size_t query_len = ctx.query.size();

  SingleGapCol & current_col = workspace.at(node_depth);
  auto get_row_min = [query_len](const SingleGapCol & col, size_t row_char_depth, int current_row_min_value) -> int {
    if(query_len > 0 && row_char_depth + 1 == query_len) {
      return col.upper;
    } else if(row_char_depth == query_len) {
      return std::min(current_row_min_value, col.diag);
    } else if(row_char_depth == query_len + 1) {
      return std::min(current_row_min_value, col.lower);
    } else {
      return NO_ALIGN;
    }
  };
  int current_row_min = get_row_min(current_col, char_depth, row_min);
  if( (current_col_min > ctx.max_distance) && (current_row_min > ctx.max_distance) ) { // case 1
    if constexpr (use_hook) return true; else return;
  } else if( (current_row_min <= ctx.max_distance) && (current_row_min <= current_col_min) ) { // case 2
    if constexpr (use_hook) {
      if(!search_add_all(ctx, path(node), current_row_min, hook)) return false;
      return true;
    } else {
      search_add_all(ctx, path(node), current_row_min, hook);
      return;
    }
  } else if(node->terminal_idx != nullidx) { // case 3
    if(current_col_min <= ctx.max_distance) {
      if constexpr (use_hook) {
        if(!search_add(ctx, path(node), current_col_min, hook)) return false;
      } else {
        search_add(ctx, path(node), current_col_min, hook);
      }
    }
  }
  const int parent_min = current_col_min;
  for(const auto& ch : node->child_nodes) {
    SingleGapCol & child_col = workspace.clone_from_parent(node_depth);
    current_col_min = parent_min;
    const branch_type & branch = ch.second->branch;
    bool max_distance_exceeded = false;
    int child_row_min = row_min;
    size_t child_char_depth = char_depth;
    for(atomic_type branchval : branch) {
      ++child_char_depth;
      current_col_min = update_col_single_gap(branchval, ctx.query, child_char_depth, child_col, gap_cost);
      child_row_min = get_row_min(child_col, child_char_depth, child_row_min);
      if( (current_col_min > ctx.max_distance) && (child_row_min > ctx.max_distance) ) { // case 1
        max_distance_exceeded = true;
        break;
      } else if( (child_row_min <= ctx.max_distance) && (child_row_min <= current_col_min) ) { // case 2
        max_distance_exceeded = true;
        if constexpr (use_hook) {
          if(!search_add_all(ctx, path(ch.second), child_row_min, hook)) return false;
        } else {
          search_add_all(ctx, path(ch.second), child_row_min, hook);
        }
        break;
      }
    }
    if(!max_distance_exceeded) {
      if constexpr (use_hook) {
        if(!single_gap_search_impl<Hook>(ch.second, node_depth + 1, child_char_depth, child_row_min, ctx, workspace, gap_cost, hook, current_col_min)) return false;
      } else {
        single_gap_search_impl<Hook>(ch.second, node_depth + 1, child_char_depth, child_row_min, ctx, workspace, gap_cost, hook, current_col_min);
      }
    }
  }
  if constexpr (use_hook) return true;
}

inline void RadixMap::mark_removed_subtree(RadixMap::weak_pointer_type node) {
  if(node == nullptr) return;
  node->parent_node = nullptr;
  for(auto ch : node->child_nodes) {
    mark_removed_subtree(ch.second);
  }
}

inline RadixMap::erase_action RadixMap::erase_impl(RadixMap::weak_pointer_type node, const span_type sequence, RadixMap::index_type & result) {
  if(sequence.size() == 0) {
    std::swap(result, node->terminal_idx); // if sequence doesn't exist, terminal_idx should be nullidx which is fine since result is initialized as nullidx
    size_t nc = node->child_nodes.size();
    if(nc == 0) return erase_action::erase;
    if(nc == 1) return erase_action::merge;
    return erase_action::keep;
  }

  // check that sequence actually exists in tree - we shouldn't assume it does
  atomic_type s = sequence[0];
  auto it = node->child_nodes.find(s);
  if(it == node->child_nodes.end()) return erase_action::keep;

  // travelling down tree
  pointer_type & child_ptr = it->second;
  const branch_type & child_branch = child_ptr->branch;
  size_t i = 0;
  for(; i < child_branch.size(); ++i) {
    if(i == sequence.size()) return erase_action::keep; // branch is longer than sequence, doesn't match
    if(child_branch[i] != sequence[i]) return erase_action::keep; // branch and sequence don't match
  }
  erase_action action = erase_impl(child_ptr, sequence.subspan(i), result); // sequence is longer or same

  // travelling back up
  if(action == erase_action::keep) {
    return erase_action::keep;
  } else if(action == erase_action::merge) {
    pointer_type removed_node = child_ptr;
    auto first_it = removed_node->child_nodes.begin();
    appendspan(removed_node->branch, first_it->second->branch);
    branch_type next_branch = std::move(removed_node->branch);
    child_ptr = first_it->second; // old intermediate node is unlinked; pool reclaims it with the tree
    child_ptr->parent_node = node;
    child_ptr->branch = std::move(next_branch);
    removed_node->parent_node = nullptr;
    return erase_action::keep;
  } else { // erase_action::erase
    mark_removed_subtree(child_ptr);
    node->child_nodes.erase(it);
    size_t nc = node->child_nodes.size();
    if((nc == 0) && (node->terminal_idx == nullidx)) {
      return erase_action::erase;
    } else if((nc == 1) && (node->terminal_idx == nullidx)) {
      return erase_action::merge;
    } else {
      return erase_action::keep;
    }
  }
}
template <typename Hook>
inline RadixMap::hook_ret<Hook> RadixMap::hamming_search_impl(RadixMap::const_weak_pointer_type node, const size_t position, const int distance, RadixMap::search_context & ctx, Hook & hook) {
  constexpr bool use_hook = !std::is_same_v<Hook, NullSearchHook>;
  if(position == ctx.query.size()) {
    if(node->terminal_idx != nullidx) {
      if constexpr (use_hook) {
        if(!search_add(ctx, path(node), distance, hook)) return false;
      } else {
        search_add(ctx, path(node), distance, hook);
      }
    }
    if constexpr (use_hook) return true;
    else return;
  }
  if(position < ctx.query.size()) {
    for(const auto& ch : node->child_nodes) {
      const branch_type & branch = ch.second->branch;
      int new_distance = distance;
      if(position + branch.size() > ctx.query.size()) continue;
      bool max_distance_exceeded = false;
      for(size_t j=0; j<branch.size(); ++j) {
        if(branch[j] != ctx.query[position+j]) new_distance++;
        if(new_distance > ctx.max_distance) {
          max_distance_exceeded = true;
          break;
        }
      }
      if(!max_distance_exceeded) {
        if constexpr (use_hook) {
          if(!hamming_search_impl<Hook>(ch.second, position + branch.size(), new_distance, ctx, hook)) return false;
        } else {
          hamming_search_impl<Hook>(ch.second, position + branch.size(), new_distance, ctx, hook);
        }
      }
    }
  }
  if constexpr (use_hook) return true;
}

inline int RadixMap::update_col_banded(const RadixMap::atomic_type branchval, const RadixMap::span_type query,
                                       std::vector<int> & col, size_t lower, size_t upper) {
  size_t query_len = query.size();
  if(upper > query_len) upper = query_len;
  if(lower > upper) {
    std::fill(col.begin(), col.end(), NO_ALIGN);
    return NO_ALIGN;
  }

  int min_element = NO_ALIGN;
  int previous_diag;
  size_t start = lower;
  if(lower == 0) {
    previous_diag = col[0];
    col[0] = col[0] + 1;
    min_element = col[0];
    start = 1;
  } else {
    previous_diag = col[lower - 1];
  }

  for(size_t i = start; i <= upper; ++i) {
    int original = col[i];
    int left = (i == lower) ? NO_ALIGN : col[i-1];
    int match_cost = previous_diag + (query[i-1] == branchval ? 0 : 1);
    int insert_cost = left + 1;
    int delete_cost = original + 1;
    int new_val = match_cost < insert_cost ? match_cost : insert_cost;
    if(delete_cost < new_val) new_val = delete_cost;
    col[i] = new_val;
    previous_diag = original;
    if(new_val < min_element) min_element = new_val;
  }

  if(lower > 0) {
    col[lower - 1] = NO_ALIGN;
  }
  if(upper + 1 < col.size()) {
    col[upper + 1] = NO_ALIGN;
  }

  return min_element;
}

inline int RadixMap::update_col_banded_ptr(const RadixMap::atomic_type branchval, const RadixMap::span_type query,
                                           int * __restrict__ col, size_t col_len, size_t lower, size_t upper) {
  size_t query_len = query.size();
  if(upper > query_len) upper = query_len;
  if(lower > upper) {
    std::fill(col, col + col_len, NO_ALIGN);
    return NO_ALIGN;
  }

  int min_element = NO_ALIGN;
  int previous_diag;
  size_t start = lower;
  if(lower == 0) {
    previous_diag = col[0];
    col[0] = col[0] + 1;
    min_element = col[0];
    start = 1;
  } else {
    previous_diag = col[lower - 1];
  }

  for(size_t i = start; i <= upper; ++i) {
    int original = col[i];
    int left = (i == lower) ? NO_ALIGN : col[i-1];
    int match_cost = previous_diag + (query[i-1] == branchval ? 0 : 1);
    int insert_cost = left + 1;
    int delete_cost = original + 1;
    int new_val = match_cost < insert_cost ? match_cost : insert_cost;
    if(delete_cost < new_val) new_val = delete_cost;
    col[i] = new_val;
    previous_diag = original;
    if(new_val < min_element) min_element = new_val;
  }

  if(lower > 0) {
    col[lower - 1] = NO_ALIGN;
  }
  if(upper + 1 < col_len) {
    col[upper + 1] = NO_ALIGN;
  }

  return min_element;
}

template <typename Hook>
inline RadixMap::hook_ret<Hook> RadixMap::global_search_impl(RadixMap::const_weak_pointer_type node,
                                                             const size_t node_depth,
                                                             const size_t char_depth,
                                                             RadixMap::search_context & ctx,
                                                             UnitFullSimpleArrayWorkspace & workspace,
                                                             Hook & hook,
                                                             int current_col_min) {
  constexpr bool use_hook = !std::is_same_v<Hook, NullSearchHook>;
  const span_type query = ctx.query;
  const size_t query_len = query.size();
  const BandLimits band = ctx.band(char_depth);
  const size_t col_len = workspace.column_len;

  int * current_col = workspace.at(node_depth);
  if(current_col_min > ctx.max_distance) { if constexpr (use_hook) return true; else return; }

  const int terminal_distance = (band.upper == query_len) ? current_col[query_len] : NO_ALIGN;
  if((node->terminal_idx != nullidx) && (terminal_distance <= ctx.max_distance)) {
    if constexpr (use_hook) {
      if(!search_add(ctx, path(node), terminal_distance, hook)) return false;
    } else {
      search_add(ctx, path(node), terminal_distance, hook);
    }
  }

  for(const auto& ch : node->child_nodes) {
    int * child_col = workspace.clone_from_parent(node_depth);
    bool prune_child = false;
    int child_col_min = current_col_min;
    size_t child_char_depth = char_depth;
    BandLimits child_band = band;
    const branch_type & branch = ch.second->branch;
    for(atomic_type branchval : branch) {
      child_band.increment(child_char_depth, ctx.band_radius, query_len);
      ++child_char_depth;
      child_col_min = update_col_banded_ptr(branchval, query, child_col, col_len, child_band.lower, child_band.upper);
      if(child_col_min > ctx.max_distance) {
        prune_child = true;
        break;
      }
    }
    if(!prune_child) {
      if constexpr (use_hook) {
        if(!global_search_impl<Hook>(ch.second, node_depth + 1, child_char_depth, ctx, workspace, hook, child_col_min)) return false;
      } else {
        global_search_impl<Hook>(ch.second, node_depth + 1, child_char_depth, ctx, workspace, hook, child_col_min);
      }
    }
  }

  if constexpr (use_hook) return true;
}

// enumerating all cases for stop conditions
// max > row > col -- if we are on a terminal, add current path (distance = col), keep going (case 3)
// max > col > row -- add all children (distance = row) and stop (case 2)
// row > max > col -- if we are on a terminal, add current path (distance = col), keep going (case 3)
// col > max > row -- add all children (distance = row) and stop (case 2)
// row > col > max -- stop (case 1)
// col > row > max -- stop (case 1)
template <typename Hook>
inline RadixMap::hook_ret<Hook> RadixMap::anchored_search_impl(RadixMap::const_weak_pointer_type node,
                                                              const size_t node_depth,
                                                              const size_t char_depth,
                                                              const int row_min,
                                                              RadixMap::search_context & ctx,
                                                              UnitFullSimpleArrayWorkspace & workspace,
                                                              Hook & hook,
                                                              int current_col_min) {
  constexpr bool use_hook = !std::is_same_v<Hook, NullSearchHook>;
  const size_t query_len = ctx.query.size();
  const BandLimits band = ctx.band(char_depth);

  int current_row_min = row_min;
  if( (current_col_min > ctx.max_distance) && (current_row_min > ctx.max_distance) ) { // case 1
    if constexpr (use_hook) return true; else return;
  } else if( (current_row_min <= ctx.max_distance) && (current_row_min <= current_col_min) ) { // case 2
    if constexpr (use_hook) {
      if(!search_add_all(ctx, path(node), current_row_min, hook)) return false;
      return true;
    } else {
      search_add_all(ctx, path(node), current_row_min, hook);
      return;
    }
  } else if(node->terminal_idx != nullidx) { // case 3
    if(current_col_min <= ctx.max_distance) {
      if constexpr (use_hook) {
        if(!search_add(ctx, path(node), current_col_min, hook)) return false;
      } else {
        search_add(ctx, path(node), current_col_min, hook);
      }
    }
  }
  const int parent_min = current_col_min;
  for(const auto& ch : node->child_nodes) {
    int * child_col = workspace.clone_from_parent(node_depth);
    current_col_min = parent_min;
    const branch_type & branch = ch.second->branch;
    bool max_distance_exceeded = false;
    int child_row_min = row_min;
    size_t child_char_depth = char_depth;
    BandLimits child_band = band;
    for(atomic_type branchval : branch) {
      child_band.increment(child_char_depth, ctx.band_radius, query_len);
      ++child_char_depth;
      current_col_min = update_col_banded_ptr(branchval, ctx.query, child_col, workspace.column_len, child_band.lower, child_band.upper);
      const int terminal_candidate = (child_band.upper == query_len) ? child_col[query_len] : NO_ALIGN;
      child_row_min = std::min(child_row_min, terminal_candidate);
      if( (current_col_min > ctx.max_distance) && (child_row_min > ctx.max_distance) ) { // case 1
        max_distance_exceeded = true;
        break;
      } else if( (child_row_min <= ctx.max_distance) && (child_row_min <= current_col_min) ) { // case 2
        max_distance_exceeded = true;
        if constexpr (use_hook) {
          if(!search_add_all(ctx, path(ch.second), child_row_min, hook)) return false;
        } else {
          search_add_all(ctx, path(ch.second), child_row_min, hook);
        }
        break;
      }
    }
    if(!max_distance_exceeded) {
      if constexpr (use_hook) {
        if(!anchored_search_impl<Hook>(ch.second, node_depth + 1, child_char_depth, child_row_min, ctx, workspace, hook, current_col_min)) return false;
      } else {
        anchored_search_impl<Hook>(ch.second, node_depth + 1, child_char_depth, child_row_min, ctx, workspace, hook, current_col_min);
      }
    }
  }
  if constexpr (use_hook) return true;
}

inline int RadixMap::update_col_linear_banded(const RadixMap::atomic_type branchval,
                                              std::vector<int> & col, size_t lower, size_t upper,
                                              const QueryCostCache & query_costs) {
  size_t query_len = query_costs.query_len;
  if(upper > query_len) upper = query_len;
  if(lower > upper) {
    std::fill(col.begin(), col.end(), NO_ALIGN);
    return NO_ALIGN;
  }

  int min_element = NO_ALIGN;
  int previous_diag;
  size_t start = lower;
  if(lower == 0) {
    previous_diag = col[0];
    col[0] = col[0] + query_costs.gap_cost;
    min_element = col[0];
    start = 1;
  } else {
    previous_diag = col[lower - 1];
  }

  const int * subst_cost = query_costs.subst_for_target(static_cast<char>(branchval));
  for(size_t i = start; i <= upper; ++i) {
    int original = col[i];
    int left = (i == lower) ? NO_ALIGN : col[i-1];
    int match_cost = previous_diag + subst_cost[i - 1];
    int gap_in_query = left + query_costs.gap_cost;
    int gap_in_target = original + query_costs.gap_cost;
    int new_val = min3(match_cost, gap_in_query, gap_in_target);
    col[i] = new_val;
    previous_diag = original;
    if(new_val < min_element) min_element = new_val;
  }

  if(lower > 0) {
    col[lower - 1] = NO_ALIGN;
  }
  if(upper + 1 < col.size()) {
    col[upper + 1] = NO_ALIGN;
  }

  return min_element;
}

inline int RadixMap::update_col_linear_banded_ptr(const RadixMap::atomic_type branchval,
                                                  int * __restrict__ col,
                                                  const size_t col_len,
                                                  size_t lower,
                                                  size_t upper,
                                                  const QueryCostCache & query_costs) {
  size_t query_len = query_costs.query_len;
  if(upper > query_len) upper = query_len;
  if(lower > upper) {
    std::fill(col, col + col_len, NO_ALIGN);
    return NO_ALIGN;
  }

  int min_element = NO_ALIGN;
  int previous_diag;
  size_t start = lower;
  if(lower == 0) {
    previous_diag = col[0];
    col[0] = col[0] + query_costs.gap_cost;
    min_element = col[0];
    start = 1;
  } else {
    previous_diag = col[lower - 1];
  }

  const int * subst_cost = query_costs.subst_for_target(static_cast<char>(branchval));
  for(size_t i = start; i <= upper; ++i) {
    int original = col[i];
    int left = (i == lower) ? NO_ALIGN : col[i-1];
    int match_cost = previous_diag + subst_cost[i - 1];
    int gap_in_query = left + query_costs.gap_cost;
    int gap_in_target = original + query_costs.gap_cost;
    int new_val = min3(match_cost, gap_in_query, gap_in_target);
    col[i] = new_val;
    previous_diag = original;
    if(new_val < min_element) min_element = new_val;
  }

  if(lower > 0) {
    col[lower - 1] = NO_ALIGN;
  }
  if(upper + 1 < col_len) {
    col[upper + 1] = NO_ALIGN;
  }

  return min_element;
}

template <typename Hook>
inline RadixMap::hook_ret<Hook> RadixMap::global_search_linear_impl(RadixMap::const_weak_pointer_type node,
                                                                   const size_t node_depth,
                                                                   const size_t char_depth,
                                                                   search_context_cached & ctx,
                                                                   UnitFullSimpleArrayWorkspace & workspace,
                                                                   Hook & hook,
                                                                   int current_col_min) {
  constexpr bool use_hook = !std::is_same_v<Hook, NullSearchHook>;
  const QueryCostCache & query_costs = ctx.query_costs;
  const size_t query_len = query_costs.query_len;
  const BandLimits band = ctx.band(char_depth);

  int * current_col = workspace.at(node_depth);
  if(current_col_min > ctx.max_distance) { if constexpr (use_hook) return true; else return; }

  const int terminal_distance = (band.upper == query_len) ? current_col[query_len] : NO_ALIGN;
  if((node->terminal_idx != nullidx) && (terminal_distance <= ctx.max_distance)) {
    if constexpr (use_hook) {
      if(!search_add(ctx, path(node), terminal_distance, hook)) return false;
    } else {
      search_add(ctx, path(node), terminal_distance, hook);
    }
  }

  for(const auto& ch : node->child_nodes) {
    int * child_col = workspace.clone_from_parent(node_depth);
    bool prune_child = false;
    int child_col_min = current_col_min;
    size_t child_char_depth = char_depth;
    BandLimits child_band = band;
    const branch_type & branch = ch.second->branch;
    for(atomic_type branchval : branch) {
      child_band.increment(child_char_depth, ctx.band_radius, query_len);
      ++child_char_depth;
      child_col_min = update_col_linear_banded_ptr(branchval, child_col, workspace.column_len, child_band.lower, child_band.upper, query_costs);
      if(child_col_min > ctx.max_distance) {
        prune_child = true;
        break;
      }
    }
    if(!prune_child) {
      if constexpr (use_hook) {
        if(!global_search_linear_impl<Hook>(ch.second, node_depth + 1, child_char_depth, ctx, workspace, hook, child_col_min)) return false;
      } else {
        global_search_linear_impl<Hook>(ch.second, node_depth + 1, child_char_depth, ctx, workspace, hook, child_col_min);
      }
    }
  }
  if constexpr (use_hook) return true;
}

template <typename Hook>
inline RadixMap::hook_ret<Hook> RadixMap::anchored_search_linear_impl(RadixMap::const_weak_pointer_type node,
                                                                     const size_t node_depth,
                                                                     const size_t char_depth,
                                                                     const int row_min,
                                                                     search_context_cached & ctx,
                                                                     UnitFullSimpleArrayWorkspace & workspace,
                                                                     Hook & hook,
                                                                     int current_col_min) {
  constexpr bool use_hook = !std::is_same_v<Hook, NullSearchHook>;
  const QueryCostCache & query_costs = ctx.query_costs;
  const size_t query_len = query_costs.query_len;
  const BandLimits band = ctx.band(char_depth);

  int current_row_min = row_min;
  if( (current_col_min > ctx.max_distance) && (current_row_min > ctx.max_distance) ) { // case 1
    if constexpr (use_hook) return true; else return;
  } else if( (current_row_min <= ctx.max_distance) && (current_row_min <= current_col_min) ) { // case 2
    if constexpr (use_hook) {
      if(!search_add_all(ctx, path(node), current_row_min, hook)) return false;
      return true;
    } else {
      search_add_all(ctx, path(node), current_row_min, hook);
      return;
    }
  } else if(node->terminal_idx != nullidx) { // case 3
    if(current_col_min <= ctx.max_distance) {
      if constexpr (use_hook) {
        if(!search_add(ctx, path(node), current_col_min, hook)) return false;
      } else {
        search_add(ctx, path(node), current_col_min, hook);
      }
    }
  }

  const int parent_min = current_col_min;
  for(const auto& ch : node->child_nodes) {
    int * child_col = workspace.clone_from_parent(node_depth);
    current_col_min = parent_min;
    bool max_distance_exceeded = false;
    int child_row_min = row_min;
    size_t child_char_depth = char_depth;
    BandLimits child_band = band;
    const branch_type & branch = ch.second->branch;
    for(atomic_type branchval : branch) {
      child_band.increment(child_char_depth, ctx.band_radius, query_len);
      ++child_char_depth;
      current_col_min = update_col_linear_banded_ptr(branchval, child_col, workspace.column_len, child_band.lower, child_band.upper, query_costs);
      const int terminal_candidate = (child_band.upper == query_len) ? child_col[query_len] : NO_ALIGN;
      child_row_min = std::min(child_row_min, terminal_candidate);
      if( (current_col_min > ctx.max_distance) && (child_row_min > ctx.max_distance) ) {
        max_distance_exceeded = true;
        break;
      } else if( (child_row_min <= ctx.max_distance) && (child_row_min <= current_col_min) ) {
        max_distance_exceeded = true;
        if constexpr (use_hook) {
          if(!search_add_all(ctx, path(ch.second), child_row_min, hook)) return false;
        } else {
          search_add_all(ctx, path(ch.second), child_row_min, hook);
        }
        break;
      }
    }
    if(!max_distance_exceeded) {
      if constexpr (use_hook) {
        if(!anchored_search_linear_impl<Hook>(ch.second, node_depth + 1, child_char_depth, child_row_min, ctx, workspace, hook, current_col_min)) return false;
      } else {
        anchored_search_linear_impl<Hook>(ch.second, node_depth + 1, child_char_depth, child_row_min, ctx, workspace, hook, current_col_min);
      }
    }
  }
  if constexpr (use_hook) return true;
}

inline int RadixMap::update_col_affine_banded(const RadixMap::atomic_type branchval,
                                              affine_col_type & col,
                                              size_t lower,
                                              size_t upper,
                                              const QueryCostCache & query_costs) {
  auto & M_col =  std::get<0>(col); // match
  auto & X_col = std::get<1>(col); // gap in query
  auto & Y_col = std::get<2>(col); // gap in target

  const size_t col_size = M_col.size();
  if(col_size == 0) return NO_ALIGN;
  if(upper >= col_size) upper = col_size - 1;

  if(lower > upper) {
    std::fill(M_col.begin(), M_col.end(), NO_ALIGN);
    std::fill(X_col.begin(), X_col.end(), NO_ALIGN);
    std::fill(Y_col.begin(), Y_col.end(), NO_ALIGN);
    return NO_ALIGN;
  }

  int min_element = NO_ALIGN;
  int previous_M_i_minus_1;
  int previous_X_i_minus_1;
  int previous_Y_i_minus_1;
  size_t start = lower;

  if(lower == 0) {
    previous_M_i_minus_1 = M_col[0];
    previous_X_i_minus_1 = X_col[0];
    previous_Y_i_minus_1 = Y_col[0];

    M_col[0] = NO_ALIGN;
    X_col[0] = previous_X_i_minus_1 == NO_ALIGN ? query_costs.gap_open_including_first_extension
                                                : previous_X_i_minus_1 + query_costs.gap_cost;
    Y_col[0] = NO_ALIGN;
    min_element = std::min(min_element, X_col[0]);
    start = 1;
  } else {
    previous_M_i_minus_1 = M_col[lower - 1];
    previous_X_i_minus_1 = X_col[lower - 1];
    previous_Y_i_minus_1 = Y_col[lower - 1];
  }

  const int * subst_cost = query_costs.subst_for_target(static_cast<char>(branchval));
  for(size_t i = start; i <= upper; ++i) {
    int original_M = M_col[i];
    int original_X = X_col[i];
    int original_Y = Y_col[i];

    int match_cost = subst_cost[i - 1];
    int M_col_i = match_cost + min3(previous_M_i_minus_1, previous_X_i_minus_1, previous_Y_i_minus_1);
    int X_col_i = min3(query_costs.gap_open_including_first_extension + original_M,
                       query_costs.gap_cost      + original_X,
                       query_costs.gap_open_including_first_extension + original_Y);
    const bool at_band_top = (lower > 0 && i == lower);
    const int prev_M_for_Y = at_band_top ? NO_ALIGN : M_col[i-1];
    const int prev_X_for_Y = at_band_top ? NO_ALIGN : X_col[i-1];
    const int prev_Y_for_Y = at_band_top ? NO_ALIGN : Y_col[i-1];
    int Y_col_i = min3(query_costs.gap_open_including_first_extension + prev_M_for_Y,
                       query_costs.gap_open_including_first_extension + prev_X_for_Y,
                       query_costs.gap_cost      + prev_Y_for_Y);

    previous_M_i_minus_1 = original_M;
    previous_X_i_minus_1 = original_X;
    previous_Y_i_minus_1 = original_Y;

    M_col[i] = M_col_i;
    X_col[i] = X_col_i;
    Y_col[i] = Y_col_i;
    int current_min = min3(M_col_i, X_col_i, Y_col_i);
    if(current_min < min_element) min_element = current_min;
  }

  if(lower > 0) {
    M_col[lower - 1] = NO_ALIGN;
    X_col[lower - 1] = NO_ALIGN;
    Y_col[lower - 1] = NO_ALIGN;
  }
  if(upper + 1 < col_size) {
    M_col[upper + 1] = NO_ALIGN;
    X_col[upper + 1] = NO_ALIGN;
    Y_col[upper + 1] = NO_ALIGN;
  }

  return min_element;
}

template <typename Hook>
inline RadixMap::hook_ret<Hook> RadixMap::global_search_affine_impl(RadixMap::const_weak_pointer_type node,
                                                                    const size_t node_depth,
                                                                    const size_t char_depth,
                                                                    search_context_cached & ctx,
                                                                    AffineWorkspace & workspace,
                                                                    Hook & hook,
                                                                    int current_col_min) {
  constexpr bool use_hook = !std::is_same_v<Hook, NullSearchHook>;
  const QueryCostCache & query_costs = ctx.query_costs;
  const size_t query_len = query_costs.query_len;
  const BandLimits band = ctx.band(char_depth);

  const affine_col_type & parent_column = workspace.at(node_depth);
  const auto & M_col = std::get<0>(parent_column);
  const auto & X_col = std::get<1>(parent_column);
  const auto & Y_col = std::get<2>(parent_column);

  if(current_col_min > ctx.max_distance) { if constexpr (use_hook) return true; else return; }

  int terminal_distance = NO_ALIGN;
  if(band.upper == query_len) {
    terminal_distance = min3(M_col[query_len], X_col[query_len], Y_col[query_len]);
  }
  if((node->terminal_idx != nullidx) && (terminal_distance <= ctx.max_distance)) {
    if constexpr (use_hook) {
      if(!search_add(ctx, path(node), terminal_distance, hook)) return false;
    } else {
      search_add(ctx, path(node), terminal_distance, hook);
    }
  }

  for(const auto& ch : node->child_nodes) {
    affine_col_type & child_col = workspace.clone_from_parent(node_depth);
    bool max_distance_exceeded = false;
    int child_col_min = current_col_min;
    size_t child_char_depth = char_depth;
    BandLimits child_band = band;
    const branch_type & branch = ch.second->branch;
    for(atomic_type branchval : branch) {
      child_band.increment(child_char_depth, ctx.band_radius, query_len);
      ++child_char_depth;
      child_col_min = update_col_affine_banded(branchval, child_col, child_band.lower, child_band.upper, query_costs);
      if(child_col_min > ctx.max_distance) {
        max_distance_exceeded = true;
        break;
      }
    }
    if(!max_distance_exceeded) {
      if constexpr (use_hook) {
        if(!global_search_affine_impl<Hook>(ch.second, node_depth + 1, child_char_depth, ctx, workspace, hook, child_col_min)) return false;
      } else {
        global_search_affine_impl<Hook>(ch.second, node_depth + 1, child_char_depth, ctx, workspace, hook, child_col_min);
      }
    }
  }
  if constexpr (use_hook) return true;
}

template <typename Hook>
inline RadixMap::hook_ret<Hook> RadixMap::anchored_search_affine_impl(RadixMap::const_weak_pointer_type node,
                                                                      const size_t node_depth,
                                                                      const size_t char_depth,
                                                                      const int row_min,
                                                                      search_context_cached & ctx,
                                                                      AffineWorkspace & workspace,
                                                                      Hook & hook,
                                                                      int current_col_min) {
  constexpr bool use_hook = !std::is_same_v<Hook, NullSearchHook>;
  const QueryCostCache & query_costs = ctx.query_costs;
  const size_t query_len = query_costs.query_len;
  const BandLimits band = ctx.band(char_depth);

  int current_row_min = row_min;
  if( (current_col_min > ctx.max_distance) && (current_row_min > ctx.max_distance) ) { // case 1
    if constexpr (use_hook) return true; else return;
  } else if( (current_row_min <= ctx.max_distance) && (current_row_min <= current_col_min) ) { // case 2
    if constexpr (use_hook) {
      if(!search_add_all(ctx, path(node), current_row_min, hook)) return false;
      return true;
    } else {
      search_add_all(ctx, path(node), current_row_min, hook);
      return;
    }
  } else if(node->terminal_idx != nullidx) { // case 3
    if(current_col_min <= ctx.max_distance) {
      if constexpr (use_hook) {
        if(!search_add(ctx, path(node), current_col_min, hook)) return false;
      } else {
        search_add(ctx, path(node), current_col_min, hook);
      }
    }
  }

  const int parent_min = current_col_min;
  for(const auto& ch : node->child_nodes) {
    affine_col_type & child_col = workspace.clone_from_parent(node_depth);
    current_col_min = parent_min;
    bool max_distance_exceeded = false;
    int child_row_min = row_min;
    size_t child_char_depth = char_depth;
    BandLimits child_band = band;
    const branch_type & branch = ch.second->branch;
    for(atomic_type branchval : branch) {
      child_band.increment(child_char_depth, ctx.band_radius, query_len);
      ++child_char_depth;
      current_col_min = update_col_affine_banded(branchval, child_col, child_band.lower, child_band.upper, query_costs);
      int current_col_back = (child_band.upper == query_len)
                               ? min3(std::get<0>(child_col)[query_len],
                                      std::get<1>(child_col)[query_len],
                                      std::get<2>(child_col)[query_len])
                               : NO_ALIGN;
      child_row_min = std::min(child_row_min, current_col_back);
      if( (current_col_min > ctx.max_distance) && (child_row_min > ctx.max_distance) ) { // case 1
        max_distance_exceeded = true;
        break;
      } else if( (child_row_min <= ctx.max_distance) && (child_row_min <= current_col_min) ) { // case 2
        max_distance_exceeded = true;
        if constexpr (use_hook) {
          if(!search_add_all(ctx, path(ch.second), child_row_min, hook)) return false;
        } else {
          search_add_all(ctx, path(ch.second), child_row_min, hook);
        }
        break;
      }
    }
    if(!max_distance_exceeded) {
      if constexpr (use_hook) {
        if(!anchored_search_affine_impl<Hook>(ch.second, node_depth + 1, child_char_depth, child_row_min, ctx, workspace, hook, current_col_min)) return false;
      } else {
        anchored_search_affine_impl<Hook>(ch.second, node_depth + 1, child_char_depth, child_row_min, ctx, workspace, hook, current_col_min);
      }
    }
  }
  if constexpr (use_hook) return true;
}

} // namespace seqtrie

inline bool operator==(const seqtrie::RadixMap::path& a, const seqtrie::RadixMap::path& b) noexcept {
  return a.m == b.m;
}

namespace std {
  template<>
  struct hash<seqtrie::RadixMap::path> {
    size_t operator()(seqtrie::RadixMap::path const& p) const noexcept {
      return std::hash<const seqtrie::RadixMap*>{}(p.m);
    }
  };
}

#endif // seqtrie_RADIXMAP_H
