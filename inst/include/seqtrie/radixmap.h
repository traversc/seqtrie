#ifndef seqtrie_RADIXMAP_H
#define seqtrie_RADIXMAP_H

#include "seqtrie/utility.h"
#include "ankerl/unordered_dense.h"
#include "simple_array/small_array.h"
#include <algorithm>
#include <limits>
#include <utility>

#ifndef SEQTRIE_SMALL_ARRAY_SIZE
#error "SEQTRIE_SMALL_ARRAY_SIZE must be defined"
#endif

// template parameters moved to macros
#define CHAR_T char
#define MAP_T ankerl::unordered_dense::map
#define BRANCH_T trqwe::small_array<CHAR_T, std::allocator<CHAR_T>, size_t, std::integral_constant<size_t, SEQTRIE_SMALL_ARRAY_SIZE>>
#define INDEX_T size_t

namespace seqtrie {

class RadixMap {
public:
  // fundamental types
  typedef CHAR_T                                   atomic_type;
  typedef BRANCH_T                                 branch_type;
  typedef INDEX_T                                  index_type;
  typedef size_t                                  size_type;
  typedef RadixMap                                self_type;
  typedef std::unique_ptr<self_type>              pointer_type;
  typedef const self_type*                        const_weak_pointer_type;
  typedef self_type*                              weak_pointer_type;
  typedef MAP_T<CHAR_T,pointer_type>              map_type;
  typedef std::pair<CHAR_T,CHAR_T>                pairchar_type;
  typedef nonstd::span<const CHAR_T>              span_type;
  typedef std::tuple<std::vector<int>,std::vector<int>,std::vector<int>> affine_col_type;

  // constants
  static constexpr INDEX_T nullidx     = std::numeric_limits<INDEX_T>::max();
  static constexpr atomic_type GAP_CHAR      = CHAR_T(0);                                      // '\0' gap for non-affine
  static constexpr atomic_type GAP_OPEN_CHAR = std::numeric_limits<atomic_type>::min();        // '\255' gap open for affine
  static constexpr atomic_type GAP_EXTN_CHAR = CHAR_T(0);                                      // '\0' gap extension for affine
  static constexpr int         NO_ALIGN     = std::numeric_limits<int>::max() / 2;             // impossible positions

private:
  map_type                child_nodes;
  branch_type             branch;
  const_weak_pointer_type parent_node;
  index_type              terminal_idx;

public:
  RadixMap() : parent_node(nullptr), terminal_idx(nullidx) {}

  struct path {
    const_weak_pointer_type m;
    path() : m(nullptr) {}
    path(const_weak_pointer_type x) : m(x) {}
    const_weak_pointer_type operator->() const { return m; }
  };

  struct search_context {
    std::vector<path> match;
    std::vector<int> distance;
    span_type query;
    int max_distance;
    search_context() {}
    search_context(span_type query, int max_distance) : query(query), max_distance(max_distance) {}
    void append(const search_context & other) {
      match.insert(match.end(), other.match.begin(), other.match.end());
      distance.insert(distance.end(), other.distance.begin(), other.distance.end());
      query = other.query;
      max_distance = other.max_distance;
    }
  };

  struct UnitWorkspace {
    std::vector<std::vector<int>> columns;

    void initialize(std::vector<int> initial, size_t depth_hint = 0) {
      columns.clear();
      if(depth_hint > 0) {
        if(columns.capacity() < depth_hint) columns.reserve(depth_hint);
      }
      columns.emplace_back(std::move(initial));
    }

    void ensure_child_slot(size_t depth) {
      size_t required_size = depth + 2;
      if(columns.size() < required_size) {
        columns.resize(required_size);
      }
    }

    std::vector<int> & at(const size_t depth) { return columns[depth]; }
    const std::vector<int> & at(const size_t depth) const { return columns[depth]; }

    std::vector<int> & clone_from_parent(size_t depth) {
      ensure_child_slot(depth);
      const auto & parent_col = columns[depth];
      std::vector<int> & child_col = columns[depth + 1];
      child_col.assign(parent_col.begin(), parent_col.end());
      return child_col;
    }
  };

  struct AffineWorkspace {
    std::vector<affine_col_type> columns;

    void initialize(affine_col_type initial, size_t depth_hint = 0) {
      columns.clear();
      if(depth_hint > 0) {
        if(columns.capacity() < depth_hint) columns.reserve(depth_hint);
      }
      columns.emplace_back(std::move(initial));
    }

    void ensure_child_slot(size_t depth) {
      size_t required_size = depth + 2;
      if(columns.size() < required_size) {
        columns.resize(required_size);
      }
    }

    affine_col_type & at(const size_t depth) { return columns[depth]; }
    const affine_col_type & at(const size_t depth) const { return columns[depth]; }

    affine_col_type & clone_from_parent(size_t depth) {
      ensure_child_slot(depth);
      affine_col_type & child_col = columns[depth + 1];
      child_col = columns[depth];
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

  struct BandLimits {
    size_t lower;
    size_t upper;

    static inline size_t linear_radius(int max_distance, int gap_cost, size_t query_len) {
      if(max_distance <= 0) return 0;
      if(gap_cost <= 0) return query_len;
      return static_cast<size_t>(max_distance / gap_cost);
    }

    static inline size_t affine_radius(int max_distance, int gap_cost, int gap_open_cost, size_t query_len) {
      if(max_distance < gap_open_cost) return 0;
      const int remaining = max_distance - gap_open_cost;
      return 1 + static_cast<size_t>(remaining / gap_cost);
    }

    static inline BandLimits around(size_t char_depth, size_t band_radius, size_t query_len) {
      size_t lower = (char_depth > band_radius) ? char_depth - band_radius : 0;
      if(lower > query_len) lower = query_len;
      const size_t upper = std::min(query_len, char_depth + band_radius);
      return {lower, upper};
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
  std::vector<path> all(size_t max_depth = -1) const;
  std::string print() const { return print_impl(0); }
  std::pair<std::vector<path>, std::vector<path>> graph(size_t max_depth = -1) const;
  path find(const span_type query) const;
  path insert(const span_type sequence, index_type idx);
  // Ensure sequence exists and return its terminal path
  path insert_get_path(const span_type sequence, index_type idx);
  index_type erase(const span_type sequence);
  std::vector<path> prefix_search(const span_type query) const;

  search_context hamming_search(const span_type query, int max_distance) const;
  search_context global_search(const span_type query, int max_distance) const;
  search_context anchored_search(const span_type query, int max_distance) const;

  search_context global_search_linear(const span_type query, int max_distance, const CostMap & cost_map) const;
  search_context anchored_search_linear(const span_type query, int max_distance, const CostMap & cost_map) const;

  search_context global_search_affine(const span_type query, int max_distance, const CostMap & cost_map) const;
  search_context anchored_search_affine(const span_type query, int max_distance, const CostMap & cost_map) const;
  search_context single_gap_search(const span_type query, int max_distance, const int gap_cost) const;

private:
  std::string print_impl(size_t depth) const;
  enum class erase_action { erase, merge, keep };
  static erase_action erase_impl(weak_pointer_type node, const span_type sequence, index_type & result);

  // Implementation functions
  template <bool ReturnPathAlways>
  path insert_impl(const span_type sequence, index_type idx);

  static void hamming_search_impl(const_weak_pointer_type node, size_t position, int distance, search_context & ctx);
  static int update_col_banded(atomic_type branchval, const span_type query, std::vector<int> & col, size_t lower, size_t upper);
  static void global_search_impl(const_weak_pointer_type node, size_t node_depth, size_t char_depth, search_context & ctx, UnitWorkspace & workspace); // unitary cost
  static void anchored_search_impl(const_weak_pointer_type node, size_t node_depth, size_t char_depth, int row_min, search_context & ctx, UnitWorkspace & workspace); // unitary cost

  static int update_col_linear_banded(atomic_type branchval, const span_type query, std::vector<int> & col, size_t lower, size_t upper, const CostMap & cost_map);
  static void global_search_linear_impl(const_weak_pointer_type node, size_t node_depth, size_t char_depth, search_context & ctx, UnitWorkspace & workspace, const CostMap & cost_map);
  static void anchored_search_linear_impl(const_weak_pointer_type node, size_t node_depth, size_t char_depth, int row_min, search_context & ctx, UnitWorkspace & workspace, const CostMap & cost_map);

  static int update_col_affine_banded(atomic_type branchval, const span_type query, affine_col_type & col, size_t lower, size_t upper, const CostMap & cost_map);
  static void global_search_affine_impl(const_weak_pointer_type node, size_t node_depth, size_t char_depth, search_context & ctx, AffineWorkspace & workspace, const CostMap & cost_map);
  static void anchored_search_affine_impl(const_weak_pointer_type node, size_t node_depth, size_t char_depth, int row_min, search_context & ctx, AffineWorkspace & workspace, const CostMap & cost_map);

  static int update_col_single_gap(atomic_type branchval, const span_type query, size_t ref_len, SingleGapCol & col, const int gap_cost);
  static void single_gap_search_impl(const_weak_pointer_type node, size_t node_depth, size_t char_depth, int row_min, search_context & ctx, SingleGapWorkspace & workspace, const int gap_cost);
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
  for(auto & ch : child_nodes) {
    size_t next_depth = (max_depth == std::numeric_limits<size_t>::max()) ? max_depth : (max_depth - 1);
    std::vector<path> x = ch.second->all(next_depth);
    appendspan(result, x);
  }
  return result;
}

inline RadixMap::size_type RadixMap::size() const {
  size_type result = terminal_idx == nullidx ? 0 : 1;
  for(auto & ch : child_nodes) { result += ch.second->size(); }
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
  for(auto & ch : child_nodes) {
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
  for(auto & ch : child_nodes) {
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
  for(auto & ch : child_nodes) {
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
    if(node->child_nodes.find(query[position]) != node->child_nodes.end()) {
      node = node->child_nodes.at(query[position]).get();
      if(position + node->branch.size() > query.size()) return path();
      for(size_t j=0; j<node->branch.size(); ++j) {
        if(node->branch[j] != query[position+j]) return path();
      }
      position += node->branch.size();
    } else {
      return path();
    }
  }
  return path(node);
}

inline RadixMap::path RadixMap::insert(const RadixMap::span_type sequence, index_type idx) {
  return insert_impl<false>(sequence, idx);
}

inline RadixMap::path RadixMap::insert_get_path(const RadixMap::span_type sequence, index_type idx) {
  return insert_impl<true>(sequence, idx);
}

template <bool ReturnPathAlways>
inline RadixMap::path RadixMap::insert_impl(const RadixMap::span_type sequence, index_type idx) {
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
  if(child_nodes.find(s) == child_nodes.end()) {
    child_nodes.emplace(s, pointer_type(new self_type));
    child_nodes[s]->parent_node = this;
    child_nodes[s]->branch = subvector<branch_type>(sequence, 0);
    child_nodes[s]->terminal_idx = idx;
    if constexpr (ReturnPathAlways) {
      return path(child_nodes[s].get());
    } else {
      return path();
    }
  }
  size_t i = 0;
  while(i < child_nodes[s]->branch.size() && i < sequence.size() && sequence[i] == child_nodes[s]->branch[i]) { ++i; }
  
  if(i == sequence.size() && i == child_nodes[s]->branch.size()) {
    if(child_nodes[s]->terminal_idx == nullidx) {
      child_nodes[s]->terminal_idx = idx;
      if constexpr (ReturnPathAlways) {
        return path(child_nodes[s].get());
      } else {
        return path();
      }
    } else {
      return path(child_nodes[s].get());
    }
  } else if(i == sequence.size()) {
    branch_type branch_prefix = subvector<branch_type>(child_nodes[s]->branch,0,i);
    branch_type branch_suffix = subvector<branch_type>(child_nodes[s]->branch,i);
    atomic_type s_insert = branch_suffix[0];
    pointer_type inserted_node(new self_type);
    inserted_node->parent_node = this;
    inserted_node->child_nodes[s_insert] = std::move(child_nodes[s]);
    inserted_node->child_nodes[s_insert]->parent_node = inserted_node.get();
    inserted_node->branch = std::move(branch_prefix);
    inserted_node->terminal_idx = idx;
    inserted_node->child_nodes[s_insert]->branch = std::move(branch_suffix);
    child_nodes[s] = std::move(inserted_node);
    if constexpr (ReturnPathAlways) {
      return path(child_nodes[s].get());
    } else {
      return path();
    }
  } else if(i == child_nodes[s]->branch.size()) {
    span_type seq_suffix = sequence.subspan(i);
    return child_nodes[s]->insert_impl<ReturnPathAlways>(seq_suffix, idx);
  } else {
    branch_type branch_prefix = subvector<branch_type>(child_nodes[s]->branch,0,i);
    branch_type branch_suffix = subvector<branch_type>(child_nodes[s]->branch,i);
    branch_type seq_suffix = subvector<branch_type>(sequence,i);
    atomic_type s_insert_branch = branch_suffix[0];
    atomic_type s_insert_seq = seq_suffix[0];
    
    pointer_type inserted_node(new self_type);
    inserted_node->parent_node = this;
    inserted_node->child_nodes[s_insert_branch] = std::move(child_nodes[s]);
    inserted_node->child_nodes[s_insert_branch]->parent_node = inserted_node.get();
    inserted_node->child_nodes[s_insert_branch]->branch = std::move(branch_suffix);
    inserted_node->child_nodes[s_insert_seq] = pointer_type(new self_type);
    inserted_node->child_nodes[s_insert_seq]->parent_node = inserted_node.get();
    inserted_node->child_nodes[s_insert_seq]->branch = std::move(seq_suffix);
    inserted_node->child_nodes[s_insert_seq]->terminal_idx = idx;
    child_nodes[s] = std::move(inserted_node);
    child_nodes[s]->branch = std::move(branch_prefix);
    if constexpr (ReturnPathAlways) {
      return path(child_nodes[s]->child_nodes[s_insert_seq].get());
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
      if(node->child_nodes.find(query[query_position]) != node->child_nodes.end()) {
        node = node->child_nodes.at(query[query_position]).get();
        branch_position = 0;
      } else {
        return std::vector<path>(0);
      }
    }
    if(node->branch[branch_position] == query[query_position]) {
      branch_position++;
      query_position++;
    } else {
      return std::vector<path>(0);
    }
  }
  return node->all();
}

inline RadixMap::search_context RadixMap::hamming_search(const span_type query, const int max_distance) const {
  search_context ctx(query, max_distance);
  hamming_search_impl(this, 0, 0, ctx);
  return ctx;
}

inline RadixMap::search_context RadixMap::global_search(const RadixMap::span_type query, const int max_distance) const {
  search_context ctx(query, max_distance);
  UnitWorkspace workspace;
  workspace.initialize(iota_range<std::vector<int>>(0, query.size() + 1), query.size() + 1);
  global_search_impl(this, 0, 0, ctx, workspace); 
  return ctx;
}

// an "anchored" search can end on the last column or col of the dynamic programming array
// unlike global which must end on the bottom right corner
// we need to keep track of the minimum value in the last row
inline RadixMap::search_context RadixMap::anchored_search(const RadixMap::span_type query, const int max_distance) const {
  search_context ctx(query, max_distance);
  UnitWorkspace workspace;
  workspace.initialize(iota_range<std::vector<int>>(0, query.size() + 1), query.size() + 1);
  anchored_search_impl(this, 0, 0, query.size(), ctx, workspace); 
  return ctx;
}

inline RadixMap::search_context RadixMap::global_search_linear(const RadixMap::span_type query, const int max_distance, const CostMap & cost_map) const {
  search_context ctx(query, max_distance);
  std::vector<int> col(query.size() + 1, 0);
  for(size_t i=1; i<col.size(); ++i) {
    col[i] = col[i-1] + cost_map.gap_cost; // gap in target
  }
  UnitWorkspace workspace;
  workspace.initialize(std::move(col), query.size() + 1);
  global_search_linear_impl(this, 0, 0, ctx, workspace, cost_map);
  return ctx;
}

inline RadixMap::search_context RadixMap::anchored_search_linear(const RadixMap::span_type query, const int max_distance, const CostMap & cost_map) const {
  search_context ctx(query, max_distance);
  std::vector<int> col(query.size() + 1, 0);
  for(size_t i=1; i<col.size(); ++i) {
    col[i] = col[i-1] + cost_map.gap_cost; // gap in target
  }
  int row_min = col.back();
  UnitWorkspace workspace;
  workspace.initialize(std::move(col), query.size() + 1);
  anchored_search_linear_impl(this, 0, 0, row_min, ctx, workspace, cost_map);
  return ctx;
}

inline RadixMap::search_context RadixMap::global_search_affine(const RadixMap::span_type query, const int max_distance, const CostMap & cost_map) const {
  search_context ctx(query, max_distance);
  size_t col_size = query.size() + 1;
  affine_col_type col = std::make_tuple(std::vector<int>(col_size, 0), std::vector<int>(col_size, 0), std::vector<int>(col_size, 0));
  auto & M_col =  std::get<0>(col);
  auto & X_col = std::get<1>(col);
  auto & Y_col = std::get<2>(col);
  M_col[0] = 0;
  X_col[0] = NO_ALIGN;
  Y_col[0] = NO_ALIGN;
  for(size_t i=1; i<col_size; ++i) {
    M_col[i] = NO_ALIGN;
    X_col[i] = NO_ALIGN;
    if(i == 1) Y_col[i] = cost_map.gap_open_cost;
    else       Y_col[i] = Y_col[i-1] + cost_map.gap_cost;
  }
  // print_vec(M_col);
  // print_vec(X_col);
  // print_vec(Y_col);
  // std::cout << std::endl;
  AffineWorkspace workspace;
  workspace.initialize(std::move(col), query.size() + 1);
  global_search_affine_impl(this, 0, 0, ctx, workspace, cost_map);
  return ctx;
}

 

inline RadixMap::search_context RadixMap::anchored_search_affine(const RadixMap::span_type query, const int max_distance, const CostMap & cost_map) const {
  search_context ctx(query, max_distance);
  size_t col_size = query.size() + 1;
  affine_col_type col = std::make_tuple(std::vector<int>(col_size, 0), std::vector<int>(col_size, 0), std::vector<int>(col_size, 0));
  auto & M_col =  std::get<0>(col);
  auto & X_col = std::get<1>(col);
  auto & Y_col = std::get<2>(col);
  M_col[0] = 0;
  X_col[0] = NO_ALIGN;
  Y_col[0] = NO_ALIGN;
  for(size_t i=1; i<col_size; ++i) {
    M_col[i] = NO_ALIGN;
    X_col[i] = NO_ALIGN;
    if(i == 1) Y_col[i] = cost_map.gap_open_cost;
    else       Y_col[i] = Y_col[i-1] + cost_map.gap_cost;
  }
  int row_min = std::min({M_col.back(), Y_col.back()});
  AffineWorkspace workspace;
  workspace.initialize(std::move(col), query.size() + 1);
  anchored_search_affine_impl(this, 0, 0, row_min, ctx, workspace, cost_map);
  return ctx;
}

inline RadixMap::search_context RadixMap::single_gap_search(const RadixMap::span_type query,
                                                            const int max_distance,
                                                            const int gap_cost) const {
  search_context ctx(query, max_distance);
  SingleGapWorkspace workspace;
  SingleGapCol initial{NO_ALIGN, 0, (query.size() > 0 ? gap_cost : 0)}; // if query is empty, lower stores the row min which is 0
  workspace.initialize(initial, query.size() + 1);
  int row_min = query.size() == 0 ? 0 :
                query.size() == 1 ? gap_cost :
                NO_ALIGN;
  single_gap_search_impl(this, 0, 0, row_min, ctx, workspace, gap_cost);
  return ctx;
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
  return std::min({col.diag, col.upper, col.lower});
}

inline void RadixMap::single_gap_search_impl(RadixMap::const_weak_pointer_type node,
                                             const size_t node_depth,
                                             const size_t char_depth,
                                             int row_min,
                                             RadixMap::search_context & ctx,
                                             SingleGapWorkspace & workspace,
                                             const int gap_cost) {
  workspace.ensure_child_slot(node_depth);
  const size_t query_len = ctx.query.size();

  SingleGapCol & current_col = workspace.at(node_depth);
  int current_col_min = std::min({current_col.diag, current_col.upper, current_col.lower});
  auto get_row_min = [query_len](const SingleGapCol & col, size_t char_depth, int row_min) -> int {
    if(char_depth == query_len - 1) {
      return col.upper;
    } else if(char_depth == query_len) {
      return std::min({row_min, col.diag});
    } else if(char_depth == query_len + 1) {
      return std::min({row_min, col.lower});
    } else {
      return NO_ALIGN;
    }
  };
  int current_row_min = get_row_min(current_col, char_depth, row_min);
  if( (current_col_min > ctx.max_distance) && (current_row_min > ctx.max_distance) ) { // case 1
    return;
  } else if( (current_row_min <= ctx.max_distance) && (current_row_min <= current_col_min) ) { // case 2
    std::vector<path> child_sequences = node->all();
    for(auto & chs : child_sequences) {
      if(chs->terminal_idx != nullidx) {
        ctx.match.push_back(chs);
        ctx.distance.push_back(current_row_min);
      }
    }
    return;
  } else if(node->terminal_idx != nullidx) { // case 3
    if(current_col_min <= ctx.max_distance) {
      ctx.match.push_back(path(node));
      ctx.distance.push_back(current_col_min);
    }
  }
  const int parent_min = current_col_min;
  for (auto & ch : node->child_nodes) {
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
        std::vector<path> grandchild_sequences = ch.second->all();
        for(auto & gchs : grandchild_sequences) {
          if(gchs->terminal_idx != nullidx) {
            ctx.match.push_back(gchs);
            ctx.distance.push_back(child_row_min);
          }
        }
        break;
      }
    }
    if(!max_distance_exceeded) {
      single_gap_search_impl(ch.second.get(), node_depth + 1, child_char_depth, child_row_min, ctx, workspace, gap_cost);
    }
  }
}

inline RadixMap::erase_action RadixMap::erase_impl(RadixMap::weak_pointer_type node, const span_type sequence, RadixMap::index_type & result) {
  if(sequence.size() == 0) {
    std::swap(result, node->terminal_idx); // if sequence doesn't exist, terminal_idx should be nullidx which is fine since result is initialized as nullidx
    size_t nc = node->child_nodes.size();
    if(nc == 0) { // no children
      return erase_action::erase;
    } else if(nc == 1) { // one child
      return erase_action::merge;
    } else { // two or more
      return erase_action::keep;
    }
  }
  
  // check that sequence actually exists in tree - we shouldn't assume it does
  atomic_type s = sequence[0];
  if(node->child_nodes.find(s) == node->child_nodes.end()) {
    return erase_action::keep;
  }
  
  // travelling down tree
  size_t i = 0;
  for(; i<node->child_nodes[s]->branch.size(); ++i) {
    if(i == sequence.size()) return erase_action::keep; // branch is longer than sequence, doesn't match
    if(node->child_nodes[s]->branch[i] != sequence[i]) return erase_action::keep; // branch and sequence don't match
  }
  erase_action action = erase_impl(node->child_nodes[s].get(), sequence.subspan(i), result); // sequence is longer or same
  
  // travelling back up
  if(action == erase_action::keep) {
    return erase_action::keep;
  } else if(action == erase_action::merge) {
    size_t next_s = 0;
    for(auto & x : node->child_nodes[s]->child_nodes) {
      next_s = x.first;
      break;
    }
    appendspan(node->child_nodes[s]->branch, node->child_nodes[s]->child_nodes[next_s]->branch);
    branch_type next_branch = std::move(node->child_nodes[s]->branch);
    node->child_nodes[s] = std::move(node->child_nodes[s]->child_nodes[next_s]);
    node->child_nodes[s]->parent_node = node;
    node->child_nodes[s]->branch = std::move(next_branch);
    return erase_action::keep;
  } else { // if(action == erase_action::erase) {
    node->child_nodes.erase(s);
    size_t nc = node->child_nodes.size();
    if((nc == 0) && (node->terminal_idx == nullidx)) { // no children and not a sequence
      return erase_action::erase;
    } else if((nc == 1) && (node->terminal_idx == nullidx)) { // one child and not a sequence
      return erase_action::merge;
    } else { // 2+ children or node is a sequence
      return erase_action::keep;
    }
  }
}

inline void RadixMap::hamming_search_impl(RadixMap::const_weak_pointer_type node, const size_t position, const int distance, RadixMap::search_context & ctx) {
  if(position == ctx.query.size()) {
    if(node->terminal_idx != nullidx) {
      ctx.match.push_back(path(node));
      ctx.distance.push_back(distance);
    }
    return;
  }
  if(position < ctx.query.size()) {
    for (auto & ch : node->child_nodes) {
      branch_type & branch = ch.second->branch;
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
        hamming_search_impl(ch.second.get(), position + branch.size(), new_distance, ctx);
      }
    }
  }
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
    int new_val = std::min({match_cost, insert_cost, delete_cost});
    col[i] = new_val;
    previous_diag = original;
    if(new_val < min_element) min_element = new_val;
  }

  if(lower > 0) {
    col[lower - 1] = NO_ALIGN;
  }
  for(size_t i = upper + 1; i < col.size(); ++i) col[i] = NO_ALIGN;

  return min_element;
}

inline void RadixMap::global_search_impl(RadixMap::const_weak_pointer_type node,
                                         const size_t node_depth,
                                         const size_t char_depth,
                                         RadixMap::search_context & ctx,
                                         UnitWorkspace & workspace) {
  workspace.ensure_child_slot(node_depth);
  const size_t query_len = ctx.query.size();
  const size_t band_radius = static_cast<size_t>(ctx.max_distance);
  const BandLimits band = BandLimits::around(char_depth, band_radius, query_len);

  std::vector<int> & current_col = workspace.at(node_depth);
  int current_col_min = NO_ALIGN;
  if(band.lower <= band.upper) {
    for(size_t i = band.lower; i <= band.upper; ++i) {
      current_col_min = std::min(current_col_min, current_col[i]);
    }
  }
  if(current_col_min > ctx.max_distance) { return; }

  const int terminal_distance = (band.upper == query_len) ? current_col.back() : NO_ALIGN;
  if((node->terminal_idx != nullidx) && (terminal_distance <= ctx.max_distance)) {
    ctx.match.push_back(path(node));
    ctx.distance.push_back(terminal_distance);
  }

  for (auto & ch : node->child_nodes) {
    std::vector<int> & child_col = workspace.clone_from_parent(node_depth);
    bool prune_child = false;
    size_t child_char_depth = char_depth;
    const branch_type & branch = ch.second->branch;
    for(atomic_type branchval : branch) {
      ++child_char_depth;
      const BandLimits child_band = BandLimits::around(child_char_depth, band_radius, query_len);
      const int current_dist = update_col_banded(branchval, ctx.query, child_col, child_band.lower, child_band.upper);
      if(current_dist > ctx.max_distance) {
        prune_child = true;
        break;
      }
    }
    if(!prune_child) {
      global_search_impl(ch.second.get(), node_depth + 1, child_char_depth, ctx, workspace);
    }
  }
}

// enumerating all cases for stop conditions
// max > row > col -- if we are on a terminal, add current path (distance = col), keep going (case 3)
// max > col > row -- add all children (distance = row) and stop (case 2)
// row > max > col -- if we are on a terminal, add current path (distance = col), keep going (case 3)
// col > max > row -- add all children (distance = row) and stop (case 2)
// row > col > max -- stop (case 1)
// col > row > max -- stop (case 1)
inline void RadixMap::anchored_search_impl(RadixMap::const_weak_pointer_type node,
                                           const size_t node_depth,
                                           const size_t char_depth,
                                           const int row_min,
                                           RadixMap::search_context & ctx,
                                           UnitWorkspace & workspace) {
  workspace.ensure_child_slot(node_depth);
  const size_t query_len = ctx.query.size();
  const size_t band_radius = static_cast<size_t>(ctx.max_distance);
  const BandLimits band = BandLimits::around(char_depth, band_radius, query_len);

  std::vector<int> & current_col = workspace.at(node_depth);
  int current_col_min = NO_ALIGN;
  if(band.lower <= band.upper) {
    for(size_t i = band.lower; i <= band.upper; ++i) {
      current_col_min = std::min(current_col_min, current_col[i]);
    }
  }

  int current_row_min = row_min;
  if( (current_col_min > ctx.max_distance) && (current_row_min > ctx.max_distance) ) { // case 1
    return;
  } else if( (current_row_min <= ctx.max_distance) && (current_row_min <= current_col_min) ) { // case 2
    std::vector<path> child_sequences = node->all();
    for(auto & chs : child_sequences) {
      if(chs->terminal_idx != nullidx) {
        ctx.match.push_back(chs);
        ctx.distance.push_back(current_row_min);
      }
    }
    return;
  } else if(node->terminal_idx != nullidx) { // case 3
    if(current_col_min <= ctx.max_distance) {
      ctx.match.push_back(path(node));
      ctx.distance.push_back(current_col_min);
    }
  }
  const int parent_min = current_col_min;
  for (auto & ch : node->child_nodes) {
    std::vector<int> & child_col = workspace.clone_from_parent(node_depth);
    current_col_min = parent_min;
    const branch_type & branch = ch.second->branch;
    bool max_distance_exceeded = false;
    int child_row_min = row_min;
    size_t child_char_depth = char_depth;
    for(atomic_type branchval : branch) {
      ++child_char_depth;
      const BandLimits child_band = BandLimits::around(child_char_depth, band_radius, query_len);
      current_col_min = update_col_banded(branchval, ctx.query, child_col, child_band.lower, child_band.upper);
      const int terminal_candidate = (child_band.upper == query_len) ? child_col.back() : NO_ALIGN;
      child_row_min = std::min(child_row_min, terminal_candidate);
      if( (current_col_min > ctx.max_distance) && (child_row_min > ctx.max_distance) ) { // case 1
        max_distance_exceeded = true;
        break;
      } else if( (child_row_min <= ctx.max_distance) && (child_row_min <= current_col_min) ) { // case 2
        max_distance_exceeded = true;
        std::vector<path> grandchild_sequences = ch.second->all();
        for(auto & gchs : grandchild_sequences) {
          if(gchs->terminal_idx != nullidx) {
            ctx.match.push_back(gchs);
            ctx.distance.push_back(child_row_min);
          }
        }
        break;
      }
    }
    if(!max_distance_exceeded) anchored_search_impl(ch.second.get(), node_depth + 1, child_char_depth, child_row_min, ctx, workspace);
  }
}

inline int RadixMap::update_col_linear_banded(const RadixMap::atomic_type branchval, const RadixMap::span_type query,
                                              std::vector<int> & col, size_t lower, size_t upper,
                                              const CostMap & cost_map) {
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
    col[0] = col[0] + cost_map.gap_cost;
    min_element = col[0];
    start = 1;
  } else {
    previous_diag = col[lower - 1];
  }

  for(size_t i = start; i <= upper; ++i) {
    int original = col[i];
    int left = (i == lower) ? NO_ALIGN : col[i-1];
    int match_cost = previous_diag + cost_map.char_cost_map.at(std::make_pair(query[i-1], static_cast<char>(branchval)));
    int gap_in_query = left + cost_map.gap_cost;
    int gap_in_target = original + cost_map.gap_cost;
    int new_val = std::min({match_cost, gap_in_query, gap_in_target});
    col[i] = new_val;
    previous_diag = original;
    if(new_val < min_element) min_element = new_val;
  }

  if(lower > 0) {
    col[lower - 1] = NO_ALIGN;
  }
  for(size_t i = upper + 1; i < col.size(); ++i) col[i] = NO_ALIGN;

  return min_element;
}

inline void RadixMap::global_search_linear_impl(RadixMap::const_weak_pointer_type node,
                                                const size_t node_depth,
                                                const size_t char_depth,
                                                search_context & ctx,
                                                UnitWorkspace & workspace,
                                                const CostMap & cost_map) {
  workspace.ensure_child_slot(node_depth);
  const size_t query_len = ctx.query.size();
  const size_t band_radius = BandLimits::linear_radius(ctx.max_distance, cost_map.gap_cost, query_len);
  const BandLimits band = BandLimits::around(char_depth, band_radius, query_len);

  std::vector<int> & current_col = workspace.at(node_depth);
  int current_col_min = NO_ALIGN;
  if(band.lower <= band.upper) {
    for(size_t i = band.lower; i <= band.upper; ++i) {
      current_col_min = std::min(current_col_min, current_col[i]);
    }
  }
  if(current_col_min > ctx.max_distance) { return; }

  const int terminal_distance = (band.upper == query_len) ? current_col.back() : NO_ALIGN;
  if((node->terminal_idx != nullidx) && (terminal_distance <= ctx.max_distance)) {
    ctx.match.push_back(path(node));
    ctx.distance.push_back(terminal_distance);
  }

  for (auto & ch : node->child_nodes) {
    std::vector<int> & child_col = workspace.clone_from_parent(node_depth);
    bool prune_child = false;
    size_t child_char_depth = char_depth;
    const branch_type & branch = ch.second->branch;
    for(atomic_type branchval : branch) {
      ++child_char_depth;
      const BandLimits child_band = BandLimits::around(child_char_depth, band_radius, query_len);
      const int current_dist = update_col_linear_banded(branchval, ctx.query, child_col, child_band.lower, child_band.upper, cost_map);
      if(current_dist > ctx.max_distance) {
        prune_child = true;
        break;
      }
    }
    if(!prune_child) {
      global_search_linear_impl(ch.second.get(), node_depth + 1, child_char_depth, ctx, workspace, cost_map);
    }
  }
}

inline void RadixMap::anchored_search_linear_impl(RadixMap::const_weak_pointer_type node,
                                                  const size_t node_depth,
                                                  const size_t char_depth,
                                                  const int row_min,
                                                  search_context & ctx,
                                                  UnitWorkspace & workspace,
                                                  const CostMap & cost_map) {
  workspace.ensure_child_slot(node_depth);
  const size_t query_len = ctx.query.size();
  const size_t band_radius = BandLimits::linear_radius(ctx.max_distance, cost_map.gap_cost, query_len);
  const BandLimits band = BandLimits::around(char_depth, band_radius, query_len);

  std::vector<int> & current_col = workspace.at(node_depth);
  int current_col_min = NO_ALIGN;
  if(band.lower <= band.upper) {
    for(size_t i = band.lower; i <= band.upper; ++i) {
      current_col_min = std::min(current_col_min, current_col[i]);
    }
  }

  int current_row_min = row_min;
  if( (current_col_min > ctx.max_distance) && (current_row_min > ctx.max_distance) ) { // case 1
    return;
  } else if( (current_row_min <= ctx.max_distance) && (current_row_min <= current_col_min) ) { // case 2
    std::vector<path> child_sequences = node->all();
    for(auto & chs : child_sequences) {
      if(chs->terminal_idx != nullidx) {
        ctx.match.push_back(chs);
        ctx.distance.push_back(current_row_min);
      }
    }
    return;
  } else if(node->terminal_idx != nullidx) { // case 3
    if(current_col_min <= ctx.max_distance) {
      ctx.match.push_back(path(node));
      ctx.distance.push_back(current_col_min);
    }
  }

  const int parent_min = current_col_min;
  for (auto & ch : node->child_nodes) {
    std::vector<int> & child_col = workspace.clone_from_parent(node_depth);
    current_col_min = parent_min;
    bool max_distance_exceeded = false;
    int child_row_min = row_min;
    size_t child_char_depth = char_depth;
    const branch_type & branch = ch.second->branch;
    for(atomic_type branchval : branch) {
      ++child_char_depth;
      const BandLimits child_band = BandLimits::around(child_char_depth, band_radius, query_len);
      current_col_min = update_col_linear_banded(branchval, ctx.query, child_col, child_band.lower, child_band.upper, cost_map);
      const int terminal_candidate = (child_band.upper == query_len) ? child_col[query_len] : NO_ALIGN;
      child_row_min = std::min(child_row_min, terminal_candidate);
      if( (current_col_min > ctx.max_distance) && (child_row_min > ctx.max_distance) ) {
        max_distance_exceeded = true;
        break;
      } else if( (child_row_min <= ctx.max_distance) && (child_row_min <= current_col_min) ) {
        max_distance_exceeded = true;
        std::vector<path> grandchild_sequences = ch.second->all();
        for(auto & gchs : grandchild_sequences) {
          if(gchs->terminal_idx != nullidx) {
            ctx.match.push_back(gchs);
            ctx.distance.push_back(child_row_min);
          }
        }
        break;
      }
    }
    if(!max_distance_exceeded) {
      anchored_search_linear_impl(ch.second.get(), node_depth + 1, child_char_depth, child_row_min, ctx, workspace, cost_map);
    }
  }
}

inline int RadixMap::update_col_affine_banded(const RadixMap::atomic_type branchval,
                                              const RadixMap::span_type query,
                                              affine_col_type & col,
                                              size_t lower,
                                              size_t upper,
                                              const CostMap & cost_map) {
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
    X_col[0] = previous_X_i_minus_1 == NO_ALIGN ? cost_map.gap_open_cost
                                                : previous_X_i_minus_1 + cost_map.gap_cost;
    Y_col[0] = NO_ALIGN;
    min_element = std::min(min_element, X_col[0]);
    start = 1;
  } else {
    previous_M_i_minus_1 = M_col[lower - 1];
    previous_X_i_minus_1 = X_col[lower - 1];
    previous_Y_i_minus_1 = Y_col[lower - 1];
    for(size_t i = 0; i < lower; ++i) {
      M_col[i] = NO_ALIGN;
      X_col[i] = NO_ALIGN;
      Y_col[i] = NO_ALIGN;
    }
  }

  for(size_t i = start; i <= upper; ++i) {
    int original_M = M_col[i];
    int original_X = X_col[i];
    int original_Y = Y_col[i];

    int match_cost = cost_map.char_cost_map.at(std::make_pair(query[i-1], static_cast<char>(branchval)));
    int M_col_i = match_cost + std::min({previous_M_i_minus_1, previous_X_i_minus_1, previous_Y_i_minus_1});
    int X_col_i = std::min({
      cost_map.gap_open_cost + original_M,
      cost_map.gap_cost      + original_X,
      cost_map.gap_open_cost + original_Y});
    const bool at_band_top = (lower > 0 && i == lower);
    const int prev_M_for_Y = at_band_top ? NO_ALIGN : M_col[i-1];
    const int prev_X_for_Y = at_band_top ? NO_ALIGN : X_col[i-1];
    const int prev_Y_for_Y = at_band_top ? NO_ALIGN : Y_col[i-1];
    int Y_col_i = std::min({
      cost_map.gap_open_cost + prev_M_for_Y,
      cost_map.gap_open_cost + prev_X_for_Y,
      cost_map.gap_cost      + prev_Y_for_Y});

    previous_M_i_minus_1 = original_M;
    previous_X_i_minus_1 = original_X;
    previous_Y_i_minus_1 = original_Y;

    M_col[i] = M_col_i;
    X_col[i] = X_col_i;
    Y_col[i] = Y_col_i;
    int current_min = std::min({M_col_i, X_col_i, Y_col_i});
    if(current_min < min_element) min_element = current_min;
  }

  for(size_t i = upper + 1; i < col_size; ++i) {
    M_col[i] = NO_ALIGN;
    X_col[i] = NO_ALIGN;
    Y_col[i] = NO_ALIGN;
  }

  return min_element;
}

inline void RadixMap::global_search_affine_impl(RadixMap::const_weak_pointer_type node,
                                                const size_t node_depth,
                                                const size_t char_depth,
                                                search_context & ctx,
                                                AffineWorkspace & workspace,
                                                const CostMap & cost_map) {
  workspace.ensure_child_slot(node_depth);
  const size_t query_len = ctx.query.size();
  const size_t band_radius = BandLimits::affine_radius(ctx.max_distance, cost_map.gap_cost, cost_map.gap_open_cost, query_len);
  const BandLimits band = BandLimits::around(char_depth, band_radius, query_len);

  const affine_col_type & parent_column = workspace.at(node_depth);
  const auto & M_col = std::get<0>(parent_column);
  const auto & X_col = std::get<1>(parent_column);
  const auto & Y_col = std::get<2>(parent_column);

  int current_col_min = NO_ALIGN;
  if(band.lower <= band.upper) {
    for(size_t i = band.lower; i <= band.upper; ++i) {
      current_col_min = std::min(current_col_min, std::min({M_col[i], X_col[i], Y_col[i]}));
    }
  }
  if(current_col_min > ctx.max_distance) { return; }

  int terminal_distance = NO_ALIGN;
  if(band.upper == query_len) {
    terminal_distance = std::min({M_col[query_len], X_col[query_len], Y_col[query_len]});
  }
  if((node->terminal_idx != nullidx) && (terminal_distance <= ctx.max_distance)) {
    ctx.match.push_back(path(node));
    ctx.distance.push_back(terminal_distance);
  }

  for (auto & ch : node->child_nodes) {
    affine_col_type & child_col = workspace.clone_from_parent(node_depth);
    bool max_distance_exceeded = false;
    size_t child_char_depth = char_depth;
    const branch_type & branch = ch.second->branch;
    for(atomic_type branchval : branch) {
      ++child_char_depth;
      const BandLimits child_band = BandLimits::around(child_char_depth, band_radius, query_len);
      int current_dist = update_col_affine_banded(branchval, ctx.query, child_col, child_band.lower, child_band.upper, cost_map);
      if(current_dist > ctx.max_distance) {
        max_distance_exceeded = true;
        break;
      }
    }
    if(!max_distance_exceeded) {
      global_search_affine_impl(ch.second.get(), node_depth + 1, child_char_depth, ctx, workspace, cost_map);
    }
  }
}

inline void RadixMap::anchored_search_affine_impl(RadixMap::const_weak_pointer_type node,
                                                  const size_t node_depth,
                                                  const size_t char_depth,
                                                  const int row_min,
                                                  search_context & ctx,
                                                  AffineWorkspace & workspace,
                                                  const CostMap & cost_map) {
  workspace.ensure_child_slot(node_depth);
  const size_t query_len = ctx.query.size();
  const size_t band_radius = BandLimits::affine_radius(ctx.max_distance, cost_map.gap_cost, cost_map.gap_open_cost, query_len);
  const BandLimits band = BandLimits::around(char_depth, band_radius, query_len);

  const affine_col_type & parent_column = workspace.at(node_depth);
  const auto & M_col = std::get<0>(parent_column);
  const auto & X_col = std::get<1>(parent_column);
  const auto & Y_col = std::get<2>(parent_column);

  int current_col_min = NO_ALIGN;
  if(band.lower <= band.upper) {
    for(size_t i = band.lower; i <= band.upper; ++i) {
      current_col_min = std::min(current_col_min, std::min({M_col[i], X_col[i], Y_col[i]}));
    }
  }

  int current_row_min = row_min;
  if( (current_col_min > ctx.max_distance) && (current_row_min > ctx.max_distance) ) { // case 1
    return;
  } else if( (current_row_min <= ctx.max_distance) && (current_row_min <= current_col_min) ) { // case 2
    std::vector<path> child_sequences = node->all();
    for(auto & chs : child_sequences) {
      if(chs->terminal_idx != nullidx) {
        ctx.match.push_back(chs);
        ctx.distance.push_back(current_row_min);
      }
    }
    return;
  } else if(node->terminal_idx != nullidx) { // case 3
    if(current_col_min <= ctx.max_distance) {
      ctx.match.push_back(path(node));
      ctx.distance.push_back(current_col_min);
    }
  }

  const int parent_min = current_col_min;
  for (auto & ch : node->child_nodes) {
    affine_col_type & child_col = workspace.clone_from_parent(node_depth);
    current_col_min = parent_min;
    bool max_distance_exceeded = false;
    int child_row_min = row_min;
    size_t child_char_depth = char_depth;
    const branch_type & branch = ch.second->branch;
    for(atomic_type branchval : branch) {
      ++child_char_depth;
      const BandLimits child_band = BandLimits::around(child_char_depth, band_radius, query_len);
      current_col_min = update_col_affine_banded(branchval, ctx.query, child_col, child_band.lower, child_band.upper, cost_map);
      int current_col_back = (child_band.upper == query_len)
                               ? std::min({
                                   std::get<0>(child_col)[query_len],
                                   std::get<1>(child_col)[query_len],
                                   std::get<2>(child_col)[query_len]})
                               : NO_ALIGN;
      child_row_min = std::min(child_row_min, current_col_back);
      if( (current_col_min > ctx.max_distance) && (child_row_min > ctx.max_distance) ) { // case 1
        max_distance_exceeded = true;
        break;
      } else if( (child_row_min <= ctx.max_distance) && (child_row_min <= current_col_min) ) { // case 2
        max_distance_exceeded = true;
        std::vector<path> grandchild_sequences = ch.second->all();
        for(auto & gchs : grandchild_sequences) {
          if(gchs->terminal_idx != nullidx) {
            ctx.match.push_back(gchs);
            ctx.distance.push_back(child_row_min);
          }
        }
        break;
      }
    }
    if(!max_distance_exceeded) {
      anchored_search_affine_impl(ch.second.get(), node_depth + 1, child_char_depth, child_row_min, ctx, workspace, cost_map);
    }
  }
}

#undef CHAR_T
#undef MAP_T
#undef BRANCH_T
#undef INDEX_T
#undef GAP_CHAR
#undef GAP_OPEN_CHAR
#undef GAP_EXTN_CHAR
#undef NO_ALIGN

} // namespace seqtrie

#endif // seqtrie_RADIXMAP_H
