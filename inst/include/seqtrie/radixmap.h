#ifndef seqtrie_RADIXMAP_H
#define seqtrie_RADIXMAP_H

#include "seqtrie/utility.h"
#include "ankerl/unordered_dense.h"
#include "simple_array/small_array.h"

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


private:
  std::string print_impl(size_t depth) const;
  enum class erase_action { erase, merge, keep };
  static erase_action erase_impl(weak_pointer_type node, const span_type sequence, index_type & result);

  // Implementation for both insert variants
  template <bool ReturnPathAlways>
  path insert_impl(const span_type sequence, index_type idx);

  static void hamming_search_impl(const_weak_pointer_type node, size_t position, int distance, search_context & ctx);
  static int update_col(atomic_type branchval, const span_type query, std::vector<int> & col);
  static void global_search_impl(const_weak_pointer_type node, const std::vector<int> & previous_col, search_context & ctx);
  static void anchored_search_impl(const_weak_pointer_type node, const std::vector<int> & previous_col, int row_min, search_context & ctx);

  static int update_col_linear(atomic_type branchval, const span_type query, std::vector<int> & col, const CostMap & cost_map);
  static void global_search_linear_impl(const_weak_pointer_type node, const std::vector<int> & previous_col, search_context & ctx, const CostMap & cost_map);
  static void anchored_search_linear_impl(const_weak_pointer_type node, const std::vector<int> & previous_col, int row_min, search_context & ctx, const CostMap & cost_map);

  static int update_col_affine(atomic_type branchval, const span_type query, affine_col_type & col, const CostMap & cost_map);
  static void global_search_affine_impl(const_weak_pointer_type node, const affine_col_type & previous_col, search_context & ctx, const CostMap & cost_map);
  static void anchored_search_affine_impl(const_weak_pointer_type node, const affine_col_type & previous_col, int row_min, search_context & ctx, const CostMap & cost_map);

  
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
  global_search_impl(this, iota_range<std::vector<int>>(0, query.size() + 1), ctx); 
  return ctx;
}

// an "anchored" search can end on the last column or col of the dynamic programming array
// unlike global which must end on the bottom right corner
// we need to keep track of the minimum value in the last row
inline RadixMap::search_context RadixMap::anchored_search(const RadixMap::span_type query, const int max_distance) const {
  search_context ctx(query, max_distance);
  anchored_search_impl(this, iota_range<std::vector<int>>(0, query.size() + 1), query.size(), ctx); 
  return ctx;
}

inline RadixMap::search_context RadixMap::global_search_linear(const RadixMap::span_type query, const int max_distance, const CostMap & cost_map) const {
  search_context ctx(query, max_distance);
  std::vector<int> col(query.size() + 1, 0);
  for(size_t i=1; i<col.size(); ++i) {
    col[i] = col[i-1] + cost_map.gap_cost; // gap in target
  }
  global_search_linear_impl(this, col, ctx, cost_map);
  return ctx;
}

inline RadixMap::search_context RadixMap::anchored_search_linear(const RadixMap::span_type query, const int max_distance, const CostMap & cost_map) const {
  search_context ctx(query, max_distance);
  std::vector<int> col(query.size() + 1, 0);
  for(size_t i=1; i<col.size(); ++i) {
    col[i] = col[i-1] + cost_map.gap_cost; // gap in target
  }
  anchored_search_linear_impl(this, col, col.back(), ctx, cost_map);
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
  global_search_affine_impl(this, col, ctx, cost_map);
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
  anchored_search_affine_impl(this, col, 
                                 std::min({M_col.back(), Y_col.back()}), // Edge case: use M_col if query is empty
                                 ctx, cost_map);
  return ctx;
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

inline int RadixMap::update_col(const RadixMap::atomic_type branchval, const RadixMap::span_type query, std::vector<int> & col) {
  int previous_col_i_minus_1 = col[0];
  col[0] = col[0] + 1;
  int min_element = col[0];
  for(size_t i=1; i<col.size(); ++i) {
    int match_cost  = previous_col_i_minus_1 + (query[i-1] == branchval ? 0 : 1);
    int insert_cost = col[i-1] + 1;
    int delete_cost = col[i] + 1;
    previous_col_i_minus_1 = col[i];
    col[i] = std::min({match_cost, insert_cost, delete_cost});
    if(col[i] < min_element) min_element = col[i];
  }
  return min_element;
}

inline void RadixMap::global_search_impl(RadixMap::const_weak_pointer_type node, const std::vector<int> & previous_col, search_context & ctx) {
  if( *std::min_element(previous_col.begin(), previous_col.end()) > ctx.max_distance ) { return; }
  if((node->terminal_idx != nullidx) && (previous_col.back() <= ctx.max_distance)) {
    ctx.match.push_back(path(node));
    ctx.distance.push_back(previous_col.back());
  }
  for (auto & ch : node->child_nodes) {
    std::vector<int> current_col = previous_col;
    branch_type & branch = ch.second->branch;
    bool max_distance_exceeded = false;
    for(size_t u=0; u<branch.size(); ++u) {
      int current_dist = update_col(branch[u], ctx.query, current_col);
      if(current_dist > ctx.max_distance) {
        max_distance_exceeded = true;
        break;
      }
    }
    if(!max_distance_exceeded) global_search_impl(ch.second.get(), current_col, ctx);
  }
}

// enumerating all cases for stop conditions
// max > row > col -- if we are on a terminal, add current path (distance = col), keep going (case 3)
// max > col > row -- add all children (distance = row) and stop (case 2)
// row > max > col -- if we are on a terminal, add current path (distance = col), keep going (case 3)
// col > max > row -- add all children (distance = row) and stop (case 2)
// row > col > max -- stop (case 1)
// col > row > max -- stop (case 1)
inline void RadixMap::anchored_search_impl(RadixMap::const_weak_pointer_type node, const std::vector<int> & previous_col, const int row_min, search_context & ctx) {
  int current_col_min = *std::min_element(previous_col.begin(), previous_col.end());
  int current_row_min = row_min;
  if( (current_col_min > ctx.max_distance) && (current_row_min > ctx.max_distance) ) { // case 1
    return;
  } else if( (current_row_min <= ctx.max_distance) && (current_row_min <= current_col_min) ) { // case 2
    // if row_min <= col_min then the search ends on the column
    // we can stop, since it's impossible to find a smaller value by continuing the search
    std::vector<path> child_sequences = node->all();
    for(auto & chs : child_sequences) { // also includes current node
      if(chs->terminal_idx != nullidx) {
        ctx.match.push_back(chs);
        ctx.distance.push_back(current_row_min);
      }
    }
    return;
  } else if(node->terminal_idx != nullidx) { // this is a terminal leaf and its best value is on the col // case 3
    if(current_col_min <= ctx.max_distance) {
      ctx.match.push_back(path(node));
      ctx.distance.push_back(current_col_min);
    }
  }
  for (auto & ch : node->child_nodes) {
    std::vector<int> current_col = previous_col;
    branch_type & branch = ch.second->branch;
    bool max_distance_exceeded = false;
    current_row_min = row_min;
    for(size_t u=0; u<branch.size(); ++u) {
      current_col_min = update_col(branch[u], ctx.query, current_col);
      current_row_min = std::min(current_row_min, current_col.back());
      if( (current_col_min > ctx.max_distance) && (current_row_min > ctx.max_distance) ) { // case 1
        max_distance_exceeded = true;
        break;
      } else if( (current_row_min <= ctx.max_distance) && (current_row_min <= current_col_min) ) { // case 2
        // if row_min <= col_min then the search ends on the column
        // we can stop, since it's impossible to find a smaller value by continuing the search
        max_distance_exceeded = true; // "exceeded" in the sense that we include all children below and don't need to keep searching
        std::vector<path> grandchild_sequences = ch.second->all();
        for(auto & gchs : grandchild_sequences) { // also includes ch node
          if(gchs->terminal_idx != nullidx) {
            ctx.match.push_back(gchs);
            ctx.distance.push_back(current_row_min);
          }
        }
        break;
      }
      // case 3 does not need to be considered since we are never on top of a terminal leaf within this loop
    }
    if(!max_distance_exceeded) anchored_search_impl(ch.second.get(), current_col, current_row_min, ctx);
  }
}

inline int RadixMap::update_col_linear(const RadixMap::atomic_type branchval, const RadixMap::span_type query, std::vector<int> & col, const CostMap & cost_map) {
  int previous_col_i_minus_1 = col[0];
  col[0] = col[0] + cost_map.gap_cost; // gap in target
  int min_element = col[0];
  for(size_t i=1; i<col.size(); ++i) {
    int match_cost   = previous_col_i_minus_1 + cost_map.char_cost_map.at(std::make_pair(query[i-1], static_cast<char>(branchval)));
    int gap_in_query = col[i]   + cost_map.gap_cost;
    int gap_in_target= col[i-1] + cost_map.gap_cost;
    previous_col_i_minus_1 = col[i];
    col[i] = std::min({match_cost, gap_in_query, gap_in_target});
    if(col[i] < min_element) min_element = col[i];
  }
  return min_element;
}

inline void RadixMap::global_search_linear_impl(RadixMap::const_weak_pointer_type node, const std::vector<int> & previous_col, search_context & ctx, const CostMap & cost_map) {
  if( *std::min_element(previous_col.begin(), previous_col.end()) > ctx.max_distance ) { return; }
  if((node->terminal_idx != nullidx) && (previous_col.back() <= ctx.max_distance)) {
    ctx.match.push_back(path(node));
    ctx.distance.push_back(previous_col.back());
  }
  for (auto & ch : node->child_nodes) {
    std::vector<int> current_col = previous_col;
    branch_type & branch = ch.second->branch;
    bool max_distance_exceeded = false;
    for(size_t u=0; u<branch.size(); ++u) {
      int current_dist = update_col_linear(branch[u], ctx.query, current_col, cost_map);
      if(current_dist > ctx.max_distance) {
        max_distance_exceeded = true;
        break;
      }
    }
    if(!max_distance_exceeded) global_search_linear_impl(ch.second.get(), current_col, ctx, cost_map);
  }
}

inline void RadixMap::anchored_search_linear_impl(RadixMap::const_weak_pointer_type node, const std::vector<int> & previous_col, const int row_min, search_context & ctx, const CostMap & cost_map) {
  int current_col_min = *std::min_element(previous_col.begin(), previous_col.end());
  int current_row_min = row_min;
  if( (current_col_min > ctx.max_distance) && (current_row_min > ctx.max_distance) ) { // case 1
    return;
  } else if( (current_row_min <= ctx.max_distance) && (current_row_min <= current_col_min) ) { // case 2
    // if row_min <= col_min then the search ends on the column
    // we can stop, since it's impossible to find a smaller value by continuing the search
    std::vector<path> child_sequences = node->all();
    for(auto & chs : child_sequences) { // also includes current node
      if(chs->terminal_idx != nullidx) {
        ctx.match.push_back(chs);
        ctx.distance.push_back(current_row_min);
      }
    }
    return;
  } else if(node->terminal_idx != nullidx) { // this is a terminal leaf and its best value is on the col // case 3
    if(current_col_min <= ctx.max_distance) {
      ctx.match.push_back(path(node));
      ctx.distance.push_back(current_col_min);
    }
  }
  for (auto & ch : node->child_nodes) {
    std::vector<int> current_col = previous_col;
    branch_type & branch = ch.second->branch;
    bool max_distance_exceeded = false;
    current_row_min = row_min;
    for(size_t u=0; u<branch.size(); ++u) {
      current_col_min = update_col_linear(branch[u], ctx.query, current_col, cost_map);
      current_row_min = std::min(current_row_min, current_col.back());
      if( (current_col_min > ctx.max_distance) && (current_row_min > ctx.max_distance) ) { // case 1
        max_distance_exceeded = true;
        break;
      } else if( (current_row_min <= ctx.max_distance) && (current_row_min <= current_col_min) ) { // case 2
        // if row_min <= col_min then the search ends on the column
        // we can stop, since it's impossible to find a smaller value by continuing the search
        max_distance_exceeded = true; // "exceeded" in the sense that we include all children below and don't need to keep searching
        std::vector<path> grandchild_sequences = ch.second->all();
        for(auto & gchs : grandchild_sequences) { // also includes ch node
          if(gchs->terminal_idx != nullidx) {
            ctx.match.push_back(gchs);
            ctx.distance.push_back(current_row_min);
          }
        }
        break;
      }
      // case 3 does not need to be considered since we are never on top of a terminal leaf within this loop
    }
    if(!max_distance_exceeded) anchored_search_linear_impl(ch.second.get(), current_col, current_row_min, ctx, cost_map);
  }
}

inline int RadixMap::update_col_affine(const RadixMap::atomic_type branchval, const RadixMap::span_type query, affine_col_type & col, const CostMap & cost_map) {
  auto & M_col =  std::get<0>(col); // match
  auto & X_col = std::get<1>(col); // gap in query
  auto & Y_col = std::get<2>(col); // gap in target

  int previous_M_i_minus_1 = M_col[0];
  int previous_X_i_minus_1 = X_col[0];
  int previous_Y_i_minus_1 = Y_col[0];

  M_col[0] = NO_ALIGN;
  X_col[0] = previous_X_i_minus_1 == NO_ALIGN ? cost_map.gap_open_cost
                                              : previous_X_i_minus_1 + cost_map.gap_cost;
  Y_col[0] = NO_ALIGN;
  int min_element = X_col[0];
  for(size_t i=1; i<M_col.size(); ++i) {
    // Update col[i] after updating {previous_i_minus_1 <- col[i]}
    int M_col_i = cost_map.char_cost_map.at(std::make_pair(query[i-1], static_cast<char>(branchval))) + std::min({
      previous_M_i_minus_1, 
      previous_X_i_minus_1, 
      previous_Y_i_minus_1});
    int X_col_i = std::min({
      cost_map.gap_open_cost + M_col[i],
      cost_map.gap_cost      + X_col[i],
      cost_map.gap_open_cost + Y_col[i]}); 
    int Y_col_i = std::min({
      cost_map.gap_open_cost + M_col[i-1],
      cost_map.gap_open_cost + X_col[i-1],
      cost_map.gap_cost      + Y_col[i-1]});
    previous_M_i_minus_1 = M_col[i];
    previous_X_i_minus_1 = X_col[i];
    previous_Y_i_minus_1 = Y_col[i];
    M_col[i] = M_col_i;
    X_col[i] = X_col_i;
    Y_col[i] = Y_col_i;
    int current_min = std::min({M_col[i], X_col[i], Y_col[i]});
    if(current_min < min_element) min_element = current_min;
  }
  // print_vec(M_col);
  // print_vec(X_col);
  // print_vec(Y_col);
  // std::cout << std::endl;
  return min_element;
}

inline void RadixMap::global_search_affine_impl(RadixMap::const_weak_pointer_type node, const affine_col_type & previous_col, search_context & ctx, const CostMap & cost_map) {
  if(
    (*std::min_element(std::get<0>(previous_col).begin(), std::get<0>(previous_col).end()) > ctx.max_distance) &&
    (*std::min_element(std::get<1>(previous_col).begin(), std::get<1>(previous_col).end()) > ctx.max_distance) &&
    (*std::min_element(std::get<2>(previous_col).begin(), std::get<2>(previous_col).end()) > ctx.max_distance)      
  ) { return; }
  int previous_col_back = std::min({
    std::get<0>(previous_col).back(),
    std::get<1>(previous_col).back(),
    std::get<2>(previous_col).back()});
  if((node->terminal_idx != nullidx) && (previous_col_back <= ctx.max_distance)) {
    ctx.match.push_back(path(node));
    ctx.distance.push_back(previous_col_back);
  }
  for (auto & ch : node->child_nodes) {
    affine_col_type current_col = previous_col;
    branch_type & branch = ch.second->branch;
    bool max_distance_exceeded = false;
    for(size_t u=0; u<branch.size(); ++u) {
      int current_dist = update_col_affine(branch[u], ctx.query, current_col, cost_map);
      if(current_dist > ctx.max_distance) {
        max_distance_exceeded = true;
        break;
      }
    }
    if(!max_distance_exceeded) global_search_affine_impl(ch.second.get(), current_col, ctx, cost_map);
  }
}

inline void RadixMap::anchored_search_affine_impl(RadixMap::const_weak_pointer_type node, const affine_col_type & previous_col, const int row_min, search_context & ctx, const CostMap & cost_map) {
  int current_col_min = std::min({
    *std::min_element(std::get<0>(previous_col).begin(), std::get<0>(previous_col).end()),
    *std::min_element(std::get<1>(previous_col).begin(), std::get<1>(previous_col).end()),
    *std::min_element(std::get<2>(previous_col).begin(), std::get<2>(previous_col).end())});
  int current_row_min = row_min;
  if( (current_col_min > ctx.max_distance) && (current_row_min > ctx.max_distance) ) { // case 1
    return;
  } else if( (current_row_min <= ctx.max_distance) && (current_row_min <= current_col_min) ) { // case 2
    // if row_min <= col_min then the search ends on the column
    // we can stop, since it's impossible to find a smaller value by continuing the search
    std::vector<path> child_sequences = node->all();
    for(auto & chs : child_sequences) { // also includes current node
      if(chs->terminal_idx != nullidx) {
        ctx.match.push_back(chs);
        ctx.distance.push_back(current_row_min);
      }
    }
    return;
  } else if(node->terminal_idx != nullidx) { // this is a terminal leaf and its best value is on the col // case 3
    if(current_col_min <= ctx.max_distance) {
      ctx.match.push_back(path(node));
      ctx.distance.push_back(current_col_min);
    }
  }
  for (auto & ch : node->child_nodes) {
    affine_col_type current_col = previous_col;
    branch_type & branch = ch.second->branch;
    bool max_distance_exceeded = false;
    current_row_min = row_min;
    for(size_t u=0; u<branch.size(); ++u) {
      current_col_min = update_col_affine(branch[u], ctx.query, current_col, cost_map);
      int current_col_back = std::min({
        std::get<0>(current_col).back(),
        std::get<1>(current_col).back(),
        std::get<2>(current_col).back()});
      current_row_min = std::min(current_row_min, current_col_back);
      if( (current_col_min > ctx.max_distance) && (current_row_min > ctx.max_distance) ) { // case 1
        max_distance_exceeded = true;
        break;
      } else if( (current_row_min <= ctx.max_distance) && (current_row_min <= current_col_min) ) { // case 2
        // if row_min <= col_min then the search ends on the column
        // we can stop, since it's impossible to find a smaller value by continuing the search
        max_distance_exceeded = true; // "exceeded" in the sense that we include all children below and don't need to keep searching
        std::vector<path> grandchild_sequences = ch.second->all();
        for(auto & gchs : grandchild_sequences) { // also includes ch node
          if(gchs->terminal_idx != nullidx) {
            ctx.match.push_back(gchs);
            ctx.distance.push_back(current_row_min);
          }
        }
        break;
      }
      // case 3 does not need to be considered since we are never on top of a terminal leaf within this loop
    }
    if(!max_distance_exceeded) anchored_search_affine_impl(ch.second.get(), current_col, current_row_min, ctx, cost_map);
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
