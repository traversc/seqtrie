#ifndef seqtrie_RADIXMAP_H
#define seqtrie_RADIXMAP_H

#include "seqtrie/utility.h"

namespace seqtrie {

inline void print_vec(const std::vector<int> & v) {
  for(auto & x : v) {
    if(x > 100000) {
      std::cout << "Z ";
    } else {
      std::cout << x << " ";
    }
  }
  std::cout << "\t";
}

#define TEMPLATE_LIST template <class A, template<typename...> class M, template<typename...> class B, class I>
#define RADIXMAP_T RadixMap<A,M,B,I>
#define GAP_CHAR 0                                                        // '\0' any gap cost for non-affine
#define GAP_OPEN_CHAR std::numeric_limits<A>::min()                       // '\255' gap open cost for affine
#define GAP_EXTN_CHAR 0                                                   // '\0' extension for affine
#define NO_ALIGN std::numeric_limits<typename MT::mapped_type>::max() / 2 // represents impossible positions in affine alignment, MT::mapped_type should be int

template <class A = char, template<typename...> class M = std::map, template<typename...> class B = std::vector, class I=size_t> class RadixMap;
template <class A = char, template<typename...> class M = std::map, template<typename...> class B = std::vector, class I=size_t> using RadixMapUPtr = std::unique_ptr<RadixMap<A,M,B,I>>;
template <class A, template<typename...> class M, template<typename...> class B, class I>
class RadixMap {
public:
  static constexpr I nullidx = std::numeric_limits<I>::max();
  typedef A atomic_type; // '\0' and '\255' are reserved character for gaps
  typedef B<A> branch_type;
  typedef I index_type;
  typedef size_t size_type;
  typedef RadixMap<A,M,B,I> self_type;
  typedef RadixMapUPtr<A,M,B,I> pointer_type;
  typedef RadixMap<A,M,B,I> const * const_weak_pointer_type;
  typedef RadixMap<A,M,B,I> * weak_pointer_type;
  typedef M<A,RadixMapUPtr<A,M,B,I>> map_type;
  typedef std::pair<A,A> pairchar_type;
  typedef nonstd::span<const A> span_type;
  typedef std::tuple<std::vector<int>, std::vector<int>, std::vector<int>> affine_col_type;
private:
  map_type child_nodes;        // 48 bytes for std::map
  branch_type branch;         // 24 bytes for std::vector
  const_weak_pointer_type parent_node; // 8 bytes
  index_type terminal_idx;    // 8 bytes
public:
  RadixMap() : parent_node(nullptr), terminal_idx(nullidx) {}
  struct path { // like an iterator, but without iteration. Later we could turn this into a real iterator. 
    const_weak_pointer_type m;
    path() : m(nullptr) {}
    path(const_weak_pointer_type x) : m(x) {}
    const_weak_pointer_type operator->() const { return m; }
  };
  struct search_context { // result struct for global and hamming search
    std::vector<path> match;
    std::vector<int> distance;
    span_type query;
    int max_distance;
    search_context() {}
    search_context(span_type query, int max_distance) : query(query), max_distance(max_distance) {}
    void append(const search_context & other) { // convienence method for processing similar results
      this->match.insert(this->match.end(), other.match.begin(), other.match.end());
      this->distance.insert(this->distance.end(), other.distance.begin(), other.distance.end());
      this->query = other.query;
      this->max_distance = other.max_distance;
    }
  };
  const map_type & get_child_nodes() const { return child_nodes; }
  const branch_type & get_branch() const { return branch; }
  const_weak_pointer_type get_parent_node() const { return parent_node; }
  index_type get_terminal_idx() const { return terminal_idx; }

  // output the sequence represented by this node, by traversing up the tree
  // ST is an array-like container type, e.g. std::vector<char>, std::string
  template <typename ST> ST sequence() const;

  size_type size() const;
  bool validate(const bool is_root = true) const; // check if tree makes sense
  std::vector<path> all(size_t max_depth = -1) const;
  std::string print() const { return print_impl(0); } // this function really only makes sense if atomic_type is char...
  std::pair<std::vector<path>, std::vector<path>> graph(size_t max_depth = -1) const;
  index_type find(const span_type query) const; // returns terminal_idx if found, nullidx if not found
  index_type insert(const span_type sequence, index_type idx); // returns terminal_idx if already exists, otherwise nullidx
  index_type erase(const span_type sequence); // returns terminal_idx if sequence existed, nullidx if it did not exist
  std::vector<path> prefix_search(const span_type query) const;

  search_context hamming_search(const span_type query, const int max_distance) const;
  search_context global_search(const span_type query, const int max_distance) const;
  search_context anchored_search(const span_type query, const int max_distance) const;

  // search using a custom edit distance cost matrix
  // MT is a map type, e.g. std::map<pairchar_type, int>
  // (pairchar_type is defined as std::pair<char, char> in utility.h)
  // the map key is the pair of characters (query, target) to compare
  // All possibile pairs of characters must be included in the map, or else it return an error
  // The map value is the cost of the edit operation and must be non-negative
  // The map must also include the special character \0, which represents a gap
  // Example:
  // map<pairchar_type, int> m = {{pairchar_type('A','A'), 0}, {pairchar_type('A','B'), 1}, {pairchar_type('A','\0'), 1},
//                           {pairchar_type('B','A'), 1}, {pairchar_type('B','B'), 0}, {pairchar_type('B','\0'), 1} };
  template <typename MT> search_context global_search_linear(const span_type query, const int max_distance, const MT & cost_map) const;
  template <typename MT> search_context anchored_search_linear(const span_type query, const int max_distance, const MT & cost_map) const;
  
  template <typename MT> search_context global_search_affine(const span_type query, const int max_distance, const MT & cost_map) const;
  template <typename MT> search_context anchored_search_affine(const span_type query, const int max_distance, const MT & cost_map) const;
  


private:
  // implementation helpers, subject to change
  std::string print_impl(size_t depth) const;
  enum class erase_action { erase, merge, keep };
  static erase_action erase_impl(weak_pointer_type node, const span_type sequence, index_type & result);

  static void hamming_search_impl(const_weak_pointer_type node, const size_t position, const int distance, search_context & ctx);

  static int update_col(const atomic_type branchval, const span_type query, std::vector<int> & col); // helper for global and anchored search
  static void global_search_impl(const_weak_pointer_type node, const std::vector<int> & previous_col, search_context & ctx);
  static void anchored_search_impl(const_weak_pointer_type node, const std::vector<int> & previous_col, const int row_min, search_context & ctx);

  template <typename MT> static int update_col_linear(const atomic_type branchval, const span_type query, std::vector<int> & col, const MT & cost_map);
  template <typename MT> static void global_search_linear_impl(const_weak_pointer_type node, const std::vector<int> & previous_col, search_context & ctx, const MT & cost_map);
  template <typename MT> static void anchored_search_linear_impl(const_weak_pointer_type node, const std::vector<int> & previous_col, const int row_min, search_context & ctx, const MT & cost_map);

  template <typename MT> static int update_col_affine(const atomic_type branchval, const span_type query, affine_col_type & col, const MT & cost_map);
  template <typename MT> static void global_search_affine_impl(const_weak_pointer_type node, const affine_col_type & previous_col, search_context & ctx, const MT & cost_map);
  template <typename MT> static void anchored_search_affine_impl(const_weak_pointer_type node, const affine_col_type & previous_col, const int row_min, search_context & ctx, const MT & cost_map);
};

TEMPLATE_LIST template <typename ST>
ST RADIXMAP_T::sequence() const {
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

TEMPLATE_LIST std::vector<typename RADIXMAP_T::path> RADIXMAP_T::all(size_t max_depth) const {
  std::vector<path> result;
  if(terminal_idx != nullidx) {
    result.push_back(this);
  }
  if(max_depth == 0) return result;
  for(auto & ch : child_nodes) {
    std::vector<path> x = ch.second->all(--max_depth);
    appendspan(result, x);
  }
  return result;
}

TEMPLATE_LIST typename RADIXMAP_T::size_type RADIXMAP_T::size() const {
  size_type result = terminal_idx == nullidx ? 0 : 1;
  for(auto & ch : child_nodes) { result += ch.second->size(); }
  return result;
}

TEMPLATE_LIST bool RADIXMAP_T::validate(const bool is_root) const {
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
    return ch.second->validate(false);
  }
  return true;
}


TEMPLATE_LIST std::string RADIXMAP_T::print_impl(size_t depth) const {
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

TEMPLATE_LIST std::pair<std::vector<typename RADIXMAP_T::path>, std::vector<typename RADIXMAP_T::path>> RADIXMAP_T::graph(size_t max_depth) const {
  std::pair<std::vector<path>, std::vector<path>> result;
  if(parent_node != nullptr) {
    result.first.push_back(path(parent_node));
    result.second.push_back(path(this));
  }
  if(max_depth == 0) return result;
  for(auto & ch : child_nodes) {
    auto x = ch.second->graph(--max_depth);
    appendspan(result.first, x.first);
    appendspan(result.second, x.second);
  }
  return result;
}

TEMPLATE_LIST typename RADIXMAP_T::index_type RADIXMAP_T::find(const typename RADIXMAP_T::span_type query) const {
  const_weak_pointer_type node = this;
  size_t position=0;
  while(position < query.size()) {
    if(node->child_nodes.find(query[position]) != node->child_nodes.end()) {
      node = node->child_nodes.at(query[position]).get();
      if(position + node->branch.size() > query.size()) return nullidx;
      for(size_t j=0; j<node->branch.size(); ++j) {
        if(node->branch[j] != query[position+j]) return nullidx;
      }
      position += node->branch.size();
    } else {
      return nullidx;
    }
  }
  return node->terminal_idx;
}

TEMPLATE_LIST typename RADIXMAP_T::index_type RADIXMAP_T::insert(const typename RADIXMAP_T::span_type sequence, index_type idx) {
  if(sequence.size() == 0) {
    // std::cout << "case -1: end of sequence" << std::endl;
    if(terminal_idx == nullidx) {
      terminal_idx = idx;
      return nullidx;
    } else {
      return terminal_idx;
    }
  }
  atomic_type s = sequence[0];
  if(child_nodes.find(s) == child_nodes.end()) {
    // std::cout << "case 0: no existing branch" << std::endl;
    child_nodes.emplace(s, pointer_type(new self_type));
    child_nodes[s]->parent_node = this;
    child_nodes[s]->branch = subvector<branch_type>(sequence, 0);
    child_nodes[s]->terminal_idx = idx;
    return nullidx;
  }
  size_t i = 0;
  while(i < child_nodes[s]->branch.size() && i < sequence.size() && sequence[i] == child_nodes[s]->branch[i]) { ++i; }
  
  if(i == sequence.size() && i == child_nodes[s]->branch.size()) {
    // std::cout << "case 1: sequence is same as branch" << std::endl;
    if(child_nodes[s]->terminal_idx == nullidx) {
      child_nodes[s]->terminal_idx = idx;
      return nullidx;
    } else {
      return child_nodes[s]->terminal_idx;
    }
  } else if(i == sequence.size()) {
    // std::cout << "case 2: sequence is shorter than branch" << std::endl;
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
    return nullidx;
  } else if(i == child_nodes[s]->branch.size()) {
    // std::cout << "case 3: sequence is longer than branch" << std::endl;
    span_type seq_suffix = sequence.subspan(i); // will remain as span
    return child_nodes[s]->insert(seq_suffix, idx);
  } else {
    // std::cout << "case 4: sequence is different than branch" << std::endl;
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
    return nullidx;
  }
}


TEMPLATE_LIST typename RADIXMAP_T::index_type RADIXMAP_T::erase(const typename RADIXMAP_T::span_type sequence) {
  index_type result = nullidx;
  erase_impl(this, sequence, result);
  return result;
}

TEMPLATE_LIST std::vector<typename RADIXMAP_T::path> RADIXMAP_T::prefix_search(const typename RADIXMAP_T::span_type query) const {
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

TEMPLATE_LIST typename RADIXMAP_T::search_context RADIXMAP_T::hamming_search(const span_type query, const int max_distance) const {
  search_context ctx(query, max_distance);
  hamming_search_impl(this, 0, 0, ctx);
  return ctx;
}

TEMPLATE_LIST typename RADIXMAP_T::search_context RADIXMAP_T::global_search(const typename RADIXMAP_T::span_type query, const int max_distance) const {
  search_context ctx(query, max_distance);
  global_search_impl(this, iota_range<std::vector<int>>(0, query.size() + 1), ctx); 
  return ctx;
}

// an "anchored" search can end on the last column or col of the dynamic programming array
// unlike global which must end on the bottom right corner
// we need to keep track of the minimum value in the last row
TEMPLATE_LIST typename RADIXMAP_T::search_context RADIXMAP_T::anchored_search(const typename RADIXMAP_T::span_type query, const int max_distance) const {
  search_context ctx(query, max_distance);
  anchored_search_impl(this, iota_range<std::vector<int>>(0, query.size() + 1), query.size(), ctx); 
  return ctx;
}

TEMPLATE_LIST template <typename MT>
typename RADIXMAP_T::search_context RADIXMAP_T::global_search_linear(const typename RADIXMAP_T::span_type query, const int max_distance, const MT & cost_map) const {
  search_context ctx(query, max_distance);
  std::vector<int> col(query.size() + 1, 0);
  for(size_t i=1; i<col.size(); ++i) {
    col[i] = col[i-1] + cost_map.at(pairchar_type(query[i-1], GAP_CHAR)); // gap in target
  }
  global_search_linear_impl<MT>(this, col, ctx, cost_map);
  return ctx;
}

TEMPLATE_LIST template <typename MT>
typename RADIXMAP_T::search_context RADIXMAP_T::anchored_search_linear(const typename RADIXMAP_T::span_type query, const int max_distance, const MT & cost_map) const {
  search_context ctx(query, max_distance);
  std::vector<int> col(query.size() + 1, 0);
  for(size_t i=1; i<col.size(); ++i) {
    col[i] = col[i-1] + cost_map.at(pairchar_type(query[i-1], GAP_CHAR)); // gap in target
  }
  anchored_search_linear_impl<MT>(this, col, col.back(), ctx, cost_map);
  return ctx;
}

TEMPLATE_LIST template <typename MT>
typename RADIXMAP_T::search_context RADIXMAP_T::global_search_affine(const typename RADIXMAP_T::span_type query, const int max_distance, const MT & cost_map) const {
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
    if(i == 1) {
      Y_col[i] = cost_map.at(pairchar_type(query[i-1], GAP_OPEN_CHAR));
    } else {
      Y_col[i] = Y_col[i-1] + cost_map.at(pairchar_type(query[i-1], GAP_EXTN_CHAR));
    }
  }
  // print_vec(M_col);
  // print_vec(X_col);
  // print_vec(Y_col);
  // std::cout << std::endl;
  global_search_affine_impl<MT>(this, col, ctx, cost_map);
  return ctx;
}

TEMPLATE_LIST template <typename MT>
typename RADIXMAP_T::search_context RADIXMAP_T::anchored_search_affine(const typename RADIXMAP_T::span_type query, const int max_distance, const MT & cost_map) const {
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
    if(i == 1) {
      Y_col[i] = cost_map.at(pairchar_type(query[i-1], GAP_OPEN_CHAR));
    } else {
      Y_col[i] = Y_col[i-1] + cost_map.at(pairchar_type(query[i-1], GAP_EXTN_CHAR));
    }
  }
  anchored_search_affine_impl<MT>(this, col, 
                                 std::min({M_col.back(), Y_col.back()}), // Edge case: use M_col if query is empty
                                 ctx, cost_map);
  return ctx;
}

TEMPLATE_LIST typename RADIXMAP_T::erase_action RADIXMAP_T::erase_impl(typename RADIXMAP_T::weak_pointer_type node, const span_type sequence, typename RADIXMAP_T::index_type & result) {
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

TEMPLATE_LIST void RADIXMAP_T::hamming_search_impl(typename RADIXMAP_T::const_weak_pointer_type node, const size_t position, const int distance, typename RADIXMAP_T::search_context & ctx) {
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

TEMPLATE_LIST int RADIXMAP_T::update_col(const typename RADIXMAP_T::atomic_type branchval, const typename RADIXMAP_T::span_type query, std::vector<int> & col) {
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

TEMPLATE_LIST void RADIXMAP_T::global_search_impl(typename RADIXMAP_T::const_weak_pointer_type node, const std::vector<int> & previous_col, search_context & ctx) {
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
TEMPLATE_LIST void RADIXMAP_T::anchored_search_impl(typename RADIXMAP_T::const_weak_pointer_type node, const std::vector<int> & previous_col, const int row_min, search_context & ctx) {
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

TEMPLATE_LIST template <typename MT>
int RADIXMAP_T::update_col_linear(const typename RADIXMAP_T::atomic_type branchval, const typename RADIXMAP_T::span_type query, std::vector<int> & col, const MT & cost_map) {
  int previous_col_i_minus_1 = col[0];
  col[0] = col[0] + cost_map.at(pairchar_type(GAP_CHAR, branchval)); // gap in target
  int min_element = col[0];
  for(size_t i=1; i<col.size(); ++i) {
    int match_cost  = previous_col_i_minus_1 + cost_map.at(pairchar_type(query[i-1], branchval));
    int gap_in_query = col[i] + cost_map.at(pairchar_type(GAP_CHAR, branchval));
    int gap_in_target = col[i-1] + cost_map.at(pairchar_type(query[i-1], GAP_CHAR));
    previous_col_i_minus_1 = col[i];
    col[i] = std::min({match_cost, gap_in_query, gap_in_target});
    if(col[i] < min_element) min_element = col[i];
  }
  return min_element;
}

TEMPLATE_LIST template <typename MT>
void RADIXMAP_T::global_search_linear_impl(typename RADIXMAP_T::const_weak_pointer_type node, const std::vector<int> & previous_col, search_context & ctx, const MT & cost_map) {
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
      int current_dist = update_col_linear<MT>(branch[u], ctx.query, current_col, cost_map);
      if(current_dist > ctx.max_distance) {
        max_distance_exceeded = true;
        break;
      }
    }
    if(!max_distance_exceeded) global_search_linear_impl<MT>(ch.second.get(), current_col, ctx, cost_map);
  }
}

TEMPLATE_LIST template <typename MT>
void RADIXMAP_T::anchored_search_linear_impl(typename RADIXMAP_T::const_weak_pointer_type node, const std::vector<int> & previous_col, const int row_min, search_context & ctx, const MT & cost_map) {
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

TEMPLATE_LIST template <typename MT>
int RADIXMAP_T::update_col_affine(const typename RADIXMAP_T::atomic_type branchval, const typename RADIXMAP_T::span_type query, affine_col_type & col, const MT & cost_map) {
  auto & M_col =  std::get<0>(col); // match
  auto & X_col = std::get<1>(col); // gap in query
  auto & Y_col = std::get<2>(col); // gap in target

  int previous_M_i_minus_1 = M_col[0];
  int previous_X_i_minus_1 = X_col[0];
  int previous_Y_i_minus_1 = Y_col[0];

  M_col[0] = NO_ALIGN;
  X_col[0] = previous_X_i_minus_1 == NO_ALIGN ? 
              cost_map.at(pairchar_type(GAP_OPEN_CHAR, branchval)) : 
              previous_X_i_minus_1 + cost_map.at(pairchar_type(GAP_EXTN_CHAR, branchval));
  Y_col[0] = NO_ALIGN;
  int min_element = X_col[0];
  for(size_t i=1; i<M_col.size(); ++i) {
    // Update col[i] after updating {previous_i_minus_1 <- col[i]}
    int M_col_i = cost_map.at(pairchar_type(query[i-1], branchval)) + std::min({
      previous_M_i_minus_1, 
      previous_X_i_minus_1, 
      previous_Y_i_minus_1});
    int X_col_i = std::min({
      cost_map.at(pairchar_type(GAP_OPEN_CHAR, branchval)) + M_col[i],
      cost_map.at(pairchar_type(GAP_EXTN_CHAR, branchval)) + X_col[i],
      cost_map.at(pairchar_type(GAP_OPEN_CHAR, branchval)) + Y_col[i]}); 
    int Y_col_i = std::min({
      cost_map.at(pairchar_type(query[i-1], GAP_OPEN_CHAR)) + M_col[i-1],
      cost_map.at(pairchar_type(query[i-1], GAP_OPEN_CHAR)) + X_col[i-1],
      cost_map.at(pairchar_type(query[i-1], GAP_EXTN_CHAR)) + Y_col[i-1]});
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

TEMPLATE_LIST template <typename MT>
void RADIXMAP_T::global_search_affine_impl(typename RADIXMAP_T::const_weak_pointer_type node, const affine_col_type & previous_col, search_context & ctx, const MT & cost_map) {
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
      int current_dist = update_col_affine<MT>(branch[u], ctx.query, current_col, cost_map);
      if(current_dist > ctx.max_distance) {
        max_distance_exceeded = true;
        break;
      }
    }
    if(!max_distance_exceeded) global_search_affine_impl<MT>(ch.second.get(), current_col, ctx, cost_map);
  }
}

TEMPLATE_LIST template <typename MT>
void RADIXMAP_T::anchored_search_affine_impl(typename RADIXMAP_T::const_weak_pointer_type node, const affine_col_type & previous_col, const int row_min, search_context & ctx, const MT & cost_map) {
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


#undef TEMPLATE_LIST
#undef RADIXMAP_T
#undef GAP_CHAR
#undef GAP_OPEN_CHAR
#undef GAP_EXTN_CHAR
#undef NO_ALIGN
} // end namespace

#endif // include guard