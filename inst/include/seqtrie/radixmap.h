#ifndef seqtrie_RADIXMAP_H
#define seqtrie_RADIXMAP_H

#include "seqtrie/utility.h"

namespace seqtrie {

template <class A = char, template<typename...> class M = std::map, template<typename...> class B = std::vector, class I=size_t> class RadixMap;
template <class A = char, template<typename...> class M = std::map, template<typename...> class B = std::vector, class I=size_t> using RadixMapUPtr = std::unique_ptr<RadixMap<A,M,B,I>>;
template <class A, template<typename...> class M, template<typename...> class B, class I>
class RadixMap {
public:
  static constexpr I nullidx = std::numeric_limits<I>::max();
  typedef A atomic_type;
  typedef B<A> branch_type;
  typedef I index_type;
  typedef size_t size_type;
  typedef RadixMap<A,M,B,I> self_type;
  typedef RadixMapUPtr<A,M,B,I> pointer_type;
  typedef RadixMap<A,M,B,I> const * const_weak_pointer_type;
  typedef RadixMap<A,M,B,I> * weak_pointer_type;
  typedef M<A,RadixMapUPtr<A,M,B,I>> map_type;
  typedef nonstd::span<const A> span_type;
private:
  map_type child_nodes;        // 48 bytes for std::map
  branch_type branch;         // 24 bytes for std::vector
  const_weak_pointer_type parent_node; // 8 bytes, worth including?
  index_type terminal_idx;    // 8 bytes
public:
  RadixMap() : parent_node(nullptr), terminal_idx(nullidx) {}
  struct path { // like an iterator, but without iteration. Later we could turn this into a real iterator. 
    const_weak_pointer_type m;
    path() : m(nullptr) {}
    path(const_weak_pointer_type x) : m(x) {}
    template <typename T> T sequence() const;
    std::vector<atomic_type> operator*() const { return sequence<std::vector<atomic_type>>(); }
    index_type index() const { return m->terminal_idx; }
  };
  struct search_context { // result struct for levenshtein and hamming search
    std::vector<path> match;
    std::vector<int> distance;
    span_type sequence;
    int max_distance;
    search_context() {}
    search_context(span_type sequence, int max_distance) : sequence(sequence), max_distance(max_distance) {}
  };
  size_type size() const;
  bool validate(const bool is_root = true) const; // check if tree makes sense
  std::vector<path> all(size_t max_depth = -1) const;
  std::string print() const { return print_impl(0); } // this function really only makes sense if atomic_type is char...
  std::pair<std::vector<path>, std::vector<path>> graph(size_t max_depth = -1) const;
  index_type find(const span_type sequence) const;
  index_type insert(const span_type sequence, index_type idx);
  index_type erase(const span_type sequence);
  std::vector<path> prefix_search(const span_type sequence) const;
  search_context levenshtein_search(const span_type sequence, const int max_distance) const;
  search_context hamming_search(const span_type sequence, const int max_distance) const;
  
private:
  // implementation helpers, subject to change
  std::string print_impl(size_t depth) const;
  enum class erase_action { erase, merge, keep };
  static erase_action erase_impl(weak_pointer_type node, const span_type sequence, index_type & result);
  static int update_row(const atomic_type branchval, const span_type sequence, std::vector<int> & row); // helper for levenshtein search
  static void levenshtein_search_impl(const_weak_pointer_type node, const std::vector<int> & previous_row, search_context & ctx);
  static void hamming_search_impl(const_weak_pointer_type node, const size_t position, const int distance, search_context & ctx);
};


#define TEMPLATE_LIST template <class A, template<typename...> class M, template<typename...> class B, class I>
#define RADIXMAP_T RadixMap<A,M,B,I>

TEMPLATE_LIST template <typename T>
T RADIXMAP_T::path::sequence() const {
  static_assert(std::is_same<typename T::value_type, atomic_type>::value, "Output sequence value_type must be the same as RadixMap::atomic_type.");
  const_weak_pointer_type current = m;
  std::vector<const_weak_pointer_type> h;
  size_t size = 0;
  while(current != nullptr) {
    h.push_back(current);
    size += current->branch.size();
    current = current->parent_node;
  }
  T result = array_allocate<T>(size);
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
  if(terminal_idx != nullidx) {
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

TEMPLATE_LIST typename RADIXMAP_T::index_type RADIXMAP_T::find(const typename RADIXMAP_T::span_type sequence) const {
  const_weak_pointer_type node = this;
  size_t position=0;
  while(position < sequence.size()) {
    if(node->child_nodes.find(sequence[position]) != node->child_nodes.end()) {
      node = node->child_nodes.at(sequence[position]).get();
      if(position + node->branch.size() > sequence.size()) return nullidx;
      for(size_t j=0; j<node->branch.size(); ++j) {
        if(node->branch[j] != sequence[position+j]) return nullidx;
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

TEMPLATE_LIST std::vector<typename RADIXMAP_T::path> RADIXMAP_T::prefix_search(const typename RADIXMAP_T::span_type sequence) const {
  const_weak_pointer_type node = this;
  size_t sequence_position = 0;
  size_t branch_position = 0;
  while(sequence_position < sequence.size()) {
    if(branch_position >= node->branch.size()) {
      if(node->child_nodes.find(sequence[sequence_position]) != node->child_nodes.end()) {
        node = node->child_nodes.at(sequence[sequence_position]).get();
        branch_position = 0;
      } else {
        return std::vector<path>(0);
      }
    }
    if(node->branch[branch_position] == sequence[sequence_position]) {
      branch_position++;
      sequence_position++;
    } else {
      return std::vector<path>(0);
    }
  }
  return node->all();
}

TEMPLATE_LIST typename RADIXMAP_T::search_context RADIXMAP_T::levenshtein_search(const typename RADIXMAP_T::span_type sequence, const int max_distance) const {
  search_context ctx(sequence, max_distance);
  levenshtein_search_impl(this, iota_range<std::vector<int>>(0, sequence.size() + 1), ctx); 
  return ctx;
}

TEMPLATE_LIST typename RADIXMAP_T::search_context RADIXMAP_T::hamming_search(const span_type sequence, const int max_distance) const {
  search_context ctx(sequence, max_distance);
  hamming_search_impl(this, 0, 0, ctx);
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

TEMPLATE_LIST int RADIXMAP_T::update_row(const typename RADIXMAP_T::atomic_type branchval, const typename RADIXMAP_T::span_type sequence, std::vector<int> & row) {
  int previous_row_i_minus_1 = row[0];
  row[0] = row[0] + 1;
  int min_element = row[0];
  for(size_t i=1; i<row.size(); ++i) {
    int match_cost  = previous_row_i_minus_1 + (sequence[i-1] == branchval ? 0 : 1);
    int insert_cost = row[i-1] + 1;
    int delete_cost = row[i] + 1;
    previous_row_i_minus_1 = row[i];
    row[i] = std::min({match_cost, insert_cost, delete_cost});
    if(row[i] < min_element) min_element = row[i];
  }
  return min_element;
}
TEMPLATE_LIST void RADIXMAP_T::levenshtein_search_impl(typename RADIXMAP_T::const_weak_pointer_type node, const std::vector<int> & previous_row, search_context & ctx) {
  if( *std::min_element(previous_row.begin(), previous_row.end()) > ctx.max_distance ) { return; }
  if((node->terminal_idx != nullidx) && (previous_row.back() <= ctx.max_distance)) {
    ctx.match.push_back(path(node));
    ctx.distance.push_back(previous_row.back());
  }
  for (auto & ch : node->child_nodes) {
    std::vector<int> current_row = previous_row;
    branch_type & branch = ch.second->branch;
    bool max_distance_exceeded = false;
    for(size_t u=0; u<branch.size(); ++u) {
      int current_dist = update_row(branch[u], ctx.sequence, current_row);
      if(current_dist > ctx.max_distance) {
        max_distance_exceeded = true;
        break;
      }
    }
    if(!max_distance_exceeded) levenshtein_search_impl(ch.second.get(), current_row, ctx);
  }
}

TEMPLATE_LIST void RADIXMAP_T::hamming_search_impl(typename RADIXMAP_T::const_weak_pointer_type node, const size_t position, const int distance, typename RADIXMAP_T::search_context & ctx) {
  if(position == ctx.sequence.size()) {
    if(node->terminal_idx != nullidx) {
      ctx.match.push_back(path(node));
      ctx.distance.push_back(distance);
    }
    return;
  }
  if(position < ctx.sequence.size()) {
    for (auto & ch : node->child_nodes) {
      branch_type & branch = ch.second->branch;
      int new_distance = distance;
      if(position + branch.size() > ctx.sequence.size()) continue;
      bool max_distance_exceeded = false;
      for(size_t j=0; j<branch.size(); ++j) {
        if(branch[j] != ctx.sequence[position+j]) new_distance++;
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


#undef TEMPLATE_LIST
#undef RADIXMAP_T
} // end namespace

#endif // include guard
