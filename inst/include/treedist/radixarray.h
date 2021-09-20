#ifndef TREEDIST_RADIXARRAY_H
#define TREEDIST_RADIXARRAY_H

#include "treedist/utility.h"

namespace treedist {

/* interface (to do: put in virtual class?)
private:
 static inline bool map_exists(const map_type & map, const uint8_t key);
 static inline void map_erase(map_type & map, const uint8_t key);
 static std::vector<index_type> all_idx(weak_pointer_type node, size_t max_depth = -1);
public:
 TRIECLASS() : terminal_idx(nullidx) {}
 static constexpr index_type get_null_idx() { return nullidx; }
 const branch_type & get_branch() const { return branch; }
 const index_type & get_terminal_idx() const { return terminal_idx; }
 const map_type & get_child_nodes() const { return child_nodes; }
 size_type size() const;
 static std::string print(const pointer_type & node, size_t depth);
 static std::tuple<std::vector<index_type> std::vector<index_type>, std::vector<uspan>> graph(const pointer_type & node, size_t depth);
 static index_type find(const pointer_type & root, const uspan sequence);
 static std::vector<index_type> find_prefix(const pointer_type & root, const uspan sequence);
 static index_type insert(pointer_type & node, const uspan sequence, index_type idx);
 static index_type erase(pointer_type & node, const uspan sequence);
 template <typename T> static T convert_sequence(const cspan x);
 template <typename T> static T recover_sequence(const cspan x)
public:
 Levenshtein
public:
 Levenshtein(const pointer_type & root, const uspan sequence, const int max_distance);
 result_type search();
 Hamming
public:
 Hamming(const pointer_type & root, const uspan sequence, const int max_distance);
 result_type search();
 */

// templates: A = Alphabet, B = branch vector type e.g. std::vector, I = index type e.g. size_t
template <class A = boost::mpl::string<'A','C','G','T'>, template<typename...> class B = std::vector, class I=size_t> class RadixArray;
template <class A = boost::mpl::string<'A','C','G','T'>, template<typename...> class B = std::vector, class I=size_t> using RadixArrayUPtr = std::unique_ptr<RadixArray<A,B,I>>;
template <class A, template<typename...> class B, class I>
class RadixArray {
public:
  static constexpr I nullidx = std::numeric_limits<I>::max();
  static constexpr size_t degree = boost::mpl::size<A>::value;
  static constexpr const char * charset_begin = &boost::mpl::c_str<A>::value[0];
  static constexpr const char * charset_end = &boost::mpl::c_str<A>::value[0] + degree;
  typedef A charset_type;
  typedef B<uint8_t> branch_type;
  typedef I index_type;
  typedef size_t size_type;
  typedef RadixArray<A,B,I> value_type;
  typedef RadixArrayUPtr<A,B,I> pointer_type;
  typedef RadixArray<A,B,I>* weak_pointer_type;
  typedef std::array<RadixArrayUPtr<A,B,I>,degree> map_type;
  typedef B<uint8_t> convert_sequence_type;
  typedef B<char> recover_sequence_type;
private:
  map_type child_nodes;      // 32 bytes
  branch_type branch;         // 24 bytes for std::vector
  index_type terminal_idx;    // 4 bytes
  static inline bool map_exists(const map_type & map, const uint8_t key) {
    return map[key].get() != nullptr;
  }
  static inline void map_erase(map_type & map, const uint8_t key) {
    map[key].reset();
  }
  static std::vector<index_type> all_idx(const weak_pointer_type node, size_t max_depth = -1) {
    std::vector<index_type> result;
    if(node->terminal_idx != nullidx) {
      result.push_back(node->terminal_idx);
    }
    if(max_depth == 0) return result;
    for(size_t i=0; i<degree; ++i) {
      if(map_exists(node->child_nodes, i)) {
        std::vector<index_type> x = all_idx(node->child_nodes[i].get(), --max_depth);
        result.insert(result.end(), x.begin(), x.end());
      }
    }
    return result;
  }
public:
  RadixArray() : terminal_idx(nullidx) {}
  static constexpr index_type get_null_idx() { return nullidx; }
  const branch_type & get_branch() const { return branch; }
  const index_type & get_terminal_idx() const { return terminal_idx; }
  const map_type & get_child_nodes() const { return child_nodes; }
  size_type size() const { // number of direct children of current node
    size_t result = 0;
    for(size_t i=0; i<degree; ++i) {
      if(child_nodes[i]) ++result;
    }
    return result;
  }
  static convert_sequence_type convert_sequence(const cspan x) {
    convert_sequence_type result(x.size());
    for(size_t i=0; i<x.size(); ++i) {
      auto it = std::find(charset_begin, charset_end, x[i]);
      if(it == charset_end) {
        static const std::string err_msg = "Sequence must be " + std::string(charset_begin);
        throw std::runtime_error(err_msg);
      }
      result[i] = std::distance(charset_begin, it);
    }
    return result;
  }
  static recover_sequence_type recover_sequence(const uspan x) {
    recover_sequence_type result(x.size());
    for(size_t i=0; i<x.size(); ++i) {
      result[i] = charset_begin[x[i]];
    }
    return result;
  }
  static std::string print(const pointer_type & node, size_t depth = 0) {
    std::string result;
    if(depth == 0) {
      result += "(root)";
    } else {
      for(size_t i=0; i<depth; ++i) result += " ";
    }
    std::string x(node->branch.size(),0);
    for(size_t j=0; j<node->branch.size(); ++j) {
      x[j] = charset_begin[node->branch[j]];
    }
    result += x;
    
    if(node->terminal_idx != nullidx) {
      result += ":";
      result += std::to_string(node->terminal_idx);
    }
    result += "\n";
    for(uint8_t i=0; i<degree; ++i) {
      if(map_exists(node->child_nodes, i)) {
        result += print(node->child_nodes[i], depth + 1);
      }
    }
    if(depth == 0) {
      result += "\n";
    }
    return result;
  }
//   struct graph_results {
//     std::vector<size_t> parent;
//     std::vector<size_t> child;
//     std::vector<uspan> branch;
//   };
//   graph_container graph(const pointer_type & node, size_t depth = -1); {
//     std::vector<weak_pointer_type> parent;
//     std::vector<weak_pointer_type> child;
//     graph(node, depth, parent, child);
//     std::unordered_map<weak_pointer_type, size_t> node_map;
//     graph_results result;
//     for(size_t i=0; i<parent.size(); ++i) {
//       if(node_map.find(parent[i]) == node_map.end()) {
//         result.parent.push_back(result.branch.size());
//         node_map.emplace(parent[i], result.branch.size());
//         result.branch.push_back(parent[i]->branch);
//       } else {
//         result.parent.push_back(node_map[parent[i]]);
//       }
//       if(node_map.find(child[i]) == node_map.end()) {
//         result.child.push_back(result.branch.size());
//         node_map.emplace(child[i], result.branch.size());
//         result.branch.push_back(child[i]->branch);
//       } else {
//         result.child.push_back(node_map[child[i]]);
//       }
//     }
//     return result;
//   }
// private:
//   void graph(const pointer_type & node, size_t depth, std::vector<weak_pointer_type> & parent, std::vector<weak_pointer_type> & child) {
//     if(depth == 0) return;
//     for(size_t i=0; i<degree; ++i) {
//       if(map_exists(node->child_nodes, i)) {
//         parent.push_back(node->get());
//         child.push_back(node->child_nodes[i].get());
//         graph(node->child_nodes[i], depth - 1, parent, child);
//       }
//     }
//   }
// public:
  static index_type find(const pointer_type & root, const uspan sequence) {
    value_type * node = root.get();
    size_t position=0;
    while(position < sequence.size()) {
      if(map_exists(node->child_nodes, sequence[position])) {
        node = node->child_nodes[sequence[position]].get();
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
  static std::vector<index_type> find_prefix(const pointer_type & root, const uspan sequence) {
    value_type * node = root.get();
    size_t sequence_position = 0;
    size_t branch_position = 0;
    while(sequence_position < sequence.size()) {
      if(branch_position >= node->branch.size()) {
        if(map_exists(node->child_nodes, sequence[sequence_position])) {
          node = node->child_nodes[sequence[sequence_position]].get();
          branch_position = 0;
        } else {
          return std::vector<index_type>(0);
        }
      }
      if(node->branch[branch_position] == sequence[sequence_position]) {
        branch_position++;
        sequence_position++;
      } else {
        return std::vector<index_type>(0);
      }
    }
    return all_idx(node);
  }
  static index_type insert(pointer_type & node, const uspan sequence, index_type idx) {
    if(sequence.size() == 0) {
      // std::cout << "case -1: end of sequence" << std::endl;
      if(node->terminal_idx == nullidx) {
        node->terminal_idx = idx;
        return nullidx;
      } else {
        return node->terminal_idx;
      }
    }
    uint8_t s = sequence[0];
    if(!map_exists(node->child_nodes,s)) {
      // std::cout << "case 0: no existing branch" << std::endl;
      node->child_nodes[s] = pointer_type(new value_type);
      node->child_nodes[s]->branch = subvector<branch_type>(sequence, 0);
      node->child_nodes[s]->terminal_idx = idx;
      return nullidx;
    }
    size_t i = 0;
    while(i < node->child_nodes[s]->branch.size() && i < sequence.size() && sequence[i] == node->child_nodes[s]->branch[i]) { ++i; }
    
    if(i == sequence.size() && i == node->child_nodes[s]->branch.size()) {
      // std::cout << "case 1: sequence is same as branch" << std::endl;
      if(node->child_nodes[s]->terminal_idx == nullidx) {
        node->child_nodes[s]->terminal_idx = idx;
        return nullidx;
      } else {
        return node->child_nodes[s]->terminal_idx;
      }
    } else if(i == sequence.size()) {
      // std::cout << "case 2: sequence is shorter than branch" << std::endl;
      branch_type branch_prefix = subvector<branch_type>(node->child_nodes[s]->branch,0,i);
      branch_type branch_suffix = subvector<branch_type>(node->child_nodes[s]->branch,i);
      uint8_t s_insert = branch_suffix[0];
      pointer_type inserted_node(new value_type);
      inserted_node->child_nodes[s_insert] = std::move(node->child_nodes[s]);
      inserted_node->branch = std::move(branch_prefix);
      inserted_node->terminal_idx = idx;
      inserted_node->child_nodes[s_insert]->branch = std::move(branch_suffix);
      node->child_nodes[s] = std::move(inserted_node);
      return nullidx;
    } else if(i == node->child_nodes[s]->branch.size()) {
      // std::cout << "case 3: sequence is longer than branch" << std::endl;
      uspan seq_suffix = subspan(sequence, i); // will remain as ustring_view
      return insert(node->child_nodes[s], seq_suffix, idx);
    } else {
      // std::cout << "case 4: sequence is different than branch" << std::endl;
      branch_type branch_prefix = subvector<branch_type>(node->child_nodes[s]->branch,0,i);
      branch_type branch_suffix = subvector<branch_type>(node->child_nodes[s]->branch,i);
      branch_type seq_suffix = subvector<branch_type>(sequence,i);
      uint8_t s_insert_branch = branch_suffix[0];
      uint8_t s_insert_seq = seq_suffix[0];
      
      pointer_type inserted_node(new value_type);
      inserted_node->child_nodes[s_insert_branch] = std::move(node->child_nodes[s]);
      inserted_node->child_nodes[s_insert_branch]->branch = std::move(branch_suffix);
      inserted_node->child_nodes[s_insert_seq] = pointer_type(new value_type);
      inserted_node->child_nodes[s_insert_seq]->branch = std::move(seq_suffix);
      inserted_node->child_nodes[s_insert_seq]->terminal_idx = idx;
      node->child_nodes[s] = std::move(inserted_node);
      node->child_nodes[s]->branch = std::move(branch_prefix);
      return nullidx;
    }
  }
  static index_type erase(pointer_type & node, const uspan sequence) {
    index_type result = nullidx;
    erase(node, sequence, result);
    return result;
  }
private:
  enum class erase_action { erase, merge, keep };
  static erase_action erase(pointer_type & node, const uspan sequence, index_type & result) {
    // end of sequence; could be internal branch
    if(sequence.size() == 0) {
      std::swap(result, node->terminal_idx); // if sequence doesn't exist, terminal_idx should be nullidx which is fine since result is initialized as nullidx
      size_t nc = node->size();
      if(nc == 0) { // no children
        return erase_action::erase;
      } else if(nc == 1) { // one child
        return erase_action::merge;
      } else { // two or more
        return erase_action::keep;
      }
    }
    
    // check that sequence actually exists in tree - we shouldn't assume it does
    uint8_t s = sequence[0];
    if(!map_exists(node->child_nodes, s)) {
      return erase_action::keep;
    }
    
    // travelling down tree
    size_t i = 0;
    for(; i<node->child_nodes[s]->branch.size(); ++i) {
      if(i == sequence.size()) return erase_action::keep; // branch is longer than sequence, doesn't match
      if(node->child_nodes[s]->branch[i] != sequence[i]) return erase_action::keep; // branch and sequence don't match
    }
    erase_action action = erase(node->child_nodes[s], subspan(sequence,i), result); // sequence is longer or same
    
    // travelling back up
    if(action == erase_action::keep) {
      return erase_action::keep;
    } else if(action == erase_action::merge) {
      uint8_t next_s = 0;
      for(; next_s < degree; ++next_s) {
        if(map_exists(node->child_nodes[s]->child_nodes, next_s)) break;
      }
      appendspan(node->child_nodes[s]->branch, node->child_nodes[s]->child_nodes[next_s]->branch);
      branch_type next_branch = std::move(node->child_nodes[s]->branch);
      node->child_nodes[s] = std::move(node->child_nodes[s]->child_nodes[next_s]);
      node->child_nodes[s]->branch = std::move(next_branch);
      return erase_action::keep;
    } else { // if(action == erase_action::erase) {
      map_erase(node->child_nodes, s);
      size_t nc = node->size();
      if((nc == 0) && (node->terminal_idx == nullidx)) { // no children and not a sequence
        return erase_action::erase;
      } else if((nc == 1) && (node->terminal_idx == nullidx)) { // one child and not a sequence
        return erase_action::merge;
      } else { // 2+ children or node is a sequence
        return erase_action::keep;
      }
    }
  }
////////////////////////////////////////////////////////////////////////////////
// levenshtein; type T is a vector type, e.g. std::vector or tbb:concurrent_vector
public:
  class Levenshtein {
  public:
    typedef std::pair<std::vector<index_type>, std::vector<int>> result_type;
  private:
    const pointer_type & root;
    std::vector<index_type> output_index;
    std::vector<int> output_distance;
    const uspan sequence;
    const int max_distance;
  public:
    Levenshtein(const pointer_type & root, const uspan sequence, const int max_distance) : root(root), sequence(sequence), max_distance(max_distance) {}
    result_type search() {
      std::vector<int> start_row(sequence.size() + 1);
      std::iota(start_row.begin(), start_row.end(), 0);
      search(root, start_row);
      return std::make_pair(std::move(output_index), std::move(output_distance));
    }
  private:
    void search(const pointer_type & node, const std::vector<int> & previous_row) {
      if( *std::min_element(previous_row.begin(), previous_row.end()) > max_distance ) {
        return;
      }
      if((node->terminal_idx != nullidx) && (previous_row.back() <= max_distance)) {
        output_index.push_back(node->terminal_idx);
        output_distance.push_back(previous_row.back());
      }
      for (size_t i=0; i<degree; ++i) {
        if(! map_exists(node->child_nodes, i)) continue;
        std::vector<int> current_row = previous_row;
        branch_type & branch = node->child_nodes[i]->branch;
        for(size_t u=0; u<branch.size(); ++u) {
          int current_dist = update_row(branch[u], current_row);
          if(current_dist > max_distance) goto continue_outer_loop;
        }
        search(node->child_nodes[i], current_row);
        continue_outer_loop:;
      }
    }
    inline int update_row(const uint8_t branchval, std::vector<int> & row) {
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
  };
////////////////////////////////////////////////////////////////////////////////
// hamming; same interface as levenshtein
  class Hamming {
  public:
    typedef std::pair<std::vector<index_type>, std::vector<int>> result_type;
  private:
    const pointer_type & root;
    std::vector<index_type> output_index;
    std::vector<int> output_distance;
    const uspan sequence;
    const int max_distance;
  public:
    Hamming(const pointer_type & root, const uspan sequence, const int max_distance) : root(root), sequence(sequence), max_distance(max_distance) {}
    result_type search() {
      if(sequence.size() == 0) { // corner case
        if(root->terminal_idx != nullidx) {
          output_index.push_back(root->terminal_idx);
          output_distance.push_back(0);
        }
      } else {
        search(root, 0, 0);
      }
      return std::make_pair(std::move(output_index), std::move(output_distance));
    }
  private:
    void search(const pointer_type & node, const size_t position, const int distance) {
      if(position < sequence.size()) {
        for (size_t i=0; i<degree; ++i) {
          if(! map_exists(node->child_nodes, i)) continue;
          branch_type & branch = node->child_nodes[i]->branch;
          int new_distance = distance;
          if(position + branch.size() > sequence.size()) continue;
          for(size_t j=0; j<branch.size(); ++j) {
            if(branch[j] != sequence[position+j]) new_distance++;
            if(new_distance > max_distance) goto continue_outer_loop; // continue out of outer for loop
          }
          size_t new_position; new_position = position + branch.size(); // goto issue https://stackoverflow.com/q/11306799/2723734 
          if(new_position == sequence.size()) {
            if(node->child_nodes[i]->terminal_idx != nullidx) {
              output_index.push_back(node->child_nodes[i]->terminal_idx);
              output_distance.push_back(new_distance);
            }
          } else {
            search(node->child_nodes[i], new_position, new_distance);
          }
          continue_outer_loop:;
        }
      }
    }
  };
};

}


#endif // include guard
