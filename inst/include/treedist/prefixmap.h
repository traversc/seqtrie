#ifndef TREEDIST_PREFIXMAP_H
#define TREEDIST_PREFIXMAP_H

#include <memory>
#include <algorithm>
#include <limits.h> // INT_MAX
#include "treedist/utility.h"

namespace treedist {

template <template<typename...> class M = std::unordered_map, typename I=uint64_t> class PrefixMap;
template <template<typename...> class M = std::unordered_map, typename I=uint64_t> using PrefixMapUPtr = std::unique_ptr<PrefixMap<M,I>>;
template <template<typename...> class M, typename I>
class PrefixMap {
public:
  typedef I index_type;
  typedef size_t size_type;
  typedef PrefixMap<M,I> value_type;
  typedef PrefixMapUPtr<M,I> pointer_type;
  typedef PrefixMap<M,I>* weak_pointer_type;
  typedef M<uint8_t, PrefixMapUPtr<M,I>> map_type;
  static constexpr I nullidx = std::numeric_limits<index_type>::max();
private:
  map_type child_nodes;        // 56 bytes for std::unordered_map
  index_type terminal_idx;    // 4 bytes
  static inline bool map_exists(const map_type & map, const uint8_t key) {
    return map.find(key) != map.end();
  }
  static inline void map_erase(map_type & map, const uint8_t key) {
    map.erase(key);
  }
  static std::vector<index_type> all_idx(const weak_pointer_type node, size_t max_depth = -1) {
    std::vector<index_type> result;
    if(node->terminal_idx != nullidx) {
      result.push_back(node->terminal_idx);
    }
    if(max_depth == 0) return result;
    for(auto & ch : node->child_nodes) {
      std::vector<index_type> x = all_idx(ch.second.get(), --max_depth);
      result.insert(result.end(), x.begin(), x.end());
    }
    return result;
  }
public:
  PrefixMap() : terminal_idx(nullidx) {}
  static constexpr index_type get_null_idx() { return nullidx; }
  const index_type & get_terminal_idx() const { return terminal_idx; }
  const map_type & get_child_nodes() const { return child_nodes; }
  size_type size() const { return child_nodes.size(); }
  static std::vector<index_type> children(pointer_type & node, size_t max_depth = -1) {
    std::vector<index_type> result;
    if(max_depth == 0) return result;
    if(node->terminal_idx != nullidx) {
      result.push_back(node->terminal_idx);
    }
    for(auto & e : node->child_nodes) {
      std::vector<index_type> & x = children(e.second, --max_depth);
      result.insert(result.end(), x.begin(), x.end());
    }
    return result;
  }
  static std::string print(const pointer_type & node, size_t depth) {
    std::string result;
    if(depth == 0) {
      result += "(root)";
      if(node->terminal_idx != nullidx) {
        result += ":";
        result += std::to_string(node->terminal_idx);
      }
      result += "\n";
    }
    depth++;
    
    std::vector<typename map_type::key_type> child_node_keys;
    for(auto & ch : node->child_nodes) {
      child_node_keys.push_back(ch.first);
    }
    std::sort(child_node_keys.begin(), child_node_keys.end()); // sort for reproducible printing
    for(auto k : child_node_keys) {
      for(size_t i=0; i<depth; ++i) result += " ";
      result += static_cast<char>(k);
      if(node->child_nodes[k]->terminal_idx != nullidx) {
        result += ":";
        result += std::to_string(node->terminal_idx);
      }
      result += "\n";
      result += print(node->child_nodes[k], depth);
    }
    return result;
  }
  static index_type find(const pointer_type & root, const uspan sequence) {
    value_type * node = root.get();
    size_t position=0;
    while(position < sequence.size()) {
      if(map_exists(node->child_nodes, sequence[position])) {
        node = node->child_nodes[sequence[position++]].get();
      } else {
        return nullidx;
      }
    }
    return node->terminal_idx;
  }
  static std::vector<index_type> find_prefix(const pointer_type & root, const uspan sequence) {
    value_type * node = root.get();
    size_t position=0;
    while(position < sequence.size()) {
      if(map_exists(node->child_nodes, sequence[position])) {
        node = node->child_nodes[sequence[position++]].get();
      } else {
        return std::vector<index_type>(0);
      }
    }
    return all_idx(node);
  }
  static index_type insert(pointer_type & root, const uspan sequence, index_type idx) {
    value_type * node = root.get();
    for(size_t i=0; i<sequence.size(); ++i) {
      if(! map_exists(node->child_nodes, sequence[i])) {
        node->child_nodes[sequence[i]] = pointer_type(new value_type);
      }
      node = node->child_nodes[sequence[i]].get();
    }
    if(node->terminal_idx == nullidx) {
      node->terminal_idx = idx;
      return nullidx;
    } else {
      return node->terminal_idx;
    }
  }
static index_type erase(pointer_type & node, const uspan sequence) {
  index_type result = nullidx;
  erase(node, sequence, result);
  return result;
}
private:
  enum class erase_action { erase, keep };
  static erase_action erase(pointer_type & node, const uspan sequence, index_type & result) {
    // end of sequence; could be internal branch
    if(sequence.size() == 0) {
      std::swap(result, node->terminal_idx); // if sequence doesn't exist, terminal_idx should be nullidx which is fine since result is initialized as nullidx
      size_t nc = node->size();
      if(nc == 0) { // no children
        return erase_action::erase;
      } else {
        return erase_action::keep;
      }
    }
    
    // check that sequence actually exists in tree - we shouldn't assume it does
    size_t s = sequence[0];
    if(!map_exists(node->child_nodes, s)) {
      return erase_action::keep;
    }
    
    // travelling down tree
    erase_action action = erase(node->child_nodes[s], subspan(sequence,1), result);
    
    // travelling back up
    if(action == erase_action::keep) {
      return erase_action::keep;
    } else { // if(action == erase_action::erase) {
      map_erase(node->child_nodes, s);
      size_t nc = node->size();
      if((nc == 0) && (node->terminal_idx == nullidx)) { // no children and not a sequence
        return erase_action::erase;
      } else {
        return erase_action::keep;
      }
    }
  }
////////////////////////////////////////////////////////////////////////////////
// type T is a vector type, e.g. std::vector or tbb:concurrent_vector
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
      for (auto & ch : node->child_nodes) {
        search(ch.second, update_row(ch.first, previous_row));
      }
    }
    inline std::vector<int> update_row(const uint8_t branchval, const std::vector<int> & previous_row) {
      std::vector<int> current_row(previous_row.size());
      current_row[0] = previous_row[0] + 1;
      for(size_t i=1; i<current_row.size(); ++i) {
        int match_cost  = previous_row[i-1] + (sequence[i-1] == branchval ? 0 : 1);
        int insert_cost = current_row[i-1] + 1;
        int delete_cost = previous_row[i] + 1;
        current_row[i] = std::min({match_cost, insert_cost, delete_cost});
      }
      return current_row;
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
        for (auto & x : node->child_nodes) {
          int new_distance = x.first == sequence[position] ? distance : distance + 1;
          size_t new_position = position + 1;
          if(new_distance > max_distance) continue;
          if(new_position == sequence.size()) {
            if(x.second->terminal_idx != nullidx) {
              output_index.push_back(x.second->terminal_idx);
              output_distance.push_back(new_distance);
            }
          } else {
            search(x.second, new_position, new_distance);
          }
        }
      }
    }
  };
};


}

#endif // include guard