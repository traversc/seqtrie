#ifndef seqtrie_RADIX_TREE_H
#define seqtrie_RADIX_TREE_H

#include "seqtrie/radixmap.h"
#include <stdexcept>
#include <unordered_set>

namespace seqtrie {

class RadixTree {
public:
  using node_type = RadixMap;
  using self_type = RadixMap::self_type;
  using atomic_type = RadixMap::atomic_type;
  using branch_type = RadixMap::branch_type;
  using index_type = RadixMap::index_type;
  using size_type = RadixMap::size_type;
  using pointer_type = RadixMap::pointer_type;
  using const_weak_pointer_type = RadixMap::const_weak_pointer_type;
  using weak_pointer_type = RadixMap::weak_pointer_type;
  using map_type = RadixMap::map_type;
  using span_type = RadixMap::span_type;
  using path = RadixMap::path;
  using search_result = RadixMap::search_result;
  using search_context = RadixMap::search_context;
  using search_context_cached = RadixMap::search_context_cached;
  using NullSearchHook = RadixMap::NullSearchHook;
  using UnitFullSimpleArrayWorkspace = RadixMap::UnitFullSimpleArrayWorkspace;
  using AffineWorkspace = RadixMap::AffineWorkspace;
  using SingleGapWorkspace = RadixMap::SingleGapWorkspace;
  using BandLimits = RadixMap::BandLimits;

  static constexpr index_type nullidx = RadixMap::nullidx;
  static constexpr atomic_type GAP_CHAR = RadixMap::GAP_CHAR;
  static constexpr atomic_type GAP_OPEN_CHAR = RadixMap::GAP_OPEN_CHAR;
  static constexpr atomic_type GAP_EXTN_CHAR = RadixMap::GAP_EXTN_CHAR;
  static constexpr int NO_ALIGN = RadixMap::NO_ALIGN;

  RadixMap root;
  NodePool<RadixMap> node_pool;
  index_type next_terminal_idx = 1;

  RadixTree() = default;
  ~RadixTree() = default;

  RadixTree(const RadixTree &) = delete;
  RadixTree & operator=(const RadixTree &) = delete;
  RadixTree(RadixTree &&) = delete;
  RadixTree & operator=(RadixTree &&) = delete;

  path insert(const span_type sequence, index_type idx) {
    return root.insert(sequence, idx, node_pool);
  }


  index_type erase(const span_type sequence) {
    return root.erase(sequence);
  }

  void compact(bool check_removed_nodes = false) {
    if(check_removed_nodes) {
      check_removed_parent_links();
    }

    NodePool<RadixMap> new_pool;
    new_pool.reserve_exact(node_pool.allocated_nodes());
    root.compact_search_order(new_pool);
    node_pool.swap(new_pool);
  }

private:
  static void collect_reachable(const RadixMap & node,
                                std::unordered_set<const RadixMap *> & reachable) {
    for(const auto ch : node.get_child_nodes()) {
      reachable.insert(ch.second);
      collect_reachable(*ch.second, reachable);
    }
  }

  void check_removed_parent_links() const {
    std::unordered_set<const RadixMap *> reachable;
    reachable.reserve(node_pool.allocated_nodes());
    collect_reachable(root, reachable);

    node_pool.for_each_node([&reachable](const RadixMap & node) {
      const bool is_reachable = reachable.find(&node) != reachable.end();
      if(is_reachable) {
        if(node.get_parent_node() == nullptr) {
          throw std::logic_error("RadixTree::compact found reachable pooled node with null parent");
        }
      } else {
        if(node.get_parent_node() != nullptr) {
          throw std::logic_error("RadixTree::compact found removed pooled node with non-null parent");
        }
      }
    });
  }
};

} // namespace seqtrie

#endif // seqtrie_RADIX_TREE_H
