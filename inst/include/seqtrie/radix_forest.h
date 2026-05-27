#ifndef seqtrie_RADIX_FOREST_H
#define seqtrie_RADIX_FOREST_H

#include "ankerl/unordered_dense.h"
#include "seqtrie/radix_tree.h"
#include <memory>

namespace seqtrie {

class RadixForest {
public:
  using tree_type = RadixTree;
  using tree_map_type = ankerl::unordered_dense::map<size_t, std::unique_ptr<RadixTree>>;
  using span_type = RadixTree::span_type;
  using path = RadixTree::path;
  using index_type = RadixTree::index_type;
  using search_result = RadixMap::search_result;
  using search_context = RadixMap::search_context;
  using search_context_cached = RadixMap::search_context_cached;
  using NullSearchHook = RadixMap::NullSearchHook;
  using UnitFullSimpleArrayWorkspace = RadixMap::UnitFullSimpleArrayWorkspace;
  using AffineWorkspace = RadixMap::AffineWorkspace;
  using CostMap = seqtrie::CostMap;
  using QueryCostCache = seqtrie::QueryCostCache;

  tree_map_type tree_map;
  index_type next_terminal_idx = 1;

  RadixForest() = default;
  ~RadixForest() = default;

  RadixForest(const RadixForest &) = delete;
  RadixForest & operator=(const RadixForest &) = delete;
  RadixForest(RadixForest &&) = delete;
  RadixForest & operator=(RadixForest &&) = delete;

  path insert(const span_type sequence, index_type idx) {
    auto it = tree_map.find(sequence.size());
    if(it == tree_map.end()) {
      it = tree_map.emplace(sequence.size(), std::make_unique<RadixTree>()).first;
    }
    return it->second->insert(sequence, idx);
  }


  index_type erase(const span_type sequence) {
    auto it = tree_map.find(sequence.size());
    if(it == tree_map.end()) return RadixMap::nullidx;

    index_type result = it->second->erase(sequence);
    const auto & root = it->second->root;
    if(root.get_terminal_idx() == RadixMap::nullidx && root.get_child_nodes().empty()) {
      tree_map.erase(it);
    }
    return result;
  }

  void compact(bool check_removed_nodes = false) {
    for(auto & x : tree_map) {
      x.second->compact(check_removed_nodes);
    }
  }

  search_result hamming_search(const span_type query, int max_distance) const {
    NullSearchHook hook;
    return hamming_search(query, max_distance, hook);
  }

  template <typename Hook>
  search_result hamming_search(const span_type query, int max_distance, Hook & hook) const {
    auto it = tree_map.find(query.size());
    if(it == tree_map.end()) return search_result{};
    return it->second->root.hamming_search(query, max_distance, hook).release_result();
  }

  search_result global_search(const span_type query, int max_distance) const {
    NullSearchHook hook;
    return global_search(query, max_distance, hook);
  }

  template <typename Hook>
  search_result global_search(const span_type query, int max_distance, Hook & hook) const {
    const size_t radius = static_cast<size_t>(max_distance);
    const size_t len = query.size();
    const size_t min_search_len = len > radius ? len - radius : 0;
    const size_t max_search_len = len + radius;

    search_context ctx(query, max_distance);
    UnitFullSimpleArrayWorkspace workspace;
    workspace.reset_global(query.size(), max_search_len + 1);

    for(size_t j = min_search_len; j <= max_search_len; ++j) {
      auto it = tree_map.find(j);
      if(it != tree_map.end()) {
        if constexpr (std::is_same_v<Hook, NullSearchHook>) {
          RadixMap::global_search_impl<Hook>(&it->second->root, 0, 0, ctx, workspace, hook, 0);
        } else {
          if(!RadixMap::global_search_impl<Hook>(&it->second->root, 0, 0, ctx, workspace, hook, 0)) break;
        }
      }
    }

    return std::move(ctx).release_result();
  }

  search_result global_search_linear(const span_type query, int max_distance, const CostMap & cost_map) const {
    NullSearchHook hook;
    return global_search_linear(query, max_distance, cost_map, hook);
  }

  template <typename Hook>
  search_result global_search_linear(const span_type query, int max_distance, const CostMap & cost_map, Hook & hook) const {
    QueryCostCache query_costs(cost_map, query);
    UnitFullSimpleArrayWorkspace workspace;

    const size_t len = query_costs.query_len;
    const size_t radius = RadixMap::BandLimits::linear_radius(max_distance, query_costs.gap_cost, len);
    search_context_cached ctx(query_costs, max_distance, radius);
    const size_t min_search_len = len > radius ? len - radius : 0;
    const size_t max_search_len = len + radius;

    workspace.reset_global(len, max_search_len + 1);
    int * root_col = workspace.at(0);
    for(size_t i = 1; i < workspace.column_len; ++i) {
      root_col[i] = root_col[i - 1] + query_costs.gap_cost; // gap in target
    }

    for(size_t j = min_search_len; j <= max_search_len; ++j) {
      auto it = tree_map.find(j);
      if(it != tree_map.end()) {
        if constexpr (std::is_same_v<Hook, NullSearchHook>) {
          RadixMap::global_search_linear_impl<Hook>(&it->second->root, 0, 0, ctx, workspace, hook, 0);
        } else {
          if(!RadixMap::global_search_linear_impl<Hook>(&it->second->root, 0, 0, ctx, workspace, hook, 0)) break;
        }
      }
    }

    return std::move(ctx).release_result();
  }

  search_result global_search_affine(const span_type query, int max_distance, const CostMap & cost_map) const {
    NullSearchHook hook;
    return global_search_affine(query, max_distance, cost_map, hook);
  }

  template <typename Hook>
  search_result global_search_affine(const span_type query, int max_distance, const CostMap & cost_map, Hook & hook) const {
    QueryCostCache query_costs(cost_map, query);

    const size_t len = query_costs.query_len;
    const size_t radius = RadixMap::BandLimits::affine_radius(max_distance,
                                                              query_costs.gap_cost,
                                                              query_costs.gap_open_including_first_extension,
                                                              len);
    search_context_cached ctx(query_costs, max_distance, radius);
    const size_t min_search_len = len > radius ? len - radius : 0;
    const size_t max_search_len = len + radius;

    const size_t col_size = len + 1;
    AffineWorkspace workspace;
    workspace.initialize(len, max_search_len + 1);
    RadixMap::affine_col_type & col = workspace.at(0);
    auto & M_col = std::get<0>(col);
    auto & X_col = std::get<1>(col);
    auto & Y_col = std::get<2>(col);
    M_col[0] = 0;
    X_col[0] = RadixMap::NO_ALIGN;
    Y_col[0] = RadixMap::NO_ALIGN;
    for(size_t i = 1; i < col_size; ++i) {
      M_col[i] = RadixMap::NO_ALIGN;
      X_col[i] = RadixMap::NO_ALIGN;
      if(i == 1) Y_col[i] = query_costs.gap_open_including_first_extension;
      else       Y_col[i] = Y_col[i - 1] + query_costs.gap_cost;
    }

    for(size_t j = min_search_len; j <= max_search_len; ++j) {
      auto it = tree_map.find(j);
      if(it != tree_map.end()) {
        if constexpr (std::is_same_v<Hook, NullSearchHook>) {
          RadixMap::global_search_affine_impl<Hook>(&it->second->root, 0, 0, ctx, workspace, hook, 0);
        } else {
          if(!RadixMap::global_search_affine_impl<Hook>(&it->second->root, 0, 0, ctx, workspace, hook, 0)) break;
        }
      }
    }

    return std::move(ctx).release_result();
  }
};

} // namespace seqtrie

#endif // seqtrie_RADIX_FOREST_H
