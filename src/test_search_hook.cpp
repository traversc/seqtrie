// source file for testing early stop hook functionality
// This is not useful at the R level, so only exported at the C++ library layer only

#include "seqtrie_types.h"

using namespace Rcpp;

namespace {

struct StopAfterFirst {
  bool operator()(const SeqTrie::search_context &) const noexcept { return false; } // stop immediately after first match is added
};

SeqTrie::RadixTreeR build_simple_tree() {
  SeqTrie::RadixTreeR root;
  const std::vector<std::string> seqs = {"AAAA", "AAAT", "AAAAT"};
  for(size_t i = 0; i < seqs.size(); ++i) {
    const std::string & s = seqs[i];
    root.insert(nonstd::span<const char>(s.data(), s.size()), static_cast<SeqTrie::RadixTreeR::index_type>(SeqTrie::posidx + i));
  }
  return root;
}

CostMap make_unit_cost_map() {
  CostMap cm;
  const std::vector<char> alphabet = {'A', 'T'};
  for(char a : alphabet) {
    for(char b : alphabet) {
      cm.char_cost_map.emplace(std::make_pair(a, b), (a == b) ? 0 : 1);
    }
  }
  cm.gap_cost = 1;
  cm.gap_open_cost = 2;
  return cm;
}

DataFrame ctx_to_df(const SeqTrie::search_context & ctx) {
  const std::string query_str(ctx.query.begin(), ctx.query.end());
  CharacterVector query(ctx.match.size(), query_str);
  CharacterVector target(ctx.match.size());
  IntegerVector distance(ctx.distance.size());
  for(size_t i = 0; i < ctx.match.size(); ++i) {
    auto seq = ctx.match[i]->template sequence<SeqTrie::array_r<char>>();
    SET_STRING_ELT(target, i, to_charsxp(seq));
    distance[i] = ctx.distance[i];
  }
  return DataFrame::create(_["query"] = query,
                           _["target"] = target,
                           _["distance"] = distance,
                           _["stringsAsFactors"] = false);
}

List run_with_hook(const SeqTrie::search_context & full_ctx, const SeqTrie::search_context & early_ctx) {
  return List::create(_["full"] = ctx_to_df(full_ctx),
                      _["early"] = ctx_to_df(early_ctx));
}

} // namespace

// [[Rcpp::export(rng = false)]]
List test_search_hook() {
  SeqTrie::RadixTreeR root = build_simple_tree();
  const std::string query_str = "AAAA";
  auto query_span = nonstd::span<const char>(query_str.data(), query_str.size());
  CostMap cm = make_unit_cost_map();
  StopAfterFirst hook;

  List out;
  out["hamming"] = run_with_hook(
    root.hamming_search(query_span, 1),
    root.hamming_search(query_span, 1, hook)
  );
  out["global_unit"] = run_with_hook(
    root.global_search(query_span, 1),
    root.global_search(query_span, 1, hook)
  );
  out["anchored_unit"] = run_with_hook(
    root.anchored_search(query_span, 1),
    root.anchored_search(query_span, 1, hook)
  );
  out["global_linear"] = run_with_hook(
    root.global_search_linear(query_span, 1, cm),
    root.global_search_linear(query_span, 1, cm, hook)
  );
  out["anchored_linear"] = run_with_hook(
    root.anchored_search_linear(query_span, 1, cm),
    root.anchored_search_linear(query_span, 1, cm, hook)
  );
  out["global_affine"] = run_with_hook(
    root.global_search_affine(query_span, 1, cm),
    root.global_search_affine(query_span, 1, cm, hook)
  );
  out["anchored_affine"] = run_with_hook(
    root.anchored_search_affine(query_span, 1, cm),
    root.anchored_search_affine(query_span, 1, cm, hook)
  );
  out["single_gap"] = run_with_hook(
    root.single_gap_search(query_span, 1, 1),
    root.single_gap_search(query_span, 1, 1, hook)
  );
  return out;
}
