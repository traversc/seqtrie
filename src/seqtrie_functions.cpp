#include <Rcpp.h>
#include <set>
#include <memory>
#include <tuple>
#include "seqtrie_types_impl.h"

using namespace Rcpp;

// [[Rcpp::export(rng = false)]]
SEXP RadixTree_create() { return Rcpp::XPtr<RadixTree>(new RadixTree, true); }

// [[Rcpp::export(rng = false)]]
double RadixTree_size(Rcpp::XPtr<RadixTree> xp) { return xp->size(); }

// [[Rcpp::export(rng = false)]]
std::string RadixTree_print(Rcpp::XPtr<RadixTree> xp) { return xp->print(); }

// [[Rcpp::export(rng = false)]]
SEXP RadixTree_graph(Rcpp::XPtr<RadixTree> xp, const double max_depth) { return xp->graph(max_depth); }

// [[Rcpp::export(rng = false)]]
NumericVector RadixTree_insert(Rcpp::XPtr<RadixTree> xp, CharacterVector sequences) { return xp->insert(sequences); }

// [[Rcpp::export(rng = false)]]
NumericVector RadixTree_erase(Rcpp::XPtr<RadixTree> xp, CharacterVector sequences) { return xp->erase(sequences); }

// [[Rcpp::export(rng = false)]]
NumericVector RadixTree_find(Rcpp::XPtr<RadixTree> xp, CharacterVector sequences) { return xp->find(sequences); }

// [[Rcpp::export(rng = false)]]
SEXP RadixTree_to_dataframe(Rcpp::XPtr<RadixTree> xp) { return xp->to_dataframe(); }

// [[Rcpp::export(rng = false)]]
SEXP RadixTree_levenshtein_search(Rcpp::XPtr<RadixTree> xp, CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress) {
  return xp->levenshtein_search(sequences, max_distance, nthreads, show_progress);
}

// [[Rcpp::export(rng = false)]]
SEXP RadixTree_hamming_search(Rcpp::XPtr<RadixTree> xp, CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress) {
  return xp->hamming_search(sequences, max_distance, nthreads, show_progress);
}

// [[Rcpp::export(rng = false)]]
SEXP RadixTree_prefix_search(Rcpp::XPtr<RadixTree> xp, CharacterVector sequences) { return xp->prefix_search(sequences); }

// [[Rcpp::export(rng = false)]]
bool RadixTree_validate(Rcpp::XPtr<RadixTree> xp) { return xp->validate(); }