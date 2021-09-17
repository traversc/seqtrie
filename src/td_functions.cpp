#include <Rcpp.h>
#include <set>
#include <memory>
#include <tuple>
#include "treedist_types_impl.h"

using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////
// create

// [[Rcpp::export(rng = false)]]
SEXP DNATree_create() { return Rcpp::XPtr<DNATree>(new DNATree, true); }

// [[Rcpp::export(rng = false)]]
SEXP RadixTree_create() { return Rcpp::XPtr<RadixTree>(new RadixTree, true); }

// [[Rcpp::export(rng = false)]]
SEXP PrefixTree_create() {return Rcpp::XPtr<PrefixTree>(new PrefixTree, true); }

////////////////////////////////////////////////////////////////////////////////
// size
// [[Rcpp::export(rng = false)]]
double DNATree_size(Rcpp::XPtr<DNATree> xp) { return static_cast<double>(xp->size()); }

// [[Rcpp::export(rng = false)]]
double RadixTree_size(Rcpp::XPtr<RadixTree> xp) { return static_cast<double>(xp->size()); }

// [[Rcpp::export(rng = false)]]
double PrefixTree_size(Rcpp::XPtr<PrefixTree> xp) { return static_cast<double>(xp->size()); }

////////////////////////////////////////////////////////////////////////////////
// print

// [[Rcpp::export(rng = false)]]
std::string DNATree_print(Rcpp::XPtr<DNATree> xp) { return xp->print(); }

// [[Rcpp::export(rng = false)]]
std::string RadixTree_print(Rcpp::XPtr<RadixTree> xp) { return xp->print(); }

// [[Rcpp::export(rng = false)]]
std::string PrefixTree_print(Rcpp::XPtr<PrefixTree> xp) { return xp->print(); }

////////////////////////////////////////////////////////////////////////////////
// insert

// [[Rcpp::export(rng = false)]]
NumericVector DNATree_insert(Rcpp::XPtr<DNATree> xp, CharacterVector sequences) { return xp->insert(sequences); }

// [[Rcpp::export(rng = false)]]
NumericVector RadixTree_insert(Rcpp::XPtr<RadixTree> xp, CharacterVector sequences) { return xp->insert(sequences); }

// [[Rcpp::export(rng = false)]]
NumericVector PrefixTree_insert(Rcpp::XPtr<PrefixTree> xp, CharacterVector sequences) { return xp->insert(sequences); }

////////////////////////////////////////////////////////////////////////////////
// erase

// [[Rcpp::export(rng = false)]]
NumericVector DNATree_erase(Rcpp::XPtr<DNATree> xp, CharacterVector sequences) { return xp->erase(sequences); }

// [[Rcpp::export(rng = false)]]
NumericVector RadixTree_erase(Rcpp::XPtr<RadixTree> xp, CharacterVector sequences) { return xp->erase(sequences); }

// [[Rcpp::export(rng = false)]]
NumericVector PrefixTree_erase(Rcpp::XPtr<PrefixTree> xp, CharacterVector sequences) { return xp->erase(sequences); }

////////////////////////////////////////////////////////////////////////////////
// find

// [[Rcpp::export(rng = false)]]
NumericVector DNATree_find(Rcpp::XPtr<DNATree> xp, CharacterVector sequences) { return xp->find(sequences); }

// [[Rcpp::export(rng = false)]]
NumericVector RadixTree_find(Rcpp::XPtr<RadixTree> xp, CharacterVector sequences) { return xp->find(sequences); }

// [[Rcpp::export(rng = false)]]
NumericVector PrefixTree_find(Rcpp::XPtr<PrefixTree> xp, CharacterVector sequences) { return xp->find(sequences); }

////////////////////////////////////////////////////////////////////////////////
// find prefix

// [[Rcpp::export(rng = false)]]
SEXP DNATree_find_prefix(Rcpp::XPtr<DNATree> xp, CharacterVector sequences) { return xp->find_prefix(sequences); }

// [[Rcpp::export(rng = false)]]
SEXP RadixTree_find_prefix(Rcpp::XPtr<RadixTree> xp, CharacterVector sequences) { return xp->find_prefix(sequences); }

// [[Rcpp::export(rng = false)]]
SEXP PrefixTree_find_prefix(Rcpp::XPtr<PrefixTree> xp, CharacterVector sequences) { return xp->find_prefix(sequences); }

////////////////////////////////////////////////////////////////////////////////
// to_dataframe

// [[Rcpp::export(rng = false)]]
SEXP DNATree_to_dataframe(Rcpp::XPtr<DNATree> xp) { return xp->to_dataframe(); }

// [[Rcpp::export(rng = false)]]
SEXP RadixTree_to_dataframe(Rcpp::XPtr<RadixTree> xp) { return xp->to_dataframe(); }

// [[Rcpp::export(rng = false)]]
SEXP PrefixTree_to_dataframe(Rcpp::XPtr<PrefixTree> xp) { return xp->to_dataframe(); }

////////////////////////////////////////////////////////////////////////////////
// levenshtein


// [[Rcpp::export(rng = false)]]
SEXP DNATree_levenshtein(Rcpp::XPtr<DNATree> xp, CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress) {
  return xp->levenshtein_search(sequences, max_distance, nthreads, show_progress);
}

// [[Rcpp::export(rng = false)]]
SEXP RadixTree_levenshtein(Rcpp::XPtr<RadixTree> xp, CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress) {
  return xp->levenshtein_search(sequences, max_distance, nthreads, show_progress);
}

// [[Rcpp::export(rng = false)]]
SEXP PrefixTree_levenshtein(Rcpp::XPtr<PrefixTree> xp, CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress) {
  return xp->levenshtein_search(sequences, max_distance, nthreads, show_progress);
}

////////////////////////////////////////////////////////////////////////////////
// hamming

// [[Rcpp::export(rng = false)]]
SEXP DNATree_hamming(Rcpp::XPtr<DNATree> xp, CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress) {
  return xp->hamming_search(sequences, max_distance, nthreads, show_progress);
}

// [[Rcpp::export(rng = false)]]
SEXP RadixTree_hamming(Rcpp::XPtr<RadixTree> xp, CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress) {
  return xp->hamming_search(sequences, max_distance, nthreads, show_progress);
}

// [[Rcpp::export(rng = false)]]
SEXP PrefixTree_hamming(Rcpp::XPtr<PrefixTree> xp, CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress) {
  return xp->hamming_search(sequences, max_distance, nthreads, show_progress);
}
