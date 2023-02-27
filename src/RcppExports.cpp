// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "seqtrie_types.h"
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// RadixTree_create
SEXP RadixTree_create();
RcppExport SEXP _seqtrie_RadixTree_create() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    rcpp_result_gen = Rcpp::wrap(RadixTree_create());
    return rcpp_result_gen;
END_RCPP
}
// RadixTree_size
double RadixTree_size(Rcpp::XPtr<RadixTreeCtxForR> xp);
RcppExport SEXP _seqtrie_RadixTree_size(SEXP xpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<RadixTreeCtxForR> >::type xp(xpSEXP);
    rcpp_result_gen = Rcpp::wrap(RadixTree_size(xp));
    return rcpp_result_gen;
END_RCPP
}
// RadixTree_print
std::string RadixTree_print(Rcpp::XPtr<RadixTreeCtxForR> xp);
RcppExport SEXP _seqtrie_RadixTree_print(SEXP xpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<RadixTreeCtxForR> >::type xp(xpSEXP);
    rcpp_result_gen = Rcpp::wrap(RadixTree_print(xp));
    return rcpp_result_gen;
END_RCPP
}
// RadixTree_graph
SEXP RadixTree_graph(Rcpp::XPtr<RadixTreeCtxForR> xp, const double max_depth);
RcppExport SEXP _seqtrie_RadixTree_graph(SEXP xpSEXP, SEXP max_depthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<RadixTreeCtxForR> >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< const double >::type max_depth(max_depthSEXP);
    rcpp_result_gen = Rcpp::wrap(RadixTree_graph(xp, max_depth));
    return rcpp_result_gen;
END_RCPP
}
// RadixTree_insert
LogicalVector RadixTree_insert(Rcpp::XPtr<RadixTreeCtxForR> xp, CharacterVector sequences);
RcppExport SEXP _seqtrie_RadixTree_insert(SEXP xpSEXP, SEXP sequencesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<RadixTreeCtxForR> >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type sequences(sequencesSEXP);
    rcpp_result_gen = Rcpp::wrap(RadixTree_insert(xp, sequences));
    return rcpp_result_gen;
END_RCPP
}
// RadixTree_erase
LogicalVector RadixTree_erase(Rcpp::XPtr<RadixTreeCtxForR> xp, CharacterVector sequences);
RcppExport SEXP _seqtrie_RadixTree_erase(SEXP xpSEXP, SEXP sequencesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<RadixTreeCtxForR> >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type sequences(sequencesSEXP);
    rcpp_result_gen = Rcpp::wrap(RadixTree_erase(xp, sequences));
    return rcpp_result_gen;
END_RCPP
}
// RadixTree_find
LogicalVector RadixTree_find(Rcpp::XPtr<RadixTreeCtxForR> xp, CharacterVector sequences);
RcppExport SEXP _seqtrie_RadixTree_find(SEXP xpSEXP, SEXP sequencesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<RadixTreeCtxForR> >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type sequences(sequencesSEXP);
    rcpp_result_gen = Rcpp::wrap(RadixTree_find(xp, sequences));
    return rcpp_result_gen;
END_RCPP
}
// RadixTree_to_vector
SEXP RadixTree_to_vector(Rcpp::XPtr<RadixTreeCtxForR> xp);
RcppExport SEXP _seqtrie_RadixTree_to_vector(SEXP xpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<RadixTreeCtxForR> >::type xp(xpSEXP);
    rcpp_result_gen = Rcpp::wrap(RadixTree_to_vector(xp));
    return rcpp_result_gen;
END_RCPP
}
// RadixTree_levenshtein_search
SEXP RadixTree_levenshtein_search(Rcpp::XPtr<RadixTreeCtxForR> xp, CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress);
RcppExport SEXP _seqtrie_RadixTree_levenshtein_search(SEXP xpSEXP, SEXP sequencesSEXP, SEXP max_distanceSEXP, SEXP nthreadsSEXP, SEXP show_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<RadixTreeCtxForR> >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type sequences(sequencesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type max_distance(max_distanceSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type show_progress(show_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(RadixTree_levenshtein_search(xp, sequences, max_distance, nthreads, show_progress));
    return rcpp_result_gen;
END_RCPP
}
// RadixTree_hamming_search
SEXP RadixTree_hamming_search(Rcpp::XPtr<RadixTreeCtxForR> xp, CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress);
RcppExport SEXP _seqtrie_RadixTree_hamming_search(SEXP xpSEXP, SEXP sequencesSEXP, SEXP max_distanceSEXP, SEXP nthreadsSEXP, SEXP show_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<RadixTreeCtxForR> >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type sequences(sequencesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type max_distance(max_distanceSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type show_progress(show_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(RadixTree_hamming_search(xp, sequences, max_distance, nthreads, show_progress));
    return rcpp_result_gen;
END_RCPP
}
// RadixTree_anchored_search
SEXP RadixTree_anchored_search(Rcpp::XPtr<RadixTreeCtxForR> xp, CharacterVector sequences, IntegerVector max_distance, const int nthreads, const bool show_progress);
RcppExport SEXP _seqtrie_RadixTree_anchored_search(SEXP xpSEXP, SEXP sequencesSEXP, SEXP max_distanceSEXP, SEXP nthreadsSEXP, SEXP show_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<RadixTreeCtxForR> >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type sequences(sequencesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type max_distance(max_distanceSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type show_progress(show_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(RadixTree_anchored_search(xp, sequences, max_distance, nthreads, show_progress));
    return rcpp_result_gen;
END_RCPP
}
// RadixTree_prefix_search
SEXP RadixTree_prefix_search(Rcpp::XPtr<RadixTreeCtxForR> xp, CharacterVector sequences);
RcppExport SEXP _seqtrie_RadixTree_prefix_search(SEXP xpSEXP, SEXP sequencesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<RadixTreeCtxForR> >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type sequences(sequencesSEXP);
    rcpp_result_gen = Rcpp::wrap(RadixTree_prefix_search(xp, sequences));
    return rcpp_result_gen;
END_RCPP
}
// RadixTree_validate
bool RadixTree_validate(Rcpp::XPtr<RadixTreeCtxForR> xp);
RcppExport SEXP _seqtrie_RadixTree_validate(SEXP xpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::XPtr<RadixTreeCtxForR> >::type xp(xpSEXP);
    rcpp_result_gen = Rcpp::wrap(RadixTree_validate(xp));
    return rcpp_result_gen;
END_RCPP
}
// distance_matrix
IntegerMatrix distance_matrix(CharacterVector query, CharacterVector target, const std::string mode, const int nthreads, const bool show_progress);
RcppExport SEXP _seqtrie_distance_matrix(SEXP querySEXP, SEXP targetSEXP, SEXP modeSEXP, SEXP nthreadsSEXP, SEXP show_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type query(querySEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const std::string >::type mode(modeSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type show_progress(show_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(distance_matrix(query, target, mode, nthreads, show_progress));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_seqtrie_RadixTree_create", (DL_FUNC) &_seqtrie_RadixTree_create, 0},
    {"_seqtrie_RadixTree_size", (DL_FUNC) &_seqtrie_RadixTree_size, 1},
    {"_seqtrie_RadixTree_print", (DL_FUNC) &_seqtrie_RadixTree_print, 1},
    {"_seqtrie_RadixTree_graph", (DL_FUNC) &_seqtrie_RadixTree_graph, 2},
    {"_seqtrie_RadixTree_insert", (DL_FUNC) &_seqtrie_RadixTree_insert, 2},
    {"_seqtrie_RadixTree_erase", (DL_FUNC) &_seqtrie_RadixTree_erase, 2},
    {"_seqtrie_RadixTree_find", (DL_FUNC) &_seqtrie_RadixTree_find, 2},
    {"_seqtrie_RadixTree_to_vector", (DL_FUNC) &_seqtrie_RadixTree_to_vector, 1},
    {"_seqtrie_RadixTree_levenshtein_search", (DL_FUNC) &_seqtrie_RadixTree_levenshtein_search, 5},
    {"_seqtrie_RadixTree_hamming_search", (DL_FUNC) &_seqtrie_RadixTree_hamming_search, 5},
    {"_seqtrie_RadixTree_anchored_search", (DL_FUNC) &_seqtrie_RadixTree_anchored_search, 5},
    {"_seqtrie_RadixTree_prefix_search", (DL_FUNC) &_seqtrie_RadixTree_prefix_search, 2},
    {"_seqtrie_RadixTree_validate", (DL_FUNC) &_seqtrie_RadixTree_validate, 1},
    {"_seqtrie_distance_matrix", (DL_FUNC) &_seqtrie_distance_matrix, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_seqtrie(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
