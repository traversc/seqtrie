#include <Rcpp.h>
#include <RcppParallel.h>

#include <set>
#include <memory>
#include <tuple>
#include "seqtrie_types.h"

#include "simple_progress/simple_progress.h"

#include <boost/numeric/ublas/matrix.hpp>
using IMatrix = boost::numeric::ublas::matrix<int>;

using namespace Rcpp;
using namespace RcppParallel;


IMatrix get_dprog_matrix(cspan query, cspan target) {
  IMatrix mat(query.size()+1, target.size()+1);
  // boundaries are the same for both levenshtein and anchored
  for(size_t i=0; i<mat.size1(); ++i) mat(i,0) = i;
  for(size_t j=1; j<mat.size2(); ++j) mat(0,j) = j;
  
  // fill it in
  for(size_t i=1; i<mat.size1(); ++i) {
    for(size_t j=1; j<mat.size2(); ++j) {
      int match_cost  = mat(i-1, j-1) + (query[i-1] == target[j-1] ? 0 : 1);
      int insert_cost = mat(i, j-1) + 1;
      int delete_cost = mat(i-1, j) + 1;
      mat(i,j) = std::min({match_cost, insert_cost, delete_cost});
    }
  }
  return mat;
}

int hamming_distance(cspan query, cspan target) {
  if(query.size() != target.size()) return NA_INTEGER;
  int distance = 0;
  for(size_t i=0; i<query.size(); ++i) {
    if(query[i] != target[i]) distance++;
  }
  return distance;
}

int levenshtein_distance(cspan query, cspan target) {
  IMatrix mat = get_dprog_matrix(query, target);
  return mat(mat.size1()-1, mat.size2()-1);
}

std::tuple<int, int, int> anchored_distance(cspan query, cspan target) {
  IMatrix mat = get_dprog_matrix(query, target);
  int distance = std::numeric_limits<int>::max();
  int query_size;
  int target_size;
  for(size_t i=0; i<mat.size1(); ++i) {
    int new_dist = mat(i, mat.size2()-1);
    if(new_dist < distance) {
      distance = new_dist;
      query_size = i;
      target_size = mat.size2()-1;
    }
  }
  for(size_t j=0; j<mat.size2(); ++j) {
    int new_dist = mat(mat.size1()-1, j);
    if(new_dist < distance) {
      distance = new_dist;
      query_size = mat.size1()-1;
      target_size = j;
    }
  }
  return std::tuple<int, int, int>(distance, query_size, target_size);
}

// if we implement more distance metrics, might want to use if constexpr (require C++17)
// Rather than defining a separate class for each
// the output matrix is parallelized over target because R matrices are column major

struct HammingDistanceWorker : public Worker {
  const std::vector<cspan> & query;
  const std::vector<cspan> & target;
  int * output; // IntegerMatrix
  trqwe::simple_progress & progress_bar;
  HammingDistanceWorker(const std::vector<cspan> & query,
                        const std::vector<cspan> & target,
                        int * output,
                        trqwe::simple_progress & progress_bar) : 
    query(query), target(target), output(output), progress_bar(progress_bar) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t j=begin; j<end; ++j) {
      for(size_t i=0; i<query.size(); ++i) {
        output[i + j*query.size()] = hamming_distance(query[i], target[j]);
      }
    }
  }
};

struct LevenshteinDistanceWorker : public Worker {
  const std::vector<cspan> & query;
  const std::vector<cspan> & target;
  int * output; // IntegerMatrix
  trqwe::simple_progress & progress_bar;
  LevenshteinDistanceWorker(const std::vector<cspan> & query,
                            const std::vector<cspan> & target,
                            int * output,
                            trqwe::simple_progress & progress_bar) : 
    query(query), target(target), output(output), progress_bar(progress_bar) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t j=begin; j<end; ++j) {
      for(size_t i=0; i<query.size(); ++i) {
        output[i + j*query.size()] = levenshtein_distance(query[i], target[j]);
      }
    }
  }
};

struct AnchoredDistanceWorker : public Worker {
  const std::vector<cspan> & query;
  const std::vector<cspan> & target;
  int * output; // IntegerMatrix
  int * query_size; // IntegerMatrix
  int * target_size; // IntegerMatrix
  trqwe::simple_progress & progress_bar;
  AnchoredDistanceWorker(const std::vector<cspan> & query,
                         const std::vector<cspan> & target,
                         int * output,
                         int * query_size,
                         int * target_size,
                         trqwe::simple_progress & progress_bar) : 
    query(query), target(target), output(output), query_size(query_size), target_size(target_size),
    progress_bar(progress_bar) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t j=begin; j<end; ++j) {
      for(size_t i=0; i<query.size(); ++i) {
        auto res = anchored_distance(query[i], target[j]);
        output[i + j*query.size()] = std::get<0>(res);
        query_size[i + j*query.size()] = std::get<1>(res);
        target_size[i + j*query.size()] = std::get<2>(res);        
      }
    }
  }
};

// [[Rcpp::export(rng = false)]]
IntegerMatrix distance_matrix(CharacterVector query, CharacterVector target, 
                              const std::string mode = "levenshtein", 
                              const int nthreads = 1, const bool show_progress = false) {
  if((mode != "levenshtein") && (mode != "hamming") && (mode != "anchored")) {
    throw std::runtime_error("Metric must be one of levenshtein, hamming or anchored");
  }
  size_t query_len = Rf_xlength(query);
  SEXP * query_ptr = STRING_PTR(query);
  size_t target_len = Rf_xlength(target);
  SEXP * target_ptr = STRING_PTR(target);
  std::vector<cspan> query_span(query_len);
  for(size_t i=0; i<query_len; ++i) query_span[i] = charsxp_to_cspan(query_ptr[i]);
  std::vector<cspan> target_span(target_len);
  for(size_t i=0; i<target_len; ++i) target_span[i] = charsxp_to_cspan(target_ptr[i]);
  
  IntegerMatrix output(query_len, target_len);
  int * output_ptr = INTEGER(output);

  trqwe::simple_progress progress_bar(target_len, show_progress);
  if(mode == "levenshtein") {
    LevenshteinDistanceWorker w(query_span, target_span, output_ptr, progress_bar);
    parallelFor(0, target_len, w, 1, nthreads);
  } else if(mode == "hamming") {
    HammingDistanceWorker w(query_span, target_span, output_ptr, progress_bar);
    parallelFor(0, target_len, w, 1, nthreads);
  } else { // anchored
    IntegerMatrix query_size(query_len, target_len);
    int * query_size_ptr = INTEGER(query_size);
    IntegerMatrix target_size(query_len, target_len);
    int * target_size_ptr = INTEGER(target_size);
    AnchoredDistanceWorker w(query_span, target_span, output_ptr, query_size_ptr, target_size_ptr, progress_bar);
    parallelFor(0, target_len, w, 1, nthreads);
    output.attr("query_size") = query_size;
    output.attr("target_size") = target_size;
  }
  return output;
}
