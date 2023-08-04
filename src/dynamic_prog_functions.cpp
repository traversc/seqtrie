#include <Rcpp.h>
#include <RcppParallel.h>

#include <set>
#include <memory>
#include <tuple>
#include <unordered_map>
#include "seqtrie_types.h"

#include "simple_progress/simple_progress.h"

#include <boost/numeric/ublas/matrix.hpp>
using IMatrix = boost::numeric::ublas::matrix<int>;

// Functions
//// get_dprog_matrix(cspan query, cspan target)
//// get_dprog_matrix(cspan query, cspan target, pairchar_map & cost_map)
//// hamming_distance(cspan query, cspan target)
//// levenshtein_distance(cspan query, cspan target)
//// levenshtein_distance(cspan query, cspan target, pairchar_map & cost_map)
//// anchored_distance(cspan query, cspan target)
//// anchored_distance(cspan query, cspan target, pairchar_map & cost_map)

// Exported functions
//// IntegerVector dist_pairwise(CharacterVector query, CharacterVector target, 
                            // const std::string mode = "levenshtein", 
                            // Rcpp::Nullable<IntegerMatrix> cost_matrix = R_NilValue,
                            // const int nthreads = 1, const bool show_progress = false)
//// IntegerMatrix dist_matrix(CharacterVector query, CharacterVector target,
                            // const std::string mode = "levenshtein", 
                            // Rcpp::Nullable<IntegerMatrix> cost_matrix = R_NilValue,
                            // const int nthreads = 1, const bool show_progress = false)


// Multithreading worker classes
// Each class corresponds to a different parameter set
// This is kind of a lot. In the future, might want to consolidate with some sort of 
// template magic or just eat the (probably small) cost of branching

// Distance Matrix workers
//// HammingDistanceWorker
//// LevenshteinDistanceWorker
//// AnchoredDistanceWorker

// Distance Matrix workers with custom cost
//// LevenshteinDistanceWorkerWithCostMap
//// AnchoredDistanceWorkerWithCostMap

// Pairwise Distance workers
//// HammingDistanceWorkerPairwise
//// LevenshteinDistanceWorkerPairwise
//// AnchoredDistanceWorkerPairwise

//// LevenshteinDistanceWorkerWithCostMapPairwise
//// AnchoredDistanceWorkerWithCostMapPairwise

// included in utility.h
// hash for unordered_map with std::pair<char, char> as key
// namespace std {
//   template <>
//   struct hash<std::pair<char, char>> {
//     size_t operator()(const std::pair<char, char> & p) const {
//       return ((p.first + 128) << 8) + (p.second + 128);
//     }
//   };
// }

using namespace Rcpp;
using namespace RcppParallel;

// boundaries are the same for both levenshtein and anchored
IMatrix get_dprog_matrix(cspan query, cspan target) {
  IMatrix mat(query.size()+1, target.size()+1);
  for(size_t j=1; j<mat.size2(); ++j) mat(0,j) = j;
  for(size_t i=0; i<mat.size1(); ++i) mat(i,0) = i;
  
  // fill it in
  for(size_t i=1; i<mat.size1(); ++i) {
    for(size_t j=1; j<mat.size2(); ++j) {
      int match_cost  = mat(i-1, j-1) + (query[i-1] == target[j-1] ? 0 : 1);
      int gap_in_query = mat(i, j-1) + 1;
      int gap_in_target = mat(i-1, j) + 1;
      mat(i,j) = std::min({match_cost, gap_in_query, gap_in_target});
    }
  }
  return mat;
}

IMatrix get_dprog_matrix(cspan query, cspan target, pairchar_map & cost_map) {
  IMatrix mat(query.size()+1, target.size()+1);
  for(size_t j=1; j<mat.size2(); ++j) mat(0,j) = cost_map.at(pairchar('\0', target[j])) * j; // gap in query
  for(size_t i=0; i<mat.size1(); ++i) mat(i,0) = cost_map.at(pairchar(query[i], '\0')) * i; // gap in target
  
  // fill it in
  for(size_t i=1; i<mat.size1(); ++i) {
    for(size_t j=1; j<mat.size2(); ++j) {
      int match_cost  = mat(i-1, j-1) + cost_map.at(pairchar(query[i-1], target[j-1]));
      int gap_in_query = mat(i, j-1) + cost_map.at(pairchar('\0', target[j-1]));
      int gap_in_target = mat(i-1, j) + cost_map.at(pairchar(query[i-1], '\0'));
      mat(i,j) = std::min({match_cost, gap_in_query, gap_in_target});
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

std::tuple<int, int, int> anchored_distance(cspan query, cspan target, pairchar_map & cost_map) {
  IMatrix mat = get_dprog_matrix(query, target, cost_map);
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

int levenshtein_distance(cspan query, cspan target, pairchar_map & cost_map) {
  IMatrix mat = get_dprog_matrix(query, target, cost_map);
  return mat(mat.size1()-1, mat.size2()-1);
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

struct LevenshteinDistanceWorkerWithCostMap : public Worker {
  const std::vector<cspan> & query;
  const std::vector<cspan> & target;
  int * output; // IntegerMatrix
  pairchar_map & cost_map;
  trqwe::simple_progress & progress_bar;
  LevenshteinDistanceWorkerWithCostMap(const std::vector<cspan> & query,
                            const std::vector<cspan> & target,
                            int * output,
                            pairchar_map & cost_map,
                            trqwe::simple_progress & progress_bar) : 
    query(query), target(target), output(output), cost_map(cost_map), progress_bar(progress_bar) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t j=begin; j<end; ++j) {
      for(size_t i=0; i<query.size(); ++i) {
        output[i + j*query.size()] = levenshtein_distance(query[i], target[j], cost_map);
      }
    }
  }
};

struct AnchoredDistanceWorkerWithCostMap : public Worker {
  const std::vector<cspan> & query;
  const std::vector<cspan> & target;
  int * output; // IntegerMatrix
  int * query_size; // IntegerMatrix
  int * target_size; // IntegerMatrix
  pairchar_map & cost_map;
  trqwe::simple_progress & progress_bar;
  AnchoredDistanceWorkerWithCostMap(const std::vector<cspan> & query,
                         const std::vector<cspan> & target,
                         int * output,
                         int * query_size,
                         int * target_size,
                         pairchar_map & cost_map,
                         trqwe::simple_progress & progress_bar) : 
    query(query), target(target), output(output), query_size(query_size), target_size(target_size),
    cost_map(cost_map), progress_bar(progress_bar) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t j=begin; j<end; ++j) {
      for(size_t i=0; i<query.size(); ++i) {
        auto res = anchored_distance(query[i], target[j], cost_map);
        output[i + j*query.size()] = std::get<0>(res);
        query_size[i + j*query.size()] = std::get<1>(res);
        target_size[i + j*query.size()] = std::get<2>(res);        
      }
    }
  }
};

// [[Rcpp::export(rng = false)]]
IntegerMatrix c_dist_matrix(CharacterVector query, CharacterVector target, 
                          const std::string mode = "levenshtein", 
                          Rcpp::Nullable<IntegerMatrix> cost_matrix = R_NilValue,
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
    if(cost_matrix.isNotNull()) {
      IntegerMatrix cost_matrix_(cost_matrix);
      pairchar_map cost_map = convert_cost_matrix(cost_matrix_);
      LevenshteinDistanceWorkerWithCostMap w(query_span, target_span, output_ptr, cost_map, progress_bar);
      parallelFor(0, target_len, w, 1, nthreads);
    } else {
      LevenshteinDistanceWorker w(query_span, target_span, output_ptr, progress_bar);
      parallelFor(0, target_len, w, 1, nthreads);
    }
  } else if(mode == "hamming") {
    HammingDistanceWorker w(query_span, target_span, output_ptr, progress_bar);
    parallelFor(0, target_len, w, 1, nthreads);
  } else { // anchored
    IntegerMatrix query_size(query_len, target_len);
    int * query_size_ptr = INTEGER(query_size);
    IntegerMatrix target_size(query_len, target_len);
    int * target_size_ptr = INTEGER(target_size);
    if(cost_matrix.isNotNull()) {
      IntegerMatrix cost_matrix_(cost_matrix);
      pairchar_map cost_map = convert_cost_matrix(cost_matrix_);
      AnchoredDistanceWorkerWithCostMap w(query_span, target_span, output_ptr,  query_size_ptr, target_size_ptr, cost_map, progress_bar);
      parallelFor(0, target_len, w, 1, nthreads);
    } else {
      AnchoredDistanceWorker w(query_span, target_span, output_ptr,  query_size_ptr, target_size_ptr, progress_bar);
      parallelFor(0, target_len, w, 1, nthreads);
    }
    output.attr("query_size") = query_size;
    output.attr("target_size") = target_size;
  }
  return output;
}

struct HammingDistanceWorkerPairwise : public Worker {
  const std::vector<cspan> & query;
  const std::vector<cspan> & target;
  int * output; // IntegerVector
  trqwe::simple_progress & progress_bar;
  HammingDistanceWorkerPairwise(const std::vector<cspan> & query,
                                const std::vector<cspan> & target,
                                int * output,
                                trqwe::simple_progress & progress_bar) : 
    query(query), target(target), output(output), progress_bar(progress_bar) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t i=begin; i<end; ++i) {
      output[i] = hamming_distance(query[i], target[i]);
    }
  }
};

struct LevenshteinDistanceWorkerPairwise : public Worker {
  const std::vector<cspan> & query;
  const std::vector<cspan> & target;
  int * output; // IntegerVector
  trqwe::simple_progress & progress_bar;
  LevenshteinDistanceWorkerPairwise(const std::vector<cspan> & query,
                                    const std::vector<cspan> & target,
                                    int * output,
                                    trqwe::simple_progress & progress_bar) : 
    query(query), target(target), output(output), progress_bar(progress_bar) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t i=begin; i<end; ++i) {
      output[i] = levenshtein_distance(query[i], target[i]);
    }
  }
};

struct AnchoredDistanceWorkerPairwise : public Worker {
  const std::vector<cspan> & query;
  const std::vector<cspan> & target;
  int * output; // IntegerVector
  int * query_size; // IntegerVector
  int * target_size; // IntegerVector
  trqwe::simple_progress & progress_bar;
  AnchoredDistanceWorkerPairwise(const std::vector<cspan> & query,
                                 const std::vector<cspan> & target,
                                 int * output,
                                 int * query_size,
                                 int * target_size,
                                 trqwe::simple_progress & progress_bar) : 
    query(query), target(target), output(output), query_size(query_size), target_size(target_size),
    progress_bar(progress_bar) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t i=begin; i<end; ++i) {
      auto res = anchored_distance(query[i], target[i]);
      output[i] = std::get<0>(res);
      query_size[i] = std::get<1>(res);
      target_size[i] = std::get<2>(res);        
    }
  }
};

struct LevenshteinDistanceWorkerWithCostMapPairwise : public Worker {
  const std::vector<cspan> & query;
  const std::vector<cspan> & target;
  int * output; // IntegerVector
  pairchar_map & cost_map;
  trqwe::simple_progress & progress_bar;
  LevenshteinDistanceWorkerWithCostMapPairwise(const std::vector<cspan> & query,
                                    const std::vector<cspan> & target,
                                    int * output,
                                    pairchar_map & cost_map,
                                    trqwe::simple_progress & progress_bar) : 
    query(query), target(target), output(output), cost_map(cost_map), progress_bar(progress_bar) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t i=begin; i<end; ++i) {
      output[i] = levenshtein_distance(query[i], target[i], cost_map);
    }
  }
};

struct AnchoredDistanceWorkerWithCostMapPairwise : public Worker {
  const std::vector<cspan> & query;
  const std::vector<cspan> & target;
  int * output; // IntegerVector
  int * query_size; // IntegerVector
  int * target_size; // IntegerVector
  pairchar_map & cost_map;
  trqwe::simple_progress & progress_bar;
  AnchoredDistanceWorkerWithCostMapPairwise(const std::vector<cspan> & query,
                                 const std::vector<cspan> & target,
                                 int * output,
                                 int * query_size,
                                 int * target_size,
                                 pairchar_map & cost_map,
                                 trqwe::simple_progress & progress_bar) : 
    query(query), target(target), output(output), query_size(query_size), target_size(target_size),
    cost_map(cost_map), progress_bar(progress_bar) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t i=begin; i<end; ++i) {
      auto res = anchored_distance(query[i], target[i], cost_map);
      output[i] = std::get<0>(res);
      query_size[i] = std::get<1>(res);
      target_size[i] = std::get<2>(res);        
    }
  }
};

// [[Rcpp::export(rng = false)]]
IntegerVector c_dist_pairwise(CharacterVector query, CharacterVector target, 
                            const std::string mode = "levenshtein", 
                            Rcpp::Nullable<IntegerMatrix> cost_matrix = R_NilValue,
                            const int nthreads = 1, const bool show_progress = false) {
  if((mode != "levenshtein") && (mode != "hamming") && (mode != "anchored")) {
    throw std::runtime_error("Metric must be one of levenshtein, hamming or anchored");
  }
  size_t N = Rf_xlength(query);
  SEXP * query_ptr = STRING_PTR(query);
  size_t target_len = Rf_xlength(target);
  SEXP * target_ptr = STRING_PTR(target);
  if(N != target_len) {
    throw std::runtime_error("query and target must be the same length");
  }

  std::vector<cspan> query_span(N);
  for(size_t i=0; i<N; ++i) query_span[i] = charsxp_to_cspan(query_ptr[i]);
  std::vector<cspan> target_span(N);
  for(size_t i=0; i<N; ++i) target_span[i] = charsxp_to_cspan(target_ptr[i]);
  
  IntegerVector output(N);
  int * output_ptr = INTEGER(output);

  trqwe::simple_progress progress_bar(N, show_progress);
  if(mode == "levenshtein") {
    if(cost_matrix.isNotNull()) {
      IntegerMatrix cost_matrix_(cost_matrix);
      pairchar_map cost_map = convert_cost_matrix(cost_matrix_);
      LevenshteinDistanceWorkerWithCostMapPairwise w(query_span, target_span, output_ptr, cost_map, progress_bar);
      parallelFor(0, N, w, 1, nthreads);
    } else {
      LevenshteinDistanceWorkerPairwise w(query_span, target_span, output_ptr, progress_bar);
      parallelFor(0, N, w, 1, nthreads);
    }
  } else if(mode == "hamming") {
    HammingDistanceWorkerPairwise w(query_span, target_span, output_ptr, progress_bar);
    parallelFor(0, N, w, 1, nthreads);
  } else { // anchored
    IntegerVector query_size(N);
    int * query_size_ptr = INTEGER(query_size);
    IntegerVector target_size(N);
    int * target_size_ptr = INTEGER(target_size);
    if(cost_matrix.isNotNull()) {
      IntegerMatrix cost_matrix_(cost_matrix);
      pairchar_map cost_map = convert_cost_matrix(cost_matrix_);
      AnchoredDistanceWorkerWithCostMapPairwise w(query_span, target_span, output_ptr, query_size_ptr, target_size_ptr, cost_map, progress_bar);
      parallelFor(0, N, w, 1, nthreads);
    } else {
      AnchoredDistanceWorkerPairwise w(query_span, target_span, output_ptr, query_size_ptr, target_size_ptr, progress_bar);
      parallelFor(0, N, w, 1, nthreads);
    }
    output.attr("query_size") = query_size;
    output.attr("target_size") = target_size;
  }
  return output;
}
