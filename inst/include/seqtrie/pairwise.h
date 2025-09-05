// Header-only, basic O(N^2) pairwise alignment algorithms using dynamic programming
// Uses boost matrix, but maybe use something else in the future
// 
// These functions should be compatibile with C++11
// int hamming_distance_unary(query, target) // Mismatch = 1, NA if not same length
// int global_distance(query, target) // Global alignment, mismatch and gap = 1 (AKA levenshtein distance)
// int levenshtein_distance(query, target) // Alias for above
// int global_distance_linear(query, target, cost_map) // Global alignment, gap penalty is linear
// int global_distance_affine(query, target, cost_map) // Global alignment, gap penalty is affine
// tuple anchored_distance(query, target) // Anchored alignment, mismatch and gap = 1, returns a tuple of distance, query_size, target_size
// tuple anchored_distance_linear(query, target, cost_map) // Anchored alignment, gap penalty is linear
// tuple anchored_distance_affine(query, target, cost_map) // Anchored alignment, gap penalty is affine

#ifndef seqtrie_PAIRWISE_H
#define seqtrie_PAIRWISE_H

#include <set>
#include <memory>
#include <tuple>
#include <unordered_map>
#include <nonstd/span.hpp>
#include "seqtrie/utility.h"
#include <limits>
#include <boost/numeric/ublas/matrix.hpp>
// #include <boost/functional/hash.hpp>

namespace pairwise {
using IMatrix = boost::numeric::ublas::matrix<int>;
using cspan = nonstd::span<const char>;

constexpr int R_NA_INTEGER = std::numeric_limits<int>::min();

constexpr int NO_ALIGN = std::numeric_limits<int>::max() / 2; // used to represent impossible affine positions; use half INT_MAX so we dont overflow

// GAP_COST = GAP_OPEN_COST + GAP_EXTN_COST * (LENGTH - 1)

// void print_matrix(IMatrix & mat) {
//   for(size_t i=0; i<mat.size1(); ++i) {
//     for(size_t j=0; j<mat.size2(); ++j) {
//       if(mat(i,j) > 100000) {
//         std::cout << "Z ";
//       } else {
//         std::cout << mat(i,j) << " ";
//       }
//     }
//     std::cout << std::endl;
//   }
// }

// void print_pairchar_map(pairchar_map_type & cost_map) {
//   for(auto & kv : cost_map) {
//     std::cout << static_cast<int>(kv.first.first) << " " << static_cast<int>(kv.first.second) << " " << kv.second << std::endl;
//   }
// }

int hamming_distance(cspan query, cspan target) {
  if(query.size() != target.size()) return R_NA_INTEGER;
  int distance = 0;
  for(size_t i=0; i<query.size(); ++i) {
    if(query[i] != target[i]) distance++;
  }
  return distance;
}

// Starting boundaries are the same for both levenshtein and anchored
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

inline IMatrix get_dprog_matrix_linear(cspan query, cspan target, const seqtrie::CostMap & cost_map) {
  IMatrix mat(query.size()+1, target.size()+1);
  mat(0,0) = 0;
  for(size_t j=1; j<mat.size2(); ++j) mat(0,j) = mat(0,j-1) + cost_map.gap_cost; // gap in query
  for(size_t i=1; i<mat.size1(); ++i) mat(i,0) = mat(i-1,0) + cost_map.gap_cost; // gap in target
  
  // fill it in
  for(size_t i=1; i<mat.size1(); ++i) {
    for(size_t j=1; j<mat.size2(); ++j) {
      int match_cost   = mat(i-1, j-1) + cost_map.char_cost_map.at(std::make_pair(query[i-1], target[j-1]));
      int gap_in_query = mat(i, j-1) + cost_map.gap_cost;
      int gap_in_target= mat(i-1, j) + cost_map.gap_cost;
      mat(i,j) = std::min({match_cost, gap_in_query, gap_in_target});
    }
  }
  return mat;
}

inline std::tuple<IMatrix, IMatrix, IMatrix> get_dprog_matrix_affine(cspan query, cspan target, const seqtrie::CostMap & cost_map) {
  size_t size1 = query.size()+1;
  size_t size2 = target.size()+1;
  std::tuple<IMatrix, IMatrix, IMatrix> mats = std::make_tuple(IMatrix(size1, size2), IMatrix(size1, size2), IMatrix(size1, size2));
  IMatrix & M = std::get<0>(mats); // match
  IMatrix & X = std::get<1>(mats); // gap in query
  IMatrix & Y = std::get<2>(mats); // gap in target
  M(0,0) = 0;
  X(0,0) = NO_ALIGN;
  Y(0,0) = NO_ALIGN;
  for(size_t j=1; j<size2; ++j) {
    M(0,j) = NO_ALIGN;
    if(j == 1) X(0,j) = cost_map.gap_open_cost;
    else       X(0,j) = X(0,j-1) + cost_map.gap_cost;
    Y(0,j) = NO_ALIGN;
  }
  for(size_t i=1; i<size1; ++i) {
    M(i,0) = NO_ALIGN;
    X(i,0) = NO_ALIGN;
    if(i == 1) Y(i,0) = cost_map.gap_open_cost;
    else       Y(i,0) = Y(i-1,0) + cost_map.gap_cost;
  }

  // fill it in
  for(size_t i=1; i<size1; ++i) {
    for(size_t j=1; j<size2; ++j) {
      M(i,j) = cost_map.char_cost_map.at(std::make_pair(query[i-1], target[j-1])) + std::min({
        M(i-1, j-1), 
        X(i-1, j-1), 
        Y(i-1, j-1)});
      X(i,j) = std::min({
        cost_map.gap_open_cost + M(i, j-1),
        cost_map.gap_cost      + X(i, j-1),
        cost_map.gap_open_cost + Y(i, j-1)
      }); 
      Y(i,j) = std::min({
        cost_map.gap_open_cost + M(i-1, j),
        cost_map.gap_open_cost + X(i-1, j),
        cost_map.gap_cost      + Y(i-1, j)
      });
    }
  }
  return mats;
}

int global_distance(cspan query, cspan target) {
  IMatrix mat = get_dprog_matrix(query, target);
  return mat(mat.size1()-1, mat.size2()-1);
}
const auto& levenshtein_distance = global_distance; // alias

inline int global_distance_linear(cspan query, cspan target, const seqtrie::CostMap & cost_map) {
  IMatrix mat = get_dprog_matrix_linear(query, target, cost_map);
  return mat(mat.size1()-1, mat.size2()-1);
}

inline int global_distance_affine(cspan query, cspan target, const seqtrie::CostMap & cost_map) {
  // print_pairchar_map(cost_map);
  auto mats = get_dprog_matrix_affine(query, target, cost_map);
  IMatrix & M = std::get<0>(mats);
  IMatrix & X = std::get<1>(mats);
  IMatrix & Y = std::get<2>(mats);
  // std::cout << "M" << std::endl;
  // print_matrix(M);
  // std::cout << "X" << std::endl;
  // print_matrix(X);
  // std::cout << "Y" << std::endl;
  // print_matrix(Y);
  size_t xlast = M.size1()-1;
  size_t ylast = M.size2()-1;
  return std::min({M(xlast, ylast), X(xlast, ylast), Y(xlast, ylast)});
}


std::tuple<int, int, int> anchored_distance(cspan query, cspan target) {
  IMatrix mat = get_dprog_matrix(query, target);
  int distance = NO_ALIGN;
  int query_size = 0;
  int target_size = 0;
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

inline std::tuple<int, int, int> anchored_distance_linear(cspan query, cspan target, const seqtrie::CostMap & cost_map) {
  IMatrix mat = get_dprog_matrix_linear(query, target, cost_map);
  int distance = NO_ALIGN;
  int query_size = 0;
  int target_size = 0;
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

inline std::tuple<int, int, int> anchored_distance_affine(cspan query, cspan target, const seqtrie::CostMap & cost_map) {
  auto mats = get_dprog_matrix_affine(query, target, cost_map);
  IMatrix & M = std::get<0>(mats);
  IMatrix & X = std::get<1>(mats);
  IMatrix & Y = std::get<2>(mats);
  int distance = NO_ALIGN;
  int query_size = 0;
  int target_size = 0;
  size_t xlast = M.size1()-1;
  size_t ylast = M.size2()-1;
  for(size_t i=0; i<=xlast; ++i) {
    int new_dist = std::min({M(i, ylast), X(i, ylast), Y(i, ylast)});
    if(new_dist < distance) {
      distance = new_dist;
      query_size = i;
      target_size = ylast;
    }
  }
  for(size_t j=0; j<=ylast; ++j) {
    int new_dist = std::min({M(xlast, j), X(xlast, j), Y(xlast, j)});
    if(new_dist < distance) {
      distance = new_dist;
      query_size = xlast;
      target_size = j;
    }
  }
  return std::tuple<int, int, int>(distance, query_size, target_size);
}

}
#endif
