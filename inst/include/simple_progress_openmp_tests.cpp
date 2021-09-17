#include <Rcpp.h>
#include <atomic>
#include <omp.h>

// [[Rcpp::plugins(openmp)]]

#include "simple_progress/simple_progress.h"

using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
double test_progress(int nthreads) {
  size_t n = 100000000;
  simple_progress prog(n);
  std::atomic<uint64_t> output{0};
  
  omp_set_num_threads(nthreads);
  
  #pragma omp parallel for schedule(dynamic)  
  for(size_t i=0; i<n; ++i) {
    output.fetch_add(i);
    prog.increment();
  }
  
  return static_cast<double>(output.load());
}

/*** R
test_progress(4)
*/