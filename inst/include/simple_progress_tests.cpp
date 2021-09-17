#include <Rcpp.h>
#include <RcppParallel.h>

#include "simple_progress/simple_progress.h"

using namespace Rcpp;
using namespace RcppParallel;

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]

struct TestWorker : public Worker {
  std::atomic<uint64_t> output;
  simple_progress & prog;
  TestWorker(uint64_t output, simple_progress & prog) : output(output), prog(prog) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t i=begin; i<end; ++i) {
      output.fetch_add(i);
      prog.increment();
    }
  }
};

// [[Rcpp::export]]
double test_progress(int nthreads) {
  size_t n = 100000000;
  simple_progress prog(n);
  TestWorker w(0, prog);
  parallelFor(0, n, w, 1, nthreads);
  return static_cast<double>(w.output.load());
}

/*** R
test_progress(4)
*/

// #include <Rcpp.h>
// #include <atomic>
// #include <RcppParallel.h>
// 
// // [[Rcpp::depends(RcppProgress)]]
// #include <progress.hpp>
// #include <progress_bar.hpp>
// 
// using namespace Rcpp;
// using namespace RcppParallel;
// 
// // [[Rcpp::depends(RcppParallel)]]
// // [[Rcpp::plugins(cpp11)]]
// 
// struct TestWorker : public Worker {
//   std::atomic<uint64_t> output;
//   Progress & prog;
//   TestWorker(uint64_t output, Progress & prog) : output(output), prog(prog) {}
//   void operator()(std::size_t begin, std::size_t end) {
//     for(size_t i=begin; i<end; ++i) {
//       output.fetch_add(i);
//       prog.increment();
//     }
//   }
// };
// 
// // [[Rcpp::export]]
// double test_progress(int nthreads) {
//   size_t n = 100000000;
//   Progress prog(n, true);
//   TestWorker w(0, prog);
//   parallelFor(0, n, w, 1, nthreads);
//   return static_cast<double>(w.output.load());
// }
