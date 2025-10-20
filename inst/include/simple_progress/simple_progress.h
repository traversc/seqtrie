#ifndef SIMPLE_PROGRESS_H
#define SIMPLE_PROGRESS_H

#include <atomic>
#include <thread>
#include <Rcpp.h>
#include <R_ext/Print.h>

namespace trqwe {

class simple_progress {
private:
  const size_t max;
  std::atomic<size_t> counter;
  size_t current_ticks;
  std::thread::id main_thread;
  const bool show_progress;
  static constexpr double max_ticks = 51;
public:
  simple_progress(const size_t max, const bool show_progress) : max(max), counter(0), current_ticks(0), main_thread(std::this_thread::get_id()), show_progress(show_progress) {
    if(show_progress) {
      REprintf("|----|----|----|----|----|----|----|----|----|----|\n");
      R_FlushConsole();
    }
  }
  
  size_t increment(size_t n = 1) {
    n = counter.fetch_add(n);
    if(show_progress && (std::this_thread::get_id() == main_thread)) {
      print(n);
    }
    return n;
  }
  
  void print(size_t n) {
    size_t new_ticks = static_cast<size_t>( static_cast<double>(n) / static_cast<double>(max) * max_ticks );
    size_t tick_diff = new_ticks - current_ticks;
    if(tick_diff > 0) {
      current_ticks = new_ticks;
      for(size_t i=0; i<tick_diff; ++i) REprintf("*");
    }
    R_FlushConsole();
  }
  
  // remove copy constructor/assignment and move constructor/assignment(implicit)
  simple_progress & operator=(const simple_progress&) = delete;
  simple_progress(const simple_progress&) = delete;
  
  // destructor needs to print the final output
  ~simple_progress() {
    if(show_progress) {
      print(counter.load());
      REprintf("\n");
    }
  }
};

} // end namespace

#endif // include guard
