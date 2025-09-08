#include "seqtrie_types.h"
#include "seqtrie/pairwise.h"
#include "simple_progress/simple_progress.h"

// All input parameters should be checked in R, so any error thrown here is an internal error
// [[Rcpp::export(rng = false)]]
IntegerMatrix c_dist_matrix(CharacterVector query, CharacterVector target, 
                          const std::string mode = "global", // global, anchored or hamming
                          Rcpp::Nullable<IntegerMatrix> cost_matrix = R_NilValue,
                          const int gap_cost = 1,
                          const int gap_open_cost = 0,
                          const int nthreads = 1, const bool show_progress = false) {
  std::vector<cspan> query_span = strsxp_to_cspan(query);
  std::vector<cspan> target_span = strsxp_to_cspan(target);

  trqwe::simple_progress progress_bar(target_span.size(), show_progress);
  IntegerMatrix output(query_span.size(), target_span.size());
  int * output_ptr = INTEGER(output);
  auto algo = decide_alignment_algo(mode, cost_matrix, gap_cost, gap_open_cost);
  if(algo == AlignmentAlgo::Hamming) {
    do_parallel_for([&query_span, &target_span, &progress_bar, output_ptr](size_t begin, size_t end) {
      for(size_t j=begin; j<end; ++j) {
        for(size_t i=0; i<query_span.size(); ++i) {
          output_ptr[i + j*query_span.size()] = pairwise::hamming_distance(query_span[i], target_span[j]);
        }
        progress_bar.increment();
      }
    }, 0, target_span.size(), 1, nthreads);
  } else if(mode == "global" || mode == "gb" || mode == "lv" || mode == "levenshtein") {
    if(algo == AlignmentAlgo::GlobalUnit) {
      // If a boolean cost_matrix with unit gaps was provided, prefer Myers; else unit DP
      bool use_myers = false;
      CostMap cost_map;
      if (!cost_matrix.isNull()) {
        cost_map = convert_cost_matrix(cost_matrix.get(), gap_cost, gap_open_cost);
        use_myers = (gap_cost == 1 && gap_open_cost == 0);
        if (use_myers) {
          for (const auto & kv : cost_map.char_cost_map) { if (kv.second != 0 && kv.second != 1) { use_myers = false; break; } }
        }
      }
      if (use_myers) {
        do_parallel_for([&query_span, &target_span, &cost_map, &progress_bar, output_ptr](size_t begin, size_t end) {
          for(size_t j=begin; j<end; ++j) {
            for(size_t i=0; i<query_span.size(); ++i) {
              output_ptr[i + j*query_span.size()] = pairwise::global_distance_myers(query_span[i], target_span[j], cost_map);
            }
            progress_bar.increment();
          }
        }, 0, target_span.size(), 1, nthreads);
      } else {
        do_parallel_for([&query_span, &target_span, &progress_bar, output_ptr](size_t begin, size_t end) {
          for(size_t j=begin; j<end; ++j) {
            for(size_t i=0; i<query_span.size(); ++i) {
              output_ptr[i + j*query_span.size()] = pairwise::global_distance(query_span[i], target_span[j]);
            }
            progress_bar.increment();
          }
        }, 0, target_span.size(), 1, nthreads);
      }
    } else if(algo == AlignmentAlgo::GlobalLinear) {
      CostMap cost_map = convert_cost_matrix(cost_matrix.get(), gap_cost, gap_open_cost);
      bool boolean_subs = true;
      for (const auto & kv : cost_map.char_cost_map) {
        if (kv.second != 0 && kv.second != 1) { boolean_subs = false; break; }
      }
      const bool unit_gaps = (gap_cost == 1) && (gap_open_cost == 0);
      if (boolean_subs && unit_gaps) {
        do_parallel_for([&query_span, &target_span, &cost_map, &progress_bar, output_ptr](size_t begin, size_t end) {
          for(size_t j=begin; j<end; ++j) {
            for(size_t i=0; i<query_span.size(); ++i) {
              output_ptr[i + j*query_span.size()] = pairwise::global_distance_myers(query_span[i], target_span[j], cost_map);
            }
            progress_bar.increment();
          }
        }, 0, target_span.size(), 1, nthreads);
      } else {
        do_parallel_for([&query_span, &target_span, &cost_map, &progress_bar, output_ptr](size_t begin, size_t end) {
          for(size_t j=begin; j<end; ++j) {
            for(size_t i=0; i<query_span.size(); ++i) {
              output_ptr[i + j*query_span.size()] = pairwise::global_distance_linear(query_span[i], target_span[j], cost_map);
            }
            progress_bar.increment();
          }
        }, 0, target_span.size(), 1, nthreads);
      }
    } else { // GlobalAffine
      CostMap cost_map = convert_cost_matrix(cost_matrix.get(), gap_cost, gap_open_cost);
      do_parallel_for([&query_span, &target_span, &cost_map, &progress_bar, output_ptr](size_t begin, size_t end) {
        for(size_t j=begin; j<end; ++j) {
          for(size_t i=0; i<query_span.size(); ++i) {
            output_ptr[i + j*query_span.size()] = pairwise::global_distance_affine(query_span[i], target_span[j], cost_map);
          }
          progress_bar.increment();
        }
      }, 0, target_span.size(), 1, nthreads);
    }
  } else { // anchored
    IntegerMatrix query_size(query_span.size(), target_span.size());
    IntegerMatrix target_size(query_span.size(), target_span.size());
    int * query_size_ptr = INTEGER(query_size);
    int * target_size_ptr = INTEGER(target_size);
    if(algo == AlignmentAlgo::AnchoredUnit) {
      bool use_myers = false;
      CostMap cost_map;
      if (!cost_matrix.isNull()) {
        cost_map = convert_cost_matrix(cost_matrix.get(), gap_cost, gap_open_cost);
        use_myers = (gap_cost == 1 && gap_open_cost == 0);
        if (use_myers) {
          for (const auto & kv : cost_map.char_cost_map) { if (kv.second != 0 && kv.second != 1) { use_myers = false; break; } }
        }
      }
      if (use_myers) {
        do_parallel_for([&query_span, &target_span, &cost_map, &progress_bar, output_ptr, query_size_ptr, target_size_ptr](size_t begin, size_t end) {
          for(size_t j=begin; j<end; ++j) {
            for(size_t i=0; i<query_span.size(); ++i) {
              auto res = pairwise::anchored_distance_myers(query_span[i], target_span[j], cost_map);
              output_ptr[i + j*query_span.size()] = std::get<0>(res);
              query_size_ptr[i + j*query_span.size()] = std::get<1>(res);
              target_size_ptr[i + j*query_span.size()] = std::get<2>(res);
            }
            progress_bar.increment();
          }
        }, 0, target_span.size(), 1, nthreads);
      } else {
        do_parallel_for([&query_span, &target_span, &progress_bar, output_ptr, query_size_ptr, target_size_ptr](size_t begin, size_t end) {
          for(size_t j=begin; j<end; ++j) {
            for(size_t i=0; i<query_span.size(); ++i) {
              auto res = pairwise::anchored_distance(query_span[i], target_span[j]);
              output_ptr[i + j*query_span.size()] = std::get<0>(res);
              query_size_ptr[i + j*query_span.size()] = std::get<1>(res);
              target_size_ptr[i + j*query_span.size()] = std::get<2>(res);
            }
            progress_bar.increment();
          }
        }, 0, target_span.size(), 1, nthreads);
      }
    } else if(algo == AlignmentAlgo::AnchoredLinear) {
      CostMap cost_map = convert_cost_matrix(cost_matrix.get(), gap_cost, gap_open_cost);
      bool boolean_subs = true;
      for (const auto & kv : cost_map.char_cost_map) { if (kv.second != 0 && kv.second != 1) { boolean_subs = false; break; } }
      const bool unit_gaps = (gap_cost == 1) && (gap_open_cost == 0);
      if (boolean_subs && unit_gaps) {
        do_parallel_for([&query_span, &target_span, &cost_map, &progress_bar, output_ptr, query_size_ptr, target_size_ptr](size_t begin, size_t end) {
          for(size_t j=begin; j<end; ++j) {
            for(size_t i=0; i<query_span.size(); ++i) {
              auto res = pairwise::anchored_distance_myers(query_span[i], target_span[j], cost_map);
              output_ptr[i + j*query_span.size()] = std::get<0>(res);
              query_size_ptr[i + j*query_span.size()] = std::get<1>(res);
              target_size_ptr[i + j*query_span.size()] = std::get<2>(res);
            }
            progress_bar.increment();
          }
        }, 0, target_span.size(), 1, nthreads);
      } else {
        do_parallel_for([&query_span, &target_span, &cost_map, &progress_bar, output_ptr, query_size_ptr, target_size_ptr](size_t begin, size_t end) {
          for(size_t j=begin; j<end; ++j) {
            for(size_t i=0; i<query_span.size(); ++i) {
              auto res = pairwise::anchored_distance_linear(query_span[i], target_span[j], cost_map);
              output_ptr[i + j*query_span.size()] = std::get<0>(res);
              query_size_ptr[i + j*query_span.size()] = std::get<1>(res);
              target_size_ptr[i + j*query_span.size()] = std::get<2>(res);
            }
            progress_bar.increment();
          }
        }, 0, target_span.size(), 1, nthreads);
      }
    } else { // AnchoredAffine
      CostMap cost_map = convert_cost_matrix(cost_matrix.get(), gap_cost, gap_open_cost);
      do_parallel_for([&query_span, &target_span, &cost_map, &progress_bar, output_ptr, query_size_ptr, target_size_ptr](size_t begin, size_t end) {
        for(size_t j=begin; j<end; ++j) {
          for(size_t i=0; i<query_span.size(); ++i) {
            auto res = pairwise::anchored_distance_affine(query_span[i], target_span[j], cost_map);
            output_ptr[i + j*query_span.size()] = std::get<0>(res);
            query_size_ptr[i + j*query_span.size()] = std::get<1>(res);
            target_size_ptr[i + j*query_span.size()] = std::get<2>(res);
          }
          progress_bar.increment();
        }
      }, 0, target_span.size(), 1, nthreads);
    }
    output.attr("query_size") = query_size;
    output.attr("target_size") = target_size;
  }
  return output;
}

// [[Rcpp::export(rng = false)]]
IntegerVector c_dist_pairwise(CharacterVector query, CharacterVector target, 
                            const std::string mode = "levenshtein", 
                            Rcpp::Nullable<IntegerMatrix> cost_matrix = R_NilValue,
                            const int gap_cost = 1,
                            const int gap_open_cost = 0,
                            const int nthreads = 1, const bool show_progress = false) {

  std::vector<cspan> query_span = strsxp_to_cspan(query);
  std::vector<cspan> target_span = strsxp_to_cspan(target);
  size_t nseqs = query_span.size();
  if(nseqs != target_span.size()) {
    throw std::runtime_error("Internal Error: query and target must be the same length");
  }

  trqwe::simple_progress progress_bar(nseqs, show_progress);
  IntegerVector output(nseqs);
  int * output_ptr = INTEGER(output);
  auto algo = decide_alignment_algo(mode, cost_matrix, gap_cost, gap_open_cost);
  if(algo == AlignmentAlgo::Hamming) {
    do_parallel_for([&query_span, &target_span, &progress_bar, output_ptr](size_t begin, size_t end) {
      for(size_t i=begin; i<end; ++i) {
        output_ptr[i] = pairwise::hamming_distance(query_span[i], target_span[i]);
        progress_bar.increment();
      }
    }, 0, nseqs, 1, nthreads);
  } else if(mode == "global" || mode == "gb" || mode == "lv" || mode == "levenshtein") {
    if(algo == AlignmentAlgo::GlobalUnit) {
      bool use_myers = false;
      CostMap cost_map;
      if (!cost_matrix.isNull()) {
        cost_map = convert_cost_matrix(cost_matrix.get(), gap_cost, gap_open_cost);
        use_myers = (gap_cost == 1 && gap_open_cost == 0);
        if (use_myers) {
          for (const auto & kv : cost_map.char_cost_map) { if (kv.second != 0 && kv.second != 1) { use_myers = false; break; } }
        }
      }
      if (use_myers) {
        do_parallel_for([&query_span, &target_span, &cost_map, &progress_bar, output_ptr](size_t begin, size_t end) {
          for(size_t i=begin; i<end; ++i) {
            output_ptr[i] = pairwise::global_distance_myers(query_span[i], target_span[i], cost_map);
            progress_bar.increment();
          }
        }, 0, nseqs, 1, nthreads);
      } else {
        do_parallel_for([&query_span, &target_span, &progress_bar, output_ptr](size_t begin, size_t end) {
          for(size_t i=begin; i<end; ++i) {
            output_ptr[i] = pairwise::global_distance(query_span[i], target_span[i]);
            progress_bar.increment();
          }
        }, 0, nseqs, 1, nthreads);
      }
    } else if(algo == AlignmentAlgo::GlobalLinear) {
      CostMap cost_map = convert_cost_matrix(cost_matrix.get(), gap_cost, gap_open_cost);
      bool boolean_subs = true;
      for (const auto & kv : cost_map.char_cost_map) {
        if (kv.second != 0 && kv.second != 1) { boolean_subs = false; break; }
      }
      const bool unit_gaps = (gap_cost == 1) && (gap_open_cost == 0);
      if (boolean_subs && unit_gaps) {
        do_parallel_for([&query_span, &target_span, &cost_map, &progress_bar, output_ptr](size_t begin, size_t end) {
          for(size_t i=begin; i<end; ++i) {
            output_ptr[i] = pairwise::global_distance_myers(query_span[i], target_span[i], cost_map);
            progress_bar.increment();
          }
        }, 0, nseqs, 1, nthreads);
      } else {
        do_parallel_for([&query_span, &target_span, &cost_map, &progress_bar, output_ptr](size_t begin, size_t end) {
          for(size_t i=begin; i<end; ++i) {
            output_ptr[i] = pairwise::global_distance_linear(query_span[i], target_span[i], cost_map);
            progress_bar.increment();
          }
        }, 0, nseqs, 1, nthreads);
      }
    } else { // GlobalAffine
      CostMap cost_map = convert_cost_matrix(cost_matrix.get(), gap_cost, gap_open_cost);
      do_parallel_for([&query_span, &target_span, &cost_map, &progress_bar, output_ptr](size_t begin, size_t end) {
        for(size_t i=begin; i<end; ++i) {
          output_ptr[i] = pairwise::global_distance_affine(query_span[i], target_span[i], cost_map);
          progress_bar.increment();
        }
      }, 0, nseqs, 1, nthreads);
    }
  } else { // anchored
    IntegerVector query_size(nseqs);
    IntegerVector target_size(nseqs);
    int * query_size_ptr = INTEGER(query_size);
    int * target_size_ptr = INTEGER(target_size);
    if(algo == AlignmentAlgo::AnchoredUnit) {
      bool use_myers = false;
      CostMap cost_map;
      if (!cost_matrix.isNull()) {
        cost_map = convert_cost_matrix(cost_matrix.get(), gap_cost, gap_open_cost);
        use_myers = (gap_cost == 1 && gap_open_cost == 0);
        if (use_myers) {
          for (const auto & kv : cost_map.char_cost_map) { if (kv.second != 0 && kv.second != 1) { use_myers = false; break; } }
        }
      }
      if (use_myers) {
        do_parallel_for([&query_span, &target_span, &cost_map, &progress_bar, output_ptr, query_size_ptr, target_size_ptr](size_t begin, size_t end) {
          for(size_t i=begin; i<end; ++i) {
            auto res = pairwise::anchored_distance_myers(query_span[i], target_span[i], cost_map);
            output_ptr[i] = std::get<0>(res);
            query_size_ptr[i] = std::get<1>(res);
            target_size_ptr[i] = std::get<2>(res);
            progress_bar.increment();
          }
        }, 0, nseqs, 1, nthreads);
      } else {
        do_parallel_for([&query_span, &target_span, &progress_bar, output_ptr, query_size_ptr, target_size_ptr](size_t begin, size_t end) {
          for(size_t i=begin; i<end; ++i) {
            auto res = pairwise::anchored_distance(query_span[i], target_span[i]);
            output_ptr[i] = std::get<0>(res);
            query_size_ptr[i] = std::get<1>(res);
            target_size_ptr[i] = std::get<2>(res);
            progress_bar.increment();
          }
        }, 0, nseqs, 1, nthreads);
      }
    } else if(algo == AlignmentAlgo::AnchoredLinear) {
      CostMap cost_map = convert_cost_matrix(cost_matrix.get(), gap_cost, gap_open_cost);
      bool boolean_subs = true;
      for (const auto & kv : cost_map.char_cost_map) { if (kv.second != 0 && kv.second != 1) { boolean_subs = false; break; } }
      const bool unit_gaps = (gap_cost == 1) && (gap_open_cost == 0);
      if (boolean_subs && unit_gaps) {
        do_parallel_for([&query_span, &target_span, &cost_map, &progress_bar, output_ptr, query_size_ptr, target_size_ptr](size_t begin, size_t end) {
          for(size_t i=begin; i<end; ++i) {
            auto res = pairwise::anchored_distance_myers(query_span[i], target_span[i], cost_map);
            output_ptr[i] = std::get<0>(res);
            query_size_ptr[i] = std::get<1>(res);
            target_size_ptr[i] = std::get<2>(res);
            progress_bar.increment();
          }
        }, 0, nseqs, 1, nthreads);
      } else {
        do_parallel_for([&query_span, &target_span, &cost_map, &progress_bar, output_ptr, query_size_ptr, target_size_ptr](size_t begin, size_t end) {
          for(size_t i=begin; i<end; ++i) {
            auto res = pairwise::anchored_distance_linear(query_span[i], target_span[i], cost_map);
            output_ptr[i] = std::get<0>(res);
            query_size_ptr[i] = std::get<1>(res);
            target_size_ptr[i] = std::get<2>(res);
            progress_bar.increment();
          }
        }, 0, nseqs, 1, nthreads);
      }
    } else { // AnchoredAffine
      CostMap cost_map = convert_cost_matrix(cost_matrix.get(), gap_cost, gap_open_cost);
      do_parallel_for([&query_span, &target_span, &cost_map, &progress_bar, output_ptr, query_size_ptr, target_size_ptr](size_t begin, size_t end) {
        for(size_t i=begin; i<end; ++i) {
          auto res = pairwise::anchored_distance_affine(query_span[i], target_span[i], cost_map);
          output_ptr[i] = std::get<0>(res);
          query_size_ptr[i] = std::get<1>(res);
          target_size_ptr[i] = std::get<2>(res);
          progress_bar.increment();
        }
      }, 0, nseqs, 1, nthreads);
    }
    output.attr("query_size") = query_size;
    output.attr("target_size") = target_size;
  }
  return output;
}
