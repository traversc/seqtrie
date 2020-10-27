#include <Rcpp.h>
#include <unordered_map>
#include <memory>
#include <tbb/concurrent_vector.h>
#include <boost/functional/hash.hpp>
#include <RcppParallel.h>
#include <tuple>
#include <limits.h> // INT_MAX

#include "treedist_types.h"

using namespace Rcpp;
using namespace RcppParallel;

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppParallel)]]

inline void parallelFor2(std::size_t begin, std::size_t end, Worker& worker, std::size_t grainSize = 1, int nthreads = 1) {
  int max_threads = tbb::task_scheduler_init::default_num_threads();
  if(nthreads > max_threads) nthreads = max_threads;
  tbb::task_arena limited(nthreads);
  tbb::task_group tg;
  limited.execute([&]{
    tg.run([&]{
      parallelFor(begin, end, worker, grainSize);
    });
  });
  limited.execute([&]{ tg.wait(); });
}

std::vector<int> int_range(int start, int end) {
  if(end < start) {
    throw std::runtime_error("end must be larger or equal start");
  }
  std::vector<int> ret(end - start + 1);
  for(int i=0; i<end - start + 1; ++i) {
    ret[i] = start + i;
  }
  return ret;
}

void update_prefix_tree(SeqNode * root, const std::string & sequence, const int seq_idx) {
  SeqNode * current_hash_ref = root;
  size_t i = 0;
  for(; i<sequence.size(); ++i) {
    if(current_hash_ref->seqmap.find(sequence[i]) == current_hash_ref->seqmap.end()) {
      current_hash_ref->seqmap.insert(std::pair<char, SeqNodeUptr >(sequence[i], std::make_unique<SeqNode>() ));
    }
    current_hash_ref = current_hash_ref->seqmap[sequence[i]].get();
  }
  current_hash_ref->idx.insert(seq_idx);
}

void update_prefix_tree(SeqNode * root, const char * const sequence, const size_t seq_len, const int seq_idx) {
  SeqNode * current_hash_ref = root;
  size_t i = 0;
  for(; i<seq_len; ++i) {
    if(current_hash_ref->seqmap.find(sequence[i]) == current_hash_ref->seqmap.end()) {
      current_hash_ref->seqmap.insert(std::pair<char, SeqNodeUptr >(sequence[i], std::make_unique<SeqNode>() ));
    }
    current_hash_ref = current_hash_ref->seqmap[sequence[i]].get();
  }
  current_hash_ref->idx.insert(seq_idx);
}

// same as above
void update_suffix_tree(SeqNode * root, const char * const sequence, const size_t seq_len, const int seq_idx) {
  SeqNode * current_hash_ref = root;
  size_t i = 0;
  for(; i<seq_len; ++i) {
    if(current_hash_ref->seqmap.find(sequence[i]) == current_hash_ref->seqmap.end()) {
      current_hash_ref->seqmap.insert(std::pair<char, SeqNodeUptr >(sequence[i], std::make_unique<SeqNode>() ));
    }
    current_hash_ref = current_hash_ref->seqmap[sequence[i]].get();
  }
  current_hash_ref->idx.insert(seq_idx);
}

/////////////////////////////////////////////////////////

// [[Rcpp::export(rng = false)]]
Rcpp::XPtr<SeqNode> c_td_prefix_tree(const std::vector<std::string> & sequences) {
  auto * root = new SeqNode;
  for(size_t i=0; i<sequences.size(); ++i) {
    update_prefix_tree(root, sequences[i], i);
  }
  return Rcpp::XPtr<SeqNode> (root, true);
}

// [[Rcpp::export(rng = false)]]
Rcpp::XPtr<SeqNode> c_td_suffix_tree(const std::vector<std::string> & sequences, size_t min_length) {
  auto * root = new SeqNode;
  for(size_t i=0; i<sequences.size(); ++i) {
    for(size_t j=0; sequences[i].size() - j >= min_length; ++j) {
      update_suffix_tree(root, sequences[i].c_str() + j, sequences[i].size() - j, i);
    }
  }
  return Rcpp::XPtr<SeqNode> (root, true);
}

///////////////////////////////////////////////////////////////
// hamming

struct HammingWorker : public Worker {
  SeqNode const * const root;
  tbb::concurrent_vector<std::tuple<int, int, int>> & output;
  const std::vector<std::string> & sequences;
  const int max_distance;
  const bool symmetric;
  HammingWorker(SeqNode const * const root, 
                    tbb::concurrent_vector<std::tuple<int, int, int>> & output,
                    const std::vector<std::string> & sequences, 
                    const int max_distance, const bool symmetric) :
    root(root), output(output), sequences(sequences), max_distance(max_distance), symmetric(symmetric) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t i=begin; i<end; ++i) {
      search(root, i, sequences[i].c_str(), 0);
    }
  }
  void search(const SeqNode * node, const int seqidx, const char * sequence, const int distance) {
    // we've exceeded max_distance
    if( distance > max_distance ) {
      return;
    }
    // we've reached the end
    if( (*sequence == 0) || node->isLeaf() ) {
      for(auto i : node->idx) {
        if( (distance <= max_distance) && ((i < seqidx) || !symmetric) ) {
          output.push_back(std::make_tuple(seqidx, i, distance));
        }
      }
      return;
    }
    auto & seqmap = node->seqmap;
    for (auto & x : seqmap) {
      if(*sequence == x.first) {
        search(x.second.get(), seqidx, sequence+1, distance);
      } else {
        search(x.second.get(), seqidx, sequence+1, distance + 1);
      }
    }
  }
};

// [[Rcpp::export(rng = false)]]
DataFrame c_td_hamming(Rcpp::XPtr<SeqNode> tree, 
                           const std::vector<std::string> & sequences, 
                           const int max_distance, const bool symmetric, const int nthreads) {
  tbb::concurrent_vector<std::tuple<int, int, int>> output;
  HammingWorker w(tree.get(), output, sequences, max_distance, symmetric);
  parallelFor2(symmetric ? 1 : 0, sequences.size(), w, 1, nthreads);
  IntegerVector query(output.size());
  IntegerVector target(output.size());
  IntegerVector distance(output.size());
  for(size_t i=0; i<output.size(); ++i) {
    query[i] = std::get<0>(output[i]) + 1;
    target[i] = std::get<1>(output[i]) + 1;
    distance[i] = std::get<2>(output[i]);
  }
  return DataFrame::create(_["query"] = query, _["target"] = target, _["distance"] = distance);
}

///////////////////////////////////////////////////////////////
// levenshtein

struct LevenshteinWorker : public Worker {
  SeqNode const * const root;
  tbb::concurrent_vector<std::tuple<int, int, int>> & output;
  const std::vector<std::string> & sequences;
  const int max_distance;
  const bool symmetric;
  LevenshteinWorker(SeqNode const * const root, 
                    tbb::concurrent_vector<std::tuple<int, int, int>> & output,
                    const std::vector<std::string> & sequences, 
                    const int max_distance, const bool symmetric) :
    root(root), output(output), sequences(sequences), max_distance(max_distance), symmetric(symmetric) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t i=begin; i<end; ++i) {
      search(root, i, int_range(0, sequences[i].size()));
    }
  }
  void search(SeqNode const * const node, const int seqidx, const std::vector<int> previous_row) {
    if( *std::min_element(previous_row.begin(), previous_row.end()) > max_distance ) {
      return;
    }
    for(auto i : node->idx) {
      int dist = previous_row.back();
      if( (dist <= max_distance) && ((i < seqidx) || !symmetric) ) {
        output.push_back(std::make_tuple(seqidx, i, dist));
      }
    }
    if(node->isLeaf()) {
      return;
    }
    auto & seqmap = node->seqmap;
    for (auto & x : seqmap) {
      std::vector<int> current_row(sequences[seqidx].size() + 1);
      for(size_t i=0; i<current_row.size(); ++i) {
        int delete_cost = previous_row[i] + 1;
        if(i > 0) {
          int match_cost  = previous_row[i-1] + (sequences[seqidx][i-1] == x.first ? 0 : 1);
          int insert_cost = current_row[i-1] + 1;
          current_row[i] = std::min(std::min(delete_cost, match_cost), insert_cost);
        } else {
          current_row[i] = delete_cost;
        }
      }
      search(x.second.get(), seqidx, current_row);
    }
  }
};

// [[Rcpp::export(rng = false)]]
DataFrame c_td_levenshtein(Rcpp::XPtr<SeqNode> tree, 
                           const std::vector<std::string> & sequences, 
                           const int max_distance, const bool symmetric, const int nthreads) {
  tbb::concurrent_vector<std::tuple<int, int, int>> output;
  LevenshteinWorker w(tree.get(), output, sequences, max_distance, symmetric);
  parallelFor2(symmetric ? 1 : 0, sequences.size(), w, 1, nthreads);
  IntegerVector query(output.size());
  IntegerVector target(output.size());
  IntegerVector distance(output.size());
  for(size_t i=0; i<output.size(); ++i) {
    query[i] = std::get<0>(output[i]) + 1;
    target[i] = std::get<1>(output[i]) + 1;
    distance[i] = std::get<2>(output[i]);
  }
  return DataFrame::create(_["query"] = query, _["target"] = target, _["distance"] = distance);
}

///////////////////////////////////////////////////////////////
// partial hamming

struct PartialHammingWorker : public Worker {
  SeqNode const * const root;
  tbb::concurrent_vector<std::tuple<int, int, int>> & output;
  const std::vector<std::string> & sequences;
  const int max_distance;
  const bool symmetric;
  PartialHammingWorker(SeqNode const * const root, 
                       tbb::concurrent_vector<std::tuple<int, int, int>> & output,
                       const std::vector<std::string> & sequences, 
                       const int max_distance, const bool symmetric) :
    root(root), output(output), sequences(sequences), max_distance(max_distance), symmetric(symmetric) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t i=begin; i<end; ++i) {
      search(root, i, sequences[i].c_str(), 0, true);
    }
  }
  void search(const SeqNode * node, const int seqidx, const char * sequence, const int distance, bool opening) {
    // we've exceeded max_distance
    if( distance > max_distance ) {
      return;
    }
    // we've reached the end of the sequence
    if( (*sequence == 0) ) {
      std::set<int> idxs = node->allChildIdx();
      for(auto i : idxs) {
        if( (distance <= max_distance) && ((i < seqidx) || !symmetric) ) {
          output.push_back(std::make_tuple(seqidx, i, distance));
        }
      }
      return;
    // we've reached the end of the tree without finishing the sequence
    } else if(node->isLeaf()) {
      return;
    }
    auto & seqmap = node->seqmap;
    for (auto & x : seqmap) {
      if(*sequence == x.first) {
        search(x.second.get(), seqidx, sequence+1, distance, false);
      } else {
        search(x.second.get(), seqidx, sequence+1, distance+1, false);
      }
      if(opening) {
        search(x.second.get(), seqidx, sequence, distance, true); // search next tree element without incrementing sequence; opening gap
      }
    }
  }
};

// [[Rcpp::export(rng = false)]]
DataFrame c_td_partial_hamming(Rcpp::XPtr<SeqNode> tree, 
                               const std::vector<std::string> & sequences, 
                               const int max_distance, const bool symmetric, const int nthreads) {
  tbb::concurrent_vector<std::tuple<int, int, int>> output;
  PartialHammingWorker w(tree.get(), output, sequences, max_distance, symmetric);
  parallelFor2(symmetric ? 1 : 0, sequences.size(), w, 1, nthreads);
  IntegerVector query(output.size());
  IntegerVector target(output.size());
  IntegerVector distance(output.size());
  for(size_t i=0; i<output.size(); ++i) {
    query[i] = std::get<0>(output[i]) + 1;
    target[i] = std::get<1>(output[i]) + 1;
    distance[i] = std::get<2>(output[i]);
  }
  return DataFrame::create(_["query"] = query, _["target"] = target, _["distance"] = distance);
}

///////////////////////////////////////////////////////////////
// Partial/anchored levenshtein
enum class Direction {left, right};
struct PartialLevenshteinWorker : public Worker {
  SeqNode const * const root;
  tbb::concurrent_vector<std::tuple<int, int, int>> & output;
  const std::vector<std::string> & sequences;
  const Direction anchor;
  const int max_distance;
  const bool symmetric;
  PartialLevenshteinWorker(SeqNode const * const root, 
                    tbb::concurrent_vector<std::tuple<int, int, int>> & output,
                    const std::vector<std::string> & sequences, 
                    const std::string & anchor,
                    const int max_distance, const bool symmetric) :
    root(root), output(output), sequences(sequences), anchor(get_direction(anchor)), max_distance(max_distance), symmetric(symmetric) {}
  static Direction get_direction(const std::string & a) {
    if(a == "left") {
      return Direction::left;
    } else if(a == "right") {
      return Direction::right;
    } else {
      throw std::runtime_error("anchor must be left or right");
    }
  }
  void operator()(std::size_t begin, std::size_t end) {
    for(size_t i=begin; i<end; ++i) {
      if(anchor == Direction::left) {
        search(root, i, int_range(0, sequences[i].size()));
      } else {
        search(root, i, std::vector<int>(sequences[i].size(),0) );
      }
    }
  }
  void search(SeqNode const * const node, const int seqidx, const std::vector<int> previous_row) {
    if( *std::min_element(previous_row.begin(), previous_row.end()) > max_distance ) {
      return;
    }
    for(auto i : node->idx) {
      int dist;
      if(anchor == Direction::right) {
        dist = previous_row.back();
      } else {
        dist = *std::min_element(previous_row.begin(), previous_row.end());
      }

      if( (dist <= max_distance) && ((i < seqidx) || !symmetric) ) {
        output.push_back(std::make_tuple(seqidx, i, dist));
      }
    }
    if(node->isLeaf()) {
      return;
    }
    auto & seqmap = node->seqmap;
    for (auto & x : seqmap) {
      std::vector<int> current_row(sequences[seqidx].size() + 1);
      for(size_t i=0; i<current_row.size(); ++i) {
        if(i == 0) {
          if(anchor == Direction::left) {
            current_row[i] = previous_row[i] + 1;
          } else {
            current_row[i] = 0;
          }
        } else if(i == current_row.size() - 1) { // last column
          int delete_cost;
          if(anchor == Direction::right) {
            delete_cost = previous_row[i] + 1;
          } else {
            delete_cost = previous_row[i];
          }
          int match_cost  = previous_row[i-1] + (sequences[seqidx][i-1] == x.first ? 0 : 1);
          int insert_cost = current_row[i-1] + 1;
          current_row[i] = std::min({delete_cost, match_cost, insert_cost});
        } else {
          int delete_cost = previous_row[i] + 1;
          int match_cost  = previous_row[i-1] + (sequences[seqidx][i-1] == x.first ? 0 : 1);
          int insert_cost = current_row[i-1] + 1;
          current_row[i] = std::min({delete_cost, match_cost, insert_cost});
        }
      }
      search(x.second.get(), seqidx, current_row);
    }
  }
};

// [[Rcpp::export(rng = false)]]
DataFrame c_td_partial_levenshtein(Rcpp::XPtr<SeqNode> tree, 
                           const std::vector<std::string> & sequences, const std::string anchor,
                           const int max_distance, const bool symmetric, const int nthreads) {
  tbb::concurrent_vector<std::tuple<int, int, int>> output;
  PartialLevenshteinWorker w(tree.get(), output, sequences, anchor, max_distance, symmetric);
  parallelFor2(symmetric ? 1 : 0, sequences.size(), w, 1, nthreads);
  IntegerVector query(output.size());
  IntegerVector target(output.size());
  IntegerVector distance(output.size());
  for(size_t i=0; i<output.size(); ++i) {
    query[i] = std::get<0>(output[i]) + 1;
    target[i] = std::get<1>(output[i]) + 1;
    distance[i] = std::get<2>(output[i]);
  }
  return DataFrame::create(_["query"] = query, _["target"] = target, _["distance"] = distance);
}
