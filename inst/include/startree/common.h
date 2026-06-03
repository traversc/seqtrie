#ifndef STARTREE_COMMON_H
#define STARTREE_COMMON_H

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace startree {

constexpr int kMaxTau = 8;
constexpr int kPad = 5;
constexpr int kEos = -1;
constexpr int kNoMatch = 6;
constexpr int kAlphabet = 6;
constexpr int kMaxSeqLen = 1023;
constexpr uint32_t kNoSequence = std::numeric_limits<uint32_t>::max();
constexpr int kAutoHeight = -1;

struct Sequence {
  std::string seq;
  std::string padded;
  int64_t count = 1;
  uint32_t min_input_id = 0;
};

struct InputData {
  std::vector<Sequence> seqs;
  int median_len = 0;
  int height = 0;
};

struct AnchoredInputData {
  std::vector<Sequence> seqs;
  std::vector<std::string> target_codes;
  std::vector<std::string> query_codes;
};

struct Hit {
  uint32_t sequence_index = kNoSequence;
  int distance = 0;
};

struct PairRecord {
  uint32_t target_index = 0;
  uint32_t query_index = 0;
  int distance = 0;
};

struct Block {
  uint32_t begin = 0;
  uint32_t end = 0;
  size_t min_len = 0;
  size_t max_len = 0;
};

struct SearchParams {
  int max_distance = 0;
  int band_radius = 0;
  int lookup_distance = 0;
  bool custom_cost = false;
  bool include_zero = false;
  int mismatch_cost = 1;
  int gap_cost = 1;
};

inline int insert_code(const unsigned char c) noexcept {
  switch(c) {
    case 'N':
    case 'n':
      return 0;
    case 'A':
    case 'a':
      return 1;
    case 'C':
    case 'c':
      return 2;
    case 'G':
    case 'g':
      return 3;
    case 'T':
    case 't':
      return 4;
    case ' ':
      return kPad;
    default:
      return 0;
  }
}

inline bool is_valid_dna(const unsigned char c) noexcept {
  switch(c) {
    case 'A':
    case 'a':
    case 'C':
    case 'c':
    case 'G':
    case 'g':
    case 'T':
    case 't':
    case 'N':
    case 'n':
      return true;
    default:
      return false;
  }
}

inline char upper_base(const unsigned char c) noexcept {
  return static_cast<char>(std::toupper(c));
}

inline void normalize_sequence_in_place(std::string& seq) {
  for(char& ch : seq) {
    ch = upper_base(static_cast<unsigned char>(ch));
  }
}

inline std::string make_anchored_target_code(const std::string& seq) {
  std::string out;
  out.reserve(seq.size());
  for(char ch : seq) {
    out.push_back(static_cast<char>(
      insert_code(static_cast<unsigned char>(ch))
    ));
  }
  return out;
}

inline std::string make_anchored_query_code(const std::string& seq) {
  std::string out;
  out.reserve(seq.size());
  for(char ch : seq) {
    if(ch == 'N' || ch == 'n') {
      out.push_back(static_cast<char>(kNoMatch));
    } else {
      out.push_back(static_cast<char>(
        insert_code(static_cast<unsigned char>(ch))
      ));
    }
  }
  return out;
}

inline std::vector<Block> make_blocks(const std::vector<Sequence>& seqs,
                                      const int threads) {
  const size_t thread_count = static_cast<size_t>(threads);
  size_t nblocks = 3 * thread_count + (threads % 2 == 0 ? 1ULL : 0ULL);
  if(seqs.size() < nblocks) {
    nblocks = 1;
  }

  std::vector<Block> blocks;
  blocks.reserve(nblocks);
  const size_t q = seqs.size() / nblocks;
  const size_t r = seqs.size() % nblocks;
  for(size_t i = 0; i < nblocks; ++i) {
    const size_t begin = q * i + std::min(i, r);
    const size_t end = q * (i + 1) + std::min(i + 1, r);
    if(begin == end) {
      continue;
    }
    Block block;
    block.begin = static_cast<uint32_t>(begin);
    block.end = static_cast<uint32_t>(end);
    block.min_len = seqs[begin].seq.size();
    block.max_len = seqs[end - 1].seq.size();
    blocks.push_back(block);
  }
  return blocks;
}

inline bool length_ranges_can_match(const Block& target,
                                    const Block& query,
                                    const int max_distance,
                                    const int gap_cost) {
  const size_t max_gap_bases = static_cast<size_t>(max_distance / gap_cost);
  if(target.min_len > query.max_len) {
    return target.min_len - query.max_len <= max_gap_bases;
  }
  if(query.min_len > target.max_len) {
    return query.min_len - target.max_len <= max_gap_bases;
  }
  return true;
}

inline int common_prefix_length(const std::string& a, const std::string& b) {
  const size_t n = std::min(a.size(), b.size());
  size_t i = 0;
  while(i < n && a[i] == b[i]) {
    ++i;
  }
  return static_cast<int>(i);
}

inline bool sequence_length_can_match(const Block& target,
                                      const size_t query_len,
                                      const int max_distance,
                                      const int gap_cost) {
  const size_t max_gap_bases = static_cast<size_t>(max_distance / gap_cost);
  if(target.min_len > query_len) {
    return target.min_len - query_len <= max_gap_bases;
  }
  if(query_len > target.max_len) {
    return query_len - target.max_len <= max_gap_bases;
  }
  return true;
}

inline SearchParams make_search_params(const int max_distance,
                                       const int mismatch_cost,
                                       const int gap_cost,
                                       const bool include_zero = false) {
  SearchParams params;
  params.max_distance = max_distance;
  params.include_zero = include_zero;
  params.mismatch_cost = mismatch_cost;
  params.gap_cost = gap_cost;
  params.custom_cost = mismatch_cost != 1 || gap_cost != 1;

  if(params.custom_cost) {
    const int min_cost = std::min(mismatch_cost, gap_cost);
    params.lookup_distance = max_distance / min_cost;
    params.band_radius = max_distance / gap_cost;
  } else {
    params.lookup_distance = max_distance;
    params.band_radius = max_distance;
  }

  return params;
}

inline InputData make_input_data_from_sequences(std::vector<Sequence> processed,
                                                const bool deduplicate,
                                                const int height_override = kAutoHeight) {
  int height = 0;
  for(Sequence& seq : processed) {
    normalize_sequence_in_place(seq.seq);
    height = std::max(height, static_cast<int>(seq.seq.size()));
  }

  std::sort(processed.begin(), processed.end(), [](const Sequence& a, const Sequence& b) {
    if(a.seq.size() != b.seq.size()) {
      return a.seq.size() < b.seq.size();
    }
    return a.seq < b.seq;
  });

  if(deduplicate) {
    std::vector<Sequence> unique;
    unique.reserve(processed.size());
    for(Sequence& seq : processed) {
      if(!unique.empty() && unique.back().seq == seq.seq) {
        unique.back().count += seq.count;
        unique.back().min_input_id = std::min(unique.back().min_input_id,
                                              seq.min_input_id);
      } else {
        unique.push_back(std::move(seq));
      }
    }
    processed = std::move(unique);
  }

  if(height_override != kAutoHeight) {
    height = std::max(height, height_override);
  }

  std::vector<size_t> length_counts(static_cast<size_t>(height + 1), 0);
  for(Sequence& seq : processed) {
    const int len = static_cast<int>(seq.seq.size());
    ++length_counts[static_cast<size_t>(len)];
    seq.padded.assign(static_cast<size_t>(height - len), static_cast<char>(kPad));
    seq.padded.reserve(static_cast<size_t>(height));
    for(char ch : seq.seq) {
      seq.padded.push_back(static_cast<char>(
        insert_code(static_cast<unsigned char>(ch))
      ));
    }
  }

  int median = 0;
  size_t cumulative = 0;
  do {
    ++median;
    cumulative += length_counts[static_cast<size_t>(median)];
  } while(median < height && cumulative < processed.size() / 2);

  InputData data;
  data.seqs = std::move(processed);
  data.median_len = median;
  data.height = height;
  return data;
}

inline InputData make_input_data(std::vector<std::string> seqs,
                                 const bool deduplicate,
                                 const int height_override = kAutoHeight) {
  std::vector<Sequence> processed;
  processed.reserve(seqs.size());
  uint32_t input_id = 0;

  for(std::string& seq : seqs) {
    Sequence item;
    item.seq = std::move(seq);
    item.min_input_id = ++input_id;
    processed.push_back(std::move(item));
  }

  return make_input_data_from_sequences(std::move(processed), deduplicate,
                                        height_override);
}

inline AnchoredInputData make_anchored_input_data_from_sequences(
    std::vector<Sequence> processed,
    const bool deduplicate,
    const bool sort_by_query_code = false) {
  struct Item {
    Sequence seq;
    std::string target_code;
    std::string query_code;
  };

  std::vector<Item> items;
  items.reserve(processed.size());
  for(Sequence& seq : processed) {
    normalize_sequence_in_place(seq.seq);
    Item item;
    item.target_code = make_anchored_target_code(seq.seq);
    item.query_code = make_anchored_query_code(seq.seq);
    item.seq = std::move(seq);
    items.push_back(std::move(item));
  }

  std::sort(items.begin(), items.end(),
            [sort_by_query_code](const Item& a, const Item& b) {
    const std::string& a_key = sort_by_query_code ? a.query_code : a.target_code;
    const std::string& b_key = sort_by_query_code ? b.query_code : b.target_code;
    if(a_key != b_key) {
      return a_key < b_key;
    }
    return a.seq.seq < b.seq.seq;
  });

  if(deduplicate) {
    std::vector<Item> unique;
    unique.reserve(items.size());
    for(Item& item : items) {
      if(!unique.empty() && unique.back().seq.seq == item.seq.seq) {
        unique.back().seq.count += item.seq.count;
        unique.back().seq.min_input_id =
          std::min(unique.back().seq.min_input_id, item.seq.min_input_id);
      } else {
        unique.push_back(std::move(item));
      }
    }
    items = std::move(unique);
  }

  AnchoredInputData data;
  data.seqs.reserve(items.size());
  data.target_codes.reserve(items.size());
  data.query_codes.reserve(items.size());
  for(Item& item : items) {
    data.target_codes.push_back(std::move(item.target_code));
    data.query_codes.push_back(std::move(item.query_code));
    data.seqs.push_back(std::move(item.seq));
  }
  return data;
}

inline AnchoredInputData make_anchored_input_data(
    std::vector<std::string> seqs,
    const bool deduplicate,
    const bool sort_by_query_code = false) {
  std::vector<Sequence> processed;
  processed.reserve(seqs.size());
  uint32_t input_id = 0;

  for(std::string& seq : seqs) {
    Sequence item;
    item.seq = std::move(seq);
    item.min_input_id = ++input_id;
    processed.push_back(std::move(item));
  }

  return make_anchored_input_data_from_sequences(
    std::move(processed), deduplicate, sort_by_query_code
  );
}

}  // namespace startree

#endif  // STARTREE_COMMON_H
