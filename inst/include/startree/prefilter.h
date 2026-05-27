#ifndef STARTREE_PREFILTER_H
#define STARTREE_PREFILTER_H

#include "startree/common.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace startree {

class UnitCostPrefilter {
 public:
  UnitCostPrefilter() = default;

  void reset(const int median_len, const int padded_len, const int tau) {
    active_ = false;
    padded_len_ = padded_len;
    kmers_ = tau + 1;
    klen_.assign(static_cast<size_t>(kmers_), 0);
    bitmap_.clear();

    if(median_len < kmers_) {
      return;
    }

    const int k = median_len / kmers_;
    int rem = tau - median_len % kmers_;
    size_t total_bytes = 0;
    for(int i = 0; i < kmers_; ++i) {
      int len = k - (rem-- > 0 ? 1 : 0);
      if(len <= 0) {
        return;
      }
      len = std::min(len, kMaxLookupK);
      klen_[static_cast<size_t>(i)] = len;
      total_bytes += bitmap_bytes(len);
      if(total_bytes > kMaxLookupBytes) {
        return;
      }
    }

    bitmap_.reserve(static_cast<size_t>(kmers_));
    for(int len : klen_) {
      bitmap_.push_back(std::vector<uint8_t>(bitmap_bytes(len), 0));
    }
    active_ = true;
  }

  void insert(const std::string& padded_seq) {
    if(!active_) {
      return;
    }

    int offset = padded_len_;
    for(int i = kmers_ - 1; i >= 0; --i) {
      const int len = klen_[static_cast<size_t>(i)];
      offset -= len;
      const int seqid = seq2id(padded_seq, offset, len);
      if(seqid >= 0) {
        std::vector<uint8_t>& bits = bitmap_[static_cast<size_t>(i)];
        bits[static_cast<size_t>(seqid / 8)] |=
          static_cast<uint8_t>(1U << (seqid % 8));
      }
    }
  }

  bool maybe_contains(const std::string& padded_query) const {
    if(!active_) {
      return true;
    }

    int offset = padded_len_;
    for(int i = kmers_ - 1; i >= 0; --i) {
      const int len = klen_[static_cast<size_t>(i)];
      offset -= len;
      const int max_shift = kmers_ - 1 - i;
      for(int shift = -max_shift; shift <= max_shift; ++shift) {
        const int seqid = seq2id(padded_query, offset + shift, len);
        if(seqid < 0) {
          continue;
        }
        const std::vector<uint8_t>& bits = bitmap_[static_cast<size_t>(i)];
        if((bits[static_cast<size_t>(seqid / 8)] >> (seqid % 8)) & 1U) {
          return true;
        }
      }
    }
    return false;
  }

 private:
  static size_t bitmap_bytes(const int kmer_len) {
    const int shift = std::max(0, 2 * kmer_len - 3);
    return size_t{1} << shift;
  }

  static int seq2id(const std::string& seq, const int start, const int len) {
    if(start < 0 || len < 0 || start + len > static_cast<int>(seq.size())) {
      return -2;
    }

    static constexpr int kCodeToLookup[6] = {1, 0, 1, 2, 3, 0};
    int seqid = 0;
    const int begin = start + std::max(0, len - 16);
    for(int i = begin; i < start + len; ++i) {
      const unsigned char c = static_cast<unsigned char>(seq[static_cast<size_t>(i)]);
      seqid += kCodeToLookup[c];
      if(i < start + len - 1) {
        seqid <<= 2;
      }
    }
    return seqid;
  }

  bool active_ = false;
  int padded_len_ = 0;
  int kmers_ = 0;
  std::vector<int> klen_;
  std::vector<std::vector<uint8_t>> bitmap_;
};

class WeightedCostPrefilter {
 public:
  WeightedCostPrefilter() = default;

  void reset(const int median_len,
             const int padded_len,
             const int lookup_distance,
             const int max_distance,
             const int min_cost,
             const int gap_cost) {
    active_ = false;
    padded_len_ = padded_len;
    kmers_ = lookup_distance + 1;
    klen_.assign(static_cast<size_t>(kmers_), 0);
    max_shift_.assign(static_cast<size_t>(kmers_), 0);
    bitmap_.clear();

    if(median_len < kmers_) {
      return;
    }

    const int k = median_len / kmers_;
    int rem = lookup_distance - median_len % kmers_;
    size_t total_bytes = 0;
    for(int i = 0; i < kmers_; ++i) {
      int len = k - (rem-- > 0 ? 1 : 0);
      if(len <= 0) {
        return;
      }
      len = std::min(len, kMaxLookupK);
      klen_[static_cast<size_t>(i)] = len;
      max_shift_[static_cast<size_t>(i)] =
        max_lookup_shift(kmers_ - 1 - i, max_distance, min_cost, gap_cost);
      total_bytes += bitmap_bytes(len);
      if(total_bytes > kMaxLookupBytes) {
        return;
      }
    }

    bitmap_.reserve(static_cast<size_t>(kmers_));
    for(int len : klen_) {
      bitmap_.push_back(std::vector<uint8_t>(bitmap_bytes(len), 0));
    }
    active_ = true;
  }

  void insert(const std::string& padded_seq) {
    if(!active_) {
      return;
    }

    int offset = padded_len_;
    for(int i = kmers_ - 1; i >= 0; --i) {
      const int len = klen_[static_cast<size_t>(i)];
      offset -= len;
      const int seqid = seq2id(padded_seq, offset, len);
      if(seqid >= 0) {
        std::vector<uint8_t>& bits = bitmap_[static_cast<size_t>(i)];
        bits[static_cast<size_t>(seqid / 8)] |=
          static_cast<uint8_t>(1U << (seqid % 8));
      }
    }
  }

  bool maybe_contains(const std::string& padded_query) const {
    if(!active_) {
      return true;
    }

    int offset = padded_len_;
    for(int i = kmers_ - 1; i >= 0; --i) {
      const int len = klen_[static_cast<size_t>(i)];
      offset -= len;
      const int max_shift = max_shift_[static_cast<size_t>(i)];
      for(int shift = -max_shift; shift <= max_shift; ++shift) {
        const int seqid = seq2id(padded_query, offset + shift, len);
        if(seqid < 0) {
          continue;
        }
        const std::vector<uint8_t>& bits = bitmap_[static_cast<size_t>(i)];
        if((bits[static_cast<size_t>(seqid / 8)] >> (seqid % 8)) & 1U) {
          return true;
        }
      }
    }
    return false;
  }

 private:
  static size_t bitmap_bytes(const int kmer_len) {
    const int shift = std::max(0, 2 * kmer_len - 3);
    return size_t{1} << shift;
  }

  static int max_lookup_shift(const int chunks_to_right,
                              const int max_distance,
                              const int min_cost,
                              const int gap_cost) {
    int max_shift = 0;
    for(int shift = 0; shift <= chunks_to_right; ++shift) {
      const int remaining_chunks = chunks_to_right - shift;
      const int64_t min_required_cost =
        static_cast<int64_t>(shift) * gap_cost +
        static_cast<int64_t>(remaining_chunks) * min_cost;
      if(min_required_cost <= max_distance) {
        max_shift = shift;
      }
    }
    return max_shift;
  }

  static int seq2id(const std::string& seq, const int start, const int len) {
    if(start < 0 || len < 0 || start + len > static_cast<int>(seq.size())) {
      return -2;
    }

    static constexpr int kCodeToLookup[6] = {1, 0, 1, 2, 3, 0};
    int seqid = 0;
    const int begin = start + std::max(0, len - 16);
    for(int i = begin; i < start + len; ++i) {
      const unsigned char c = static_cast<unsigned char>(seq[static_cast<size_t>(i)]);
      seqid += kCodeToLookup[c];
      if(i < start + len - 1) {
        seqid <<= 2;
      }
    }
    return seqid;
  }

  bool active_ = false;
  int padded_len_ = 0;
  int kmers_ = 0;
  std::vector<int> klen_;
  std::vector<int> max_shift_;
  std::vector<std::vector<uint8_t>> bitmap_;
};

}  // namespace startree

#endif  // STARTREE_PREFILTER_H
