#ifndef STARTREE_PREFILTER_H
#define STARTREE_PREFILTER_H

#include "startree/common.h"
#include "ankerl/unordered_dense.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace startree {

namespace prefilter_detail {

// The segment key is a wyhash of the raw padded codes, so it is already
// avalanching. Hand it to ankerl unchanged -- the passthrough hash prevents a
// redundant second mix while still letting ankerl select buckets from the high
// bits.
struct PassthroughHash {
  using is_avalanching = void;
  uint64_t operator()(const uint64_t x) const noexcept { return x; }
};

using KeySet = ankerl::unordered_dense::set<uint64_t, PassthroughHash>;

// Key for seq[start, start + len): wyhash over the full segment of raw padded
// codes. Returns false when the window runs off the string. Padding (kPad) and
// N are their own symbols -- no folding into ACGT is needed (that was only ever
// a consequence of the old 2-bit bitmap slot). Correctness requires only that
// insert and probe derive the same key from identical content, which holds
// because both call this on the canonical padded sequence.
inline bool segment_key(const std::string& seq, const int start, const int len,
                        uint64_t* out) {
  if(start < 0 || len < 0 ||
     start + len > static_cast<int>(seq.size())) {
    return false;
  }
  *out = ankerl::unordered_dense::detail::wyhash::hash(
    seq.data() + static_cast<size_t>(start), static_cast<size_t>(len));
  return true;
}

}  // namespace prefilter_detail

// Pigeonhole k-mer prefilter for unit-cost (Levenshtein) search. The padded
// sequence is split into tau + 1 contiguous segments; a target within tau of a
// query must share at least one segment (allowing for indel drift, covered by
// the shifted probes). Each segment is keyed by wyhash and stored in a hash set,
// so it stays active for any segment length up to kMaxSeqLen.
class UnitCostPrefilter {
 public:
  UnitCostPrefilter() = default;

  void reset(const int median_len, const int padded_len, const int tau) {
    active_ = false;
    padded_len_ = padded_len;
    kmers_ = tau + 1;
    klen_.assign(static_cast<size_t>(kmers_), 0);
    sets_.clear();

    if(median_len < kmers_) {
      return;
    }

    const int k = median_len / kmers_;
    int rem = tau - median_len % kmers_;
    for(int i = 0; i < kmers_; ++i) {
      const int len = k - (rem-- > 0 ? 1 : 0);
      if(len <= 0) {
        return;
      }
      klen_[static_cast<size_t>(i)] = len;
    }

    sets_.resize(static_cast<size_t>(kmers_));
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
      uint64_t id = 0;
      if(prefilter_detail::segment_key(padded_seq, offset, len, &id)) {
        sets_[static_cast<size_t>(i)].insert(id);
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
      const prefilter_detail::KeySet& set = sets_[static_cast<size_t>(i)];
      for(int shift = -max_shift; shift <= max_shift; ++shift) {
        uint64_t id = 0;
        if(prefilter_detail::segment_key(padded_query, offset + shift, len, &id) &&
           set.find(id) != set.end()) {
          return true;
        }
      }
    }
    return false;
  }

 private:
  bool active_ = false;
  int padded_len_ = 0;
  int kmers_ = 0;
  std::vector<int> klen_;
  std::vector<prefilter_detail::KeySet> sets_;
};

// Pigeonhole prefilter for weighted (custom-cost) search. Same segment scheme,
// but the number of segments and the per-segment probe shift are derived from
// the lookup distance and the mismatch/gap costs rather than tau directly.
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
    sets_.clear();

    if(median_len < kmers_) {
      return;
    }

    const int k = median_len / kmers_;
    int rem = lookup_distance - median_len % kmers_;
    for(int i = 0; i < kmers_; ++i) {
      const int len = k - (rem-- > 0 ? 1 : 0);
      if(len <= 0) {
        return;
      }
      klen_[static_cast<size_t>(i)] = len;
      max_shift_[static_cast<size_t>(i)] =
        max_lookup_shift(kmers_ - 1 - i, max_distance, min_cost, gap_cost);
    }

    sets_.resize(static_cast<size_t>(kmers_));
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
      uint64_t id = 0;
      if(prefilter_detail::segment_key(padded_seq, offset, len, &id)) {
        sets_[static_cast<size_t>(i)].insert(id);
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
      const prefilter_detail::KeySet& set = sets_[static_cast<size_t>(i)];
      for(int shift = -max_shift; shift <= max_shift; ++shift) {
        uint64_t id = 0;
        if(prefilter_detail::segment_key(padded_query, offset + shift, len, &id) &&
           set.find(id) != set.end()) {
          return true;
        }
      }
    }
    return false;
  }

 private:
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

  bool active_ = false;
  int padded_len_ = 0;
  int kmers_ = 0;
  std::vector<int> klen_;
  std::vector<int> max_shift_;
  std::vector<prefilter_detail::KeySet> sets_;
};

}  // namespace startree

#endif  // STARTREE_PREFILTER_H
