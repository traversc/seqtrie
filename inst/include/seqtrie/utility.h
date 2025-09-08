#ifndef seqtrie_UTILITY_H
#define seqtrie_UTILITY_H

#include <type_traits>
#include <vector>
#include <map> // or unordered_map
#include <algorithm>
#include <string>
#include <array>
#include <tuple>
#include <memory>
#include <limits.h> // INT_MAX
#include <iterator>
#include <utility>
#include <cstdint>
#include <cstddef>

#include "ankerl/unordered_dense.h"

#ifdef span_CONFIG_CONTRACT_VIOLATION_TERMINATES
#undef span_CONFIG_CONTRACT_VIOLATION_TERMINATES
#endif
// #define span_CONFIG_CONTRACT_VIOLATION_THROWS 1 // turn this on for testing
#define span_CONFIG_CONTRACT_LEVEL_OFF 1 // turn this on for production
#define span_FEATURE_NON_MEMBER_FIRST_LAST_SUB 1
#define span_FEATURE_MAKE_SPAN 1
#include "nonstd/span.hpp"

// requires boost
// #include <boost/mpl/string.hpp>
// #include <boost/mpl/for_each.hpp>
// #include <boost/mpl/range_c.hpp>


// https://stackoverflow.com/questions/16260033/reinterpret-cast-between-char-and-stduint8-t-safe
// reinterpret casts between char, unsigned char and uint8_t are allowed provided they are all 1 byte
// this should basically be always true for R platforms
static_assert(std::is_same<std::uint8_t, char>::value ||
              std::is_same<std::uint8_t, unsigned char>::value,
              "std::uint8_t should be char or unsigned char.");


namespace seqtrie {

// Common pair type for substitution costs
using pairchar_type = std::pair<char, char>;

// Unified cost map: substitution table + uniform gap costs
struct CostMap {
  ankerl::unordered_dense::map<pairchar_type, int> char_cost_map; // substitution costs
  int gap_cost;        // linear gap cost and affine extension cost
  int gap_open_cost;   // affine gap opening cost
};

// constexpr test for std::unique_ptr
template <class T> struct is_std_unique_ptr : std::false_type {};
template <class T> struct is_std_unique_ptr<std::unique_ptr<T>> : std::true_type {};

// constexpr test for std::array
template<class T> struct is_std_array : std::false_type {};
template<class T, std::size_t N> struct is_std_array<std::array<T,N>> : std::true_type {};

// print a branch for debug purposes
// template <typename T> void print_branch(const T & branch) {
//   for(size_t i=0; i<branch.size(); ++i) {
//     std::cout << static_cast<int>(branch[i]) << " ";
//   }
//   std::cout << std::endl;
// }

// subspan -- aAlready implemented in nonstd::span
// inline uspan subspan(const uspan x, const size_t start, const size_t len = -1) {
//   size_t rlen = std::min(len, x.size() - start);
//   return uspan(x.data() + start, rlen);
// }
// inline cspan subspan(const cspan x, const size_t start, const size_t len = -1) {
//   size_t rlen = std::min(len, x.size() - start);
//   return cspan(x.data() + start, rlen);
// }

// create a vector or array type, we use this approach because std::string doesn't have a constructor with just size
// Just in case we want to allow template return values (e.g. as string),
// right now return values are hard-coded as std::vector<char>
// also nonstd::span doesn't play nice with std::string, so it's probably not worth it
template <typename T> inline T array_allocate(const size_t size) { return T(size); }
template <> inline std::string array_allocate(const size_t size) { return std::string(size, 0); }
// clang 18 deprecated
// template <> inline std::basic_string<uint8_t> array_allocate(const size_t size) { return std::basic_string<uint8_t>(size, 0); }

template <typename T> inline typename T::value_type * array_data(T & x) { return x.data(); }
template <> inline char * array_data(std::string & x) { return &x[0]; }
// clang 18 deprecated
// template <> inline uint8_t * array_data(std::basic_string<uint8_t> & x) { return &x[0]; }
  
// subvector
template <typename T, typename F> inline T subvector(const F & x, const size_t start, const size_t len = -1) {
  size_t rlen = std::min(len, x.size() - start);
  T result(rlen);
  std::copy(x.data() + start, x.data() + start + rlen, result.data());
  return result;
}

// appendspan -- append span to vector
template <typename T, typename S> inline void appendspan(T & x, const S & y) {
  static_assert(std::is_same<typename T::value_type, typename S::value_type>::value, "appendspan x and y value_type must be the same");
  size_t xs = x.size();
  x.resize(xs + y.size());
  std::copy(y.data(), y.data() + y.size(), x.data() + xs);
}

template <typename T> inline T iota_range(const typename T::value_type value, const size_t len) {
  T result(len);
  std::iota(result.begin(), result.end(), value);
  return result;
}

// basic_string<uint8_t> specialization
// template <> inline void appendspan<std::basic_string<uint8_t>>(std::basic_string<uint8_t> & x, const uspan y) {
//   size_t xs = x.size();
//   x.resize(xs + y.size());
//   std::copy(y.data(), y.data() + y.size(), &x[0] + xs);
// }
// string specialization
// template <> inline void appendspan<std::string>(std::string & x, const uspan y) {
//   size_t xs = x.size();
//   x.resize(xs + y.size());
//   std::copy(y.data(), y.data() + y.size(), &x[0] + xs);
// }
}

inline std::string ptr_tostring(const void * ptr) {
  std::stringstream ss;
  ss << ptr;
  return ss.str();
}

// ==== Myers Bit-Vector (shared types) ====
namespace seqtrie {
  struct MyersPattern {
    using span_type = nonstd::span<const char>;
    size_t m = 0;
    size_t words = 0;
    uint64_t last_mask = 0;
    uint64_t top_bit = 0;
    std::vector<uint64_t> valid_mask; // size words; last masked
    std::vector<uint64_t> Peq[256];   // per byte, size words

    MyersPattern() = default;
    MyersPattern(span_type pattern, const CostMap & cost_map) { init(pattern, cost_map, false); }
    MyersPattern(span_type pattern, const CostMap & cost_map, bool reverse_pairs) { init(pattern, cost_map, reverse_pairs); }

   private:
    inline void init(span_type pattern, const CostMap & cost_map, bool reverse_pairs) {
      m = pattern.size();
      words = (m + 63) / 64;
      if (m == 0) {
        last_mask = 0; top_bit = 0; valid_mask.clear();
        for (int c=0;c<256;++c) Peq[c].clear();
        return;
      }
      const size_t last_bits = ((m - 1) & 63) + 1;
      last_mask = (last_bits == 64) ? ~0ull : ((1ull << last_bits) - 1ull);
      top_bit   = 1ull << ((m - 1) & 63);
      valid_mask.assign(words, ~0ull);
      valid_mask[words - 1] = last_mask;
      for (int c = 0; c < 256; ++c) Peq[c].assign(words, 0ull);

      // Precompute positions for each byte in the pattern
      std::vector<size_t> pos[256];
      for (size_t i = 0; i < m; ++i) pos[static_cast<unsigned char>(pattern[i])].push_back(i);
      if (!reverse_pairs) {
        // For every zero-cost (a,b), set bits in Peq[b] where pattern[i]==a
        for (const auto & kv : cost_map.char_cost_map) {
          const char a = kv.first.first;
          const char b = kv.first.second;
          const int  v = kv.second;
          if (v != 0) continue;
          const auto & plist = pos[static_cast<unsigned char>(a)];
          if (plist.empty()) continue;
          auto & mask_vec = Peq[static_cast<unsigned char>(b)];
          for (size_t idx : plist) mask_vec[idx >> 6] |= (1ull << (idx & 63));
        }
      } else {
        // Reverse mapping: use zero-cost (a,b) but set bits in Peq[a] where pattern[i]==b
        for (const auto & kv : cost_map.char_cost_map) {
          const char a = kv.first.first;  // will be text char
          const char b = kv.first.second; // will be pattern char
          const int  v = kv.second;
          if (v != 0) continue;
          const auto & plist = pos[static_cast<unsigned char>(b)];
          if (plist.empty()) continue;
          auto & mask_vec = Peq[static_cast<unsigned char>(a)];
          for (size_t idx : plist) mask_vec[idx >> 6] |= (1ull << (idx & 63));
        }
      }
      // Mask last block
      for (int c = 0; c < 256; ++c) Peq[c][words - 1] &= last_mask;
    }
  };

  struct MyersState {
    std::vector<uint64_t> VP, VN;
    std::vector<uint64_t> D0, HP, HN;
    std::vector<uint64_t> HPs, HNs; // shifted buffers to avoid per-step allocation
    int score = 0;

    MyersState() = default;
    explicit MyersState(const MyersPattern & pat)
      : VP(pat.words, ~0ull), VN(pat.words, 0ull),
        D0(pat.words), HP(pat.words), HN(pat.words), HPs(pat.words), HNs(pat.words),
        score(static_cast<int>(pat.m)) {}

    inline void step(unsigned char tc, const MyersPattern & pat) {
      // D0 across blocks with carry: D0 = (((X & VP) + VP) ^ VP) | X
      uint64_t carry = 0ull;
      for (size_t w = 0; w < pat.words; ++w) {
        const uint64_t Eq = pat.Peq[tc][w];
        const uint64_t X = Eq | VN[w];
        const uint64_t a = X & VP[w];
        const uint64_t b = VP[w];
        uint64_t sum = a + b;
        const uint64_t carry1 = (sum < a);
        const uint64_t sum2 = sum + carry;
        const uint64_t carry2 = (sum2 < sum);
        carry = carry1 | carry2;
        D0[w] = (sum2 ^ VP[w]) | X;
      }
      // HP/HN and clamp
      for (size_t w = 0; w < pat.words; ++w) {
        HP[w] = VN[w] | ~(D0[w] | VP[w]);
        HN[w] = VP[w] & D0[w];
        HP[w] &= pat.valid_mask[w];
        HN[w] &= pat.valid_mask[w];
      }
      // Score update based on highest pattern bit
      if (HP[pat.words - 1] & pat.top_bit) ++score;
      else if (HN[pat.words - 1] & pat.top_bit) --score;
      // Shift left by 1 with carries
      uint64_t carry_hp = 0ull, carry_hn = 0ull;
      for (size_t w = 0; w < pat.words; ++w) {
        const uint64_t hp = HP[w];
        const uint64_t hn = HN[w];
        HPs[w] = (hp << 1) | carry_hp;
        HNs[w] = (hn << 1) | carry_hn;
        carry_hp = (hp >> 63) & 1ull;
        carry_hn = (hn >> 63) & 1ull;
      }
      // Add 1 to LSB of HP shift
      HPs[0] |= 1ull;
      // Update VP, VN using shifted values
      for (size_t w = 0; w < pat.words; ++w) {
        VP[w] = HNs[w] | ~(D0[w] | HPs[w]);
        VN[w] = HPs[w] & D0[w];
        VP[w] &= pat.valid_mask[w];
        VN[w] &= pat.valid_mask[w];
      }
    }
  };
}

#endif // include guard
