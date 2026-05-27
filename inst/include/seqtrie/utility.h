#ifndef seqtrie_UTILITY_H
#define seqtrie_UTILITY_H

#include <type_traits>
#include <vector>
#include <algorithm>
#include <string>
#include <array>
#include <tuple>
#include <memory>
#include <limits>
#include <iterator>
#include <utility>
#include <cstdint>
#include <cstddef>

#include "simple_array/simple_array.h"

#ifdef span_CONFIG_CONTRACT_VIOLATION_TERMINATES
#undef span_CONFIG_CONTRACT_VIOLATION_TERMINATES
#endif
// #define span_CONFIG_CONTRACT_VIOLATION_THROWS 1 // turn this on for testing
#define span_CONFIG_CONTRACT_LEVEL_OFF 1 // turn this on for production
#define span_FEATURE_NON_MEMBER_FIRST_LAST_SUB 1
#define span_FEATURE_MAKE_SPAN 1
#include "nonstd/span.hpp"

#include <numeric>
#include <sstream>


// https://stackoverflow.com/questions/16260033/reinterpret-cast-between-char-and-stduint8-t-safe
// reinterpret casts between char, unsigned char and uint8_t are allowed provided they are all 1 byte
// this should basically be always true for R platforms
static_assert(std::is_same<std::uint8_t, char>::value ||
              std::is_same<std::uint8_t, unsigned char>::value,
              "std::uint8_t should be char or unsigned char.");


namespace seqtrie {

// Common pair type for substitution costs
using pairchar_type = std::pair<char, char>;

// Unified cost map: transposed dense substitution table + uniform gap costs.
// Stored as [target char][query char], matching radix traversal where the
// target/radix character is fixed while query position varies.
struct CostMap {
  static constexpr size_t alphabet_size = 256;
  std::array<int, alphabet_size * alphabet_size> substitution_cost{};
  std::array<unsigned char, alphabet_size> target_chars{};
  std::array<unsigned char, alphabet_size> target_char_present{};
  size_t target_char_count = 0;
  int gap_cost;        // linear gap cost and affine extension cost
  int gap_open_including_first_extension;   // affine open cost plus first extension

  static inline size_t byte_index(char x) {
    return static_cast<size_t>(static_cast<unsigned char>(x));
  }

  inline void add_target_char(char target_char) {
    const size_t idx = byte_index(target_char);
    if(target_char_present[idx] == 0) {
      target_char_present[idx] = 1;
      target_chars[target_char_count++] = static_cast<unsigned char>(target_char);
    }
  }

  inline void set_subst(char query_char, char target_char, int cost) {
    add_target_char(target_char);
    substitution_cost[byte_index(target_char) * alphabet_size + byte_index(query_char)] = cost;
  }

  inline const int * subst_for_target(char target_char) const {
    return substitution_cost.data() + byte_index(target_char) * alphabet_size;
  }

  inline int subst(char query_char, char target_char) const {
    return subst_for_target(target_char)[byte_index(query_char)];
  }
};

// Per-query substitution cache used by radix custom-cost search.
// Rows are stored by target/radix character and columns by query position:
//   cost[target_char][query_position]
// This means the DP update loop does not need the query string at all.
struct QueryCostCache {
  std::array<const int *, CostMap::alphabet_size> substitution_for_target{};
  trqwe::simple_array<int> substitution_cost;
  size_t query_len = 0;
  int gap_cost = 0;
  int gap_open_including_first_extension = 0;

  QueryCostCache() {
    substitution_for_target.fill(nullptr);
  }

  QueryCostCache(const CostMap & cost_map, nonstd::span<const char> query) {
    reset(cost_map, query);
  }

  inline void reset(const CostMap & cost_map, nonstd::span<const char> query) {
    substitution_for_target.fill(nullptr);
    query_len = query.size();
    gap_cost = cost_map.gap_cost;
    gap_open_including_first_extension = cost_map.gap_open_including_first_extension;

    const size_t n_targets = cost_map.target_char_count;
    substitution_cost.reset(n_targets * query_len);
    int * out = substitution_cost.data();

    for(size_t target_idx = 0; target_idx < n_targets; ++target_idx) {
      const char target_char = static_cast<char>(cost_map.target_chars[target_idx]);
      int * row = out + target_idx * query_len;
      substitution_for_target[CostMap::byte_index(target_char)] = row;
      for(size_t query_idx = 0; query_idx < query_len; ++query_idx) {
        row[query_idx] = cost_map.subst(query[query_idx], target_char);
      }
    }
  }

  inline const int * subst_for_target(char target_char) const {
    return substitution_for_target[CostMap::byte_index(target_char)];
  }
};

// constexpr test for std::unique_ptr
template <class T> struct is_std_unique_ptr : std::false_type {};
template <class T> struct is_std_unique_ptr<std::unique_ptr<T>> : std::true_type {};

// constexpr test for std::array
template<class T> struct is_std_array : std::false_type {};
template<class T, std::size_t N> struct is_std_array<std::array<T,N>> : std::true_type {};


// create a vector or array type, we use this approach because std::string doesn't have a constructor with just size
// Just in case we want to allow template return values (e.g. as string),
// right now return values are hard-coded as std::vector<char>
// also nonstd::span doesn't play nice with std::string, so it's probably not worth it
template <typename T> inline T array_allocate(const size_t size) { return T(size); }
template <> inline std::string array_allocate(const size_t size) { return std::string(size, 0); }

template <typename T> inline typename T::value_type * array_data(T & x) { return x.data(); }
template <> inline char * array_data(std::string & x) { return &x[0]; }
  
// subvector
template <typename T, typename F> inline T subvector(const F & x, const size_t start, const size_t len = std::numeric_limits<size_t>::max()) {
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

}

inline std::string ptr_tostring(const void * ptr) {
  std::stringstream ss;
  ss << ptr;
  return ss.str();
}


#endif // include guard
