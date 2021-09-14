#ifndef TREEDIST_UTILITY_H
#define TREEDIST_UTILITY_H

#include <type_traits>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <string>

#ifdef span_CONFIG_CONTRACT_VIOLATION_TERMINATES
#undef span_CONFIG_CONTRACT_VIOLATION_TERMINATES
#endif
// #define span_CONFIG_CONTRACT_VIOLATION_THROWS 1 // turn this on for testing
#define span_CONFIG_CONTRACT_LEVEL_OFF 1 // turn this on for production
#include "nonstd/span.hpp"

#include "simple_array/small_array.h"
#include "simple_array/nullable_array.h"
#include "simple_array/simple_array.h"


// https://stackoverflow.com/questions/16260033/reinterpret-cast-between-char-and-stduint8-t-safe
// reinterpret casts between char, unsigned char and uint8_t are allowed provided they are all 1 byte
// this should basically be always true for R platforms
static_assert(std::is_same<std::uint8_t, char>::value ||
              std::is_same<std::uint8_t, unsigned char>::value,
              "std::uint8_t should be char or unsigned char.");


namespace treedist {

// constexpr test for std::unique_ptr
template <class T> struct is_std_unique_ptr : std::false_type {};
template <class T> struct is_std_unique_ptr<std::unique_ptr<T>> : std::true_type {};

// constexpr test for std::array
template<class T> struct is_std_array : std::false_type {};
template<class T, std::size_t N> struct is_std_array<std::array<T,N>> : std::true_type {};

// typedef the vector used for radix branch; 
// seperate typedefs for easy switching and benchmarking
// alternative, boost::small_vector<T,44> but should be further benched; folly::small_vector has too much overhead
template <typename T> using rvector = std::vector<T>;
template <typename K, typename V> using runordered_map = std::unordered_map<K, V>;

// typedefs to use generically
template <typename T> using std_vector = std::vector<T>;
template <typename K, typename V> using std_unordered_map = std::unordered_map<K, V>;

using uspan = nonstd::span<const uint8_t>;
using cspan = nonstd::span<const char>;

template <typename T> void print_branch(const T & branch) {
  for(size_t i=0; i<branch.size(); ++i) {
    std::cout << static_cast<int>(branch[i]) << " ";
  }
  std::cout << std::endl;
}

inline uspan subspan(const uspan x, const size_t start, const size_t len = -1) {
  size_t rlen = std::min(len, x.size() - start);
  return uspan(x.data() + start, rlen);
}

template <typename T> inline T subvector(const uspan x, const size_t start, const size_t len = -1) {
  size_t rlen = std::min(len, x.size() - start);
  T result(rlen);
  std::copy(x.data() + start, x.data() + start + rlen, result.data());
  // std::memcpy(result.data(), x.data() + start, rlen);
  return result;
}

template <> inline std::basic_string<uint8_t> subvector<std::basic_string<uint8_t>>(const uspan x, const size_t start, const size_t len) {
  size_t rlen = std::min(len, x.size() - start);
  std::basic_string<uint8_t> result(rlen, 0); // if T is std::basic_string, there's no constructor that takes just length so we need this specialization
  std::copy(x.data() + start, x.data() + start + rlen, &result[0]); // data returns const pointer for basic_string
  // std::memcpy(result.data(), x.data() + start, rlen);
  return result;
}

template <typename T> inline void appendspan(T & x, const uspan y) {
  size_t xs = x.size();
  x.resize(xs + y.size());
  std::copy(y.data(), y.data() + y.size(), x.data() + xs);
  // std::memcpy(x.data() + xs, y.data(), y.size());
}

template <> inline void appendspan<std::basic_string<uint8_t>>(std::basic_string<uint8_t> & x, const uspan y) {
  size_t xs = x.size();
  x.resize(xs + y.size());
  std::copy(y.data(), y.data() + y.size(), &x[0] + xs);
  // std::memcpy(x.data() + xs, y.data(), y.size());
}


}


#endif // include guard