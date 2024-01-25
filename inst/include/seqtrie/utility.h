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
template <typename T, typename F> T subvector(const F & x, const size_t start, const size_t len = -1) {
  size_t rlen = std::min(len, x.size() - start);
  T result(rlen);
  std::copy(x.data() + start, x.data() + start + rlen, result.data());
  return result;
}

// appendspan -- append span to vector
template <typename T, typename S> void appendspan(T & x, const S & y) {
  static_assert(std::is_same<typename T::value_type, typename S::value_type>::value, "appendspan x and y value_type must be the same");
  size_t xs = x.size();
  x.resize(xs + y.size());
  std::copy(y.data(), y.data() + y.size(), x.data() + xs);
}

template <typename T> T iota_range(const typename T::value_type value, const size_t len) {
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

// hash for unordered_map with std::pair<char, char> as key
namespace std {
  template <>
  struct hash<std::pair<char, char>> {
    size_t operator()(const std::pair<char, char> & p) const {
      return ((p.first + 128) << 8) + (p.second + 128);
    }
  };
}

#endif // include guard
