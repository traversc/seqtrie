#ifndef STARTREE_BITS_H
#define STARTREE_BITS_H

#include <limits>

#if defined(_MSC_VER)
#  include <intrin.h>
#endif

namespace startree {

inline int ctz_u(unsigned x) noexcept {
  if(x == 0u) {
    return std::numeric_limits<unsigned>::digits;
  }
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_ctz(x);
#elif defined(_MSC_VER)
  unsigned long i;
  _BitScanForward(&i, x);
  return static_cast<int>(i);
#else
  int n = 0;
  while((x & 1u) == 0u) {
    x >>= 1;
    ++n;
  }
  return n;
#endif
}

}  // namespace startree

#endif  // STARTREE_BITS_H
