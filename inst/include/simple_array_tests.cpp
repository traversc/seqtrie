#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <type_traits>
#include <utility>

#include "simple_array/nullable_array.h"
#include "simple_array/simple_array.h"
#include "simple_array/small_array.h"
#include "simple_array/small_nullable_array.h"

namespace {

template<uint32_t N>
using stack_size = std::integral_constant<uint32_t, N>;

void expect(bool condition, const char * const message) {
  if(!condition) {
    std::cerr << "FAILED: " << message << std::endl;
    std::exit(1);
  }
}

template<class Array>
void fill_linear(Array & x, const int offset = 0) {
  for(typename Array::size_type i = 0; i < x.size(); ++i) {
    x[i] = offset + static_cast<int>(i);
  }
}

template<class Array>
void expect_linear(Array const & x, const int offset, const char * const message) {
  for(typename Array::size_type i = 0; i < x.size(); ++i) {
    if(x[i] != offset + static_cast<int>(i)) {
      std::cerr << "FAILED: " << message << " at index " << i << std::endl;
      std::exit(1);
    }
  }
}

template<class Array>
void expect_prefix(Array const & x, typename Array::size_type prefix_size, const int offset,
                   const char * const message) {
  for(typename Array::size_type i = 0; i < prefix_size; ++i) {
    if(x[i] != offset + static_cast<int>(i)) {
      std::cerr << "FAILED: " << message << " at index " << i << std::endl;
      std::exit(1);
    }
  }
}

template<class Array>
void check_const_api() {
  static_assert(
    std::is_same<
      decltype(std::declval<Array const &>().data()),
      typename Array::const_pointer_type
    >::value,
    "const data() should return const_pointer_type"
  );
  static_assert(
    std::is_same<
      decltype(std::declval<Array const &>().begin()),
      typename Array::const_pointer_type
    >::value,
    "const begin() should return const_pointer_type"
  );
  static_assert(
    std::is_same<
      decltype(std::declval<Array const &>().end()),
      typename Array::const_pointer_type
    >::value,
    "const end() should return const_pointer_type"
  );
  static_assert(
    std::is_same<
      decltype(std::declval<Array const &>()[typename Array::size_type()]),
      typename Array::const_reference_type
    >::value,
    "const operator[] should return const_reference_type"
  );
}

void small_array_test() {
  using array_t = trqwe::small_array<int, std::allocator<int>, uint32_t, stack_size<5>>;
  check_const_api<array_t>();

  array_t stack_values(4);
  fill_linear(stack_values, 10);
  expect(stack_values.is_stack(), "small_array should use stack storage for small sizes");
  expect_linear(stack_values, 10, "small_array stack contents");

  array_t copied(stack_values);
  expect(copied.size() == 4, "small_array copy should preserve size");
  expect(copied.is_stack(), "small_array copy should preserve stack storage");
  expect_linear(copied, 10, "small_array copy contents");

  array_t moved(std::move(stack_values));
  expect(moved.size() == 4, "small_array move should preserve size");
  expect(moved.is_stack(), "small_array move should preserve stack storage");
  expect_linear(moved, 10, "small_array move contents");

  moved.resize(8);
  expect(!moved.is_stack(), "small_array should move to heap when resized past stack capacity");
  expect(moved.size() == 8, "small_array resize should update size");
  expect_prefix(moved, 4, 10, "small_array resize should preserve prefix");

  moved.resize(3);
  expect(moved.is_stack(), "small_array should return to stack storage when resized small");
  expect(moved.size() == 3, "small_array shrink should update size");
  expect_linear(moved, 10, "small_array shrink should preserve prefix");

  moved.reset(12);
  expect(moved.size() == 12, "small_array reset should update size");
  expect(!moved.is_stack(), "small_array reset should allocate heap storage when needed");

  array_t filled(3, 7);
  for(uint32_t i = 0; i < filled.size(); ++i) {
    expect(filled[i] == 7, "small_array value constructor should fill elements");
  }
}

void simple_array_test() {
  using array_t = trqwe::simple_array<int, std::allocator<int>, uint32_t>;
  check_const_api<array_t>();

  array_t empty;
  expect(empty.size() == 0, "simple_array default constructor should create empty array");
  expect(empty.data() == nullptr, "simple_array empty array should use nullptr data");

  array_t values(4);
  fill_linear(values, 20);
  expect_linear(values, 20, "simple_array contents");

  array_t copied = values;
  expect(copied.size() == 4, "simple_array copy should preserve size");
  expect_linear(copied, 20, "simple_array copy contents");

  array_t moved(std::move(values));
  expect(moved.size() == 4, "simple_array move should preserve size");
  expect_linear(moved, 20, "simple_array move contents");

  moved.resize(7);
  expect(moved.size() == 7, "simple_array resize should update size");
  expect_prefix(moved, 4, 20, "simple_array resize should preserve prefix");

  moved.reset(0);
  expect(moved.size() == 0, "simple_array reset should allow zero size");
  expect(moved.data() == nullptr, "simple_array zero-size reset should clear data pointer");

  array_t filled(3, -2);
  for(uint32_t i = 0; i < filled.size(); ++i) {
    expect(filled[i] == -2, "simple_array value constructor should fill elements");
  }
}

void nullable_array_test() {
  using array_t = trqwe::nullable_array<int, std::allocator<int>, uint32_t>;
  check_const_api<array_t>();

  array_t x;
  expect(x.is_null(), "nullable_array default constructor should create null array");
  expect(x.data() == nullptr, "nullable_array null array should expose nullptr data");
  expect(x.begin() == nullptr, "nullable_array null array should expose nullptr begin");
  expect(x.end() == nullptr, "nullable_array null array should expose nullptr end");

  x.resize(4);
  expect(!x.is_null(), "nullable_array resize should materialize storage");
  fill_linear(x, 30);
  expect_linear(x, 30, "nullable_array resize contents");

  array_t copied = x;
  expect(!copied.is_null(), "nullable_array copy should preserve non-null state");
  expect_linear(copied, 30, "nullable_array copy contents");

  array_t moved(std::move(x));
  expect(!moved.is_null(), "nullable_array move should preserve non-null state");
  expect_linear(moved, 30, "nullable_array move contents");

  moved.resize(2);
  expect(moved.size() == 2, "nullable_array shrink should update size");
  expect_linear(moved, 30, "nullable_array shrink should preserve prefix");

  moved.nullify();
  expect(moved.is_null(), "nullable_array nullify should restore null state");
  expect(moved.data() == nullptr, "nullable_array nullify should clear data pointer");

  moved.reset(5);
  expect(!moved.is_null(), "nullable_array reset should materialize storage");
  expect(moved.size() == 5, "nullable_array reset should update size");

  moved.resize(array_t::nullsize);
  expect(moved.is_null(), "nullable_array resize(nullsize) should nullify");
  expect(moved.data() == nullptr, "nullable_array resize(nullsize) should clear data pointer");
}

void small_nullable_array_test() {
  using array_t = trqwe::small_nullable_array<int, std::allocator<int>, uint32_t, stack_size<5>>;
  check_const_api<array_t>();

  array_t x;
  expect(!x.is_null(), "small_nullable_array default constructor should not be null");
  expect(x.is_stack(), "small_nullable_array default constructor should use stack storage");
  expect(x.size() == 0, "small_nullable_array default constructor should create empty array");

  x.resize(4);
  fill_linear(x, 40);
  expect(!x.is_null(), "small_nullable_array resize should preserve non-null state");
  expect(x.is_stack(), "small_nullable_array small resize should use stack storage");
  expect_linear(x, 40, "small_nullable_array stack contents");

  array_t copied = x;
  expect(!copied.is_null(), "small_nullable_array copy should preserve non-null state");
  expect(copied.is_stack(), "small_nullable_array copy should preserve stack storage");
  expect_linear(copied, 40, "small_nullable_array copy contents");

  array_t moved(std::move(x));
  expect(!moved.is_null(), "small_nullable_array move should preserve non-null state");
  expect(moved.is_stack(), "small_nullable_array move should preserve stack storage");
  expect_linear(moved, 40, "small_nullable_array move contents");

  moved.resize(9);
  expect(!moved.is_null(), "small_nullable_array heap resize should stay non-null");
  expect(!moved.is_stack(), "small_nullable_array should move to heap when resized past stack capacity");
  expect_prefix(moved, 4, 40, "small_nullable_array heap resize should preserve prefix");

  moved.resize(3);
  expect(moved.is_stack(), "small_nullable_array should return to stack storage when resized small");
  expect_linear(moved, 40, "small_nullable_array shrink should preserve prefix");

  moved.nullify();
  expect(moved.is_null(), "small_nullable_array nullify should set null state");
  expect(moved.is_stack(), "small_nullable_array null state should still use stack backing");
  expect(moved.data() == nullptr, "small_nullable_array null state should expose nullptr data");

  array_t moved_null(std::move(moved));
  expect(moved_null.is_null(), "small_nullable_array move should preserve null state");
  expect(moved_null.data() == nullptr, "small_nullable_array moved null state should expose nullptr data");

  moved_null.reset(2);
  expect(!moved_null.is_null(), "small_nullable_array reset should materialize storage");
  expect(moved_null.is_stack(), "small_nullable_array small reset should use stack storage");
  moved_null.resize(array_t::nullsize);
  expect(moved_null.is_null(), "small_nullable_array resize(nullsize) should nullify");
  expect(moved_null.data() == nullptr, "small_nullable_array resize(nullsize) should clear data pointer");
}

} // end anonymous namespace

int main() {
  small_array_test();
  simple_array_test();
  nullable_array_test();
  small_nullable_array_test();
  std::cout << "simple_array tests passed" << std::endl;
  return 0;
}
