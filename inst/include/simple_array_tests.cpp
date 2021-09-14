#include <iostream>
#include "simple_array/small_nullable_array.h"
#include "simple_array/simple_array.h"
#include "simple_array/small_array.h"
#include "simple_array/nullable_array.h"

void small_array_test1() {
  std::cout << "small_array test1" << std::endl;
  trqwe::small_array<int,std::allocator<int>, uint32_t,5> x(2);
  trqwe::small_array<int,std::allocator<int>, uint32_t,5> x2(4);
  std::cout << "sizeof should be 32: " << sizeof(x) << std::endl;
  std::cout << "sizeof should be 32: " << sizeof(x2) << std::endl;
  std::cout << "is_stack should be true: " << x2.is_stack() << std::endl;
  x2 = x;
  trqwe::small_array<int,std::allocator<int>, uint32_t,5> x3 = std::move(x2);
  std::cout << "is_stack should be true: " << x3.is_stack() << std::endl;
  x3.resize(16);
  std::cout << "is_stack should be false: " << x3.is_stack() << std::endl;
  x3 = trqwe::small_array<int, std::allocator<int>, uint32_t,5>(20);
  std::cout << "is_stack should be false: " << x3.is_stack() << std::endl;
  x3[0] = -1;
  x3[x3.size() - 1] = 5;
  std::cout << "x3[0] should be -1: " << (int)x3[0] << std::endl;
  std::cout << "x3[size - 1] should be 5: " << (int)x3[x3.size() - 1] << std::endl;
  x3[x3.size() - 1] = 254;
  std::cout << "x3[size - 1] should be 254: " << (int)x3[x3.size() - 1] << std::endl;
}

void nullable_array_test1() {
  std::cout << "nullable_array test1" << std::endl;
  trqwe::nullable_array<int,std::allocator<int>, uint32_t> x(2);
  trqwe::nullable_array<int,std::allocator<int>, uint32_t> x2(4);
  std::cout << "sizeof should be 16: " << sizeof(x) << std::endl;
  std::cout << "sizeof should be 16: " << sizeof(x2) << std::endl;
  x2.nullify();
  std::cout << "is_null should be true: " << x2.is_null() << std::endl;
  x2 = x;
  trqwe::nullable_array<int,std::allocator<int>, uint32_t> x3 = std::move(x2);
  std::cout << "is_null should be true: " << x3.is_null() << std::endl;
  x3.resize(16);
  std::cout << "is_null should be false: " << x3.is_null() << std::endl;
  x3 = trqwe::nullable_array<int, std::allocator<int>, uint32_t>(20);
  std::cout << "is_null should be false: " << x3.is_null() << std::endl;
  x3[0] = -1;
  x3[x3.size() - 1] = 5;
  std::cout << "x3[0] should be -1: " << (int)x3[0] << std::endl;
  std::cout << "x3[size - 1] should be 5: " << (int)x3[x3.size() - 1] << std::endl;
  x3[x3.size() - 1] = 254;
  std::cout << "x3[size - 1] should be 254: " << (int)x3[x3.size() - 1] << std::endl;
}

void simple_array_test1() {
  std::cout << "simple_array test1" << std::endl;
  trqwe::simple_array<int,std::allocator<int>, uint32_t> x(2);
  trqwe::simple_array<int,std::allocator<int>, uint32_t> x2(4);
  std::cout << "sizeof should be 16: " << sizeof(x) << std::endl;
  std::cout << "sizeof should be 16: " << sizeof(x2) << std::endl;
  x2 = x;
  trqwe::simple_array<int,std::allocator<int>, uint32_t> x3 = std::move(x2);
  x3.resize(16);
  x3 = trqwe::simple_array<int, std::allocator<int>, uint32_t>(20);
  x3[0] = -1;
  x3[x3.size() - 1] = 5;
  std::cout << "x3[0] should be -1: " << (int)x3[0] << std::endl;
  std::cout << "x3[size - 1] should be 5: " << (int)x3[x3.size() - 1] << std::endl;
  x3[x3.size() - 1] = 254;
  std::cout << "x3[size - 1] should be 254: " << (int)x3[x3.size() - 1] << std::endl;
}

void small_nullable_array_test1() {
  std::cout << "small_nullable_array test1" << std::endl;
  trqwe::small_nullable_array<int,std::allocator<int>, uint32_t,5> x(2);
  trqwe::small_nullable_array<int,std::allocator<int>, uint32_t,5> x2(4);
  std::cout << "sizeof should be 32: " << sizeof(x) << std::endl;
  std::cout << "sizeof should be 32: " << sizeof(x2) << std::endl;
  std::cout << "is_stack should be true: " << x2.is_stack() << std::endl;
  x2.nullify();
  std::cout << "is_null should be true: " << x2.is_null() << std::endl;
  x2 = x;
  trqwe::small_nullable_array<int,std::allocator<int>, uint32_t,5> x3 = std::move(x2);
  std::cout << "is_stack should be true: " << x3.is_stack() << std::endl;
  std::cout << "is_null should be true: " << x3.is_null() << std::endl;
  x3.resize(16);
  std::cout << "is_stack should be false: " << x3.is_stack() << std::endl;
  std::cout << "is_null should be false: " << x3.is_null() << std::endl;
  x3 = trqwe::small_nullable_array<int, std::allocator<int>, uint32_t,5>(20);
  std::cout << "is_stack should be false: " << x3.is_stack() << std::endl;
  std::cout << "is_null should be false: " << x3.is_null() << std::endl;
  x3[0] = -1;
  x3[x3.size() - 1] = 5;
  std::cout << "x3[0] should be -1: " << (int)x3[0] << std::endl;
  std::cout << "x3[size - 1] should be 5: " << (int)x3[x3.size() - 1] << std::endl;
  x3[x3.size() - 1] = 254;
  std::cout << "x3[size - 1] should be 254: " << (int)x3[x3.size() - 1] << std::endl;
}

int main() {
  small_array_test1();
  nullable_array_test1();
  simple_array_test1();
  small_nullable_array_test1();
  return 0;
}
