#include <iostream>

#include "birch.h"
#include "birch_util.h"
#include "testInteger.h"

std::ostream & operator<<(std::ostream & os, const Z128 & z)
{
  os << birch_util::convert_Integer<Z128, Z>(z);
  return os;
}

int main()
{
  testInteger<Z64> test1;
  testInteger<Z> test2;
  testInteger<Z128> test3;
  
  return 0;
}
