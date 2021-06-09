#include <iostream>

#include "birch.h"
#include "birch_util.h"
#include "testInteger.h"
#include "testRational.h"
#include "FpElement.h"

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

  testRational<Z64> test4;
  testRational<Z> test5;
  testRational<Z128> test6;
  
  W64 seed = 1;
  std::shared_ptr< W16_Fp > GF = std::make_shared(3, seed);
  W16_FpElement a(GF);
  
  return 0;
}
