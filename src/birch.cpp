#include <iostream>

#include "birch.h"
#include "testInteger.h"

int main()
{
  testInteger<Z64> test1();
  testInteger<Z> test2();
  testInteger<Z128> test3();
  
  return 0;
}
