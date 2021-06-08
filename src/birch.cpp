#include <iostream>

#include "birch.h"
#include "Integer.h"

int main()
{
  Integer<Z> a = 240;
  Integer<Z> b = 46;

  EuvlideanDomain< Integer<Z> >::XGcdRes d_s_t = a.xgcd(b);
  
  std::cout << d_s_t << std::endl;
  
  return 0;
}
