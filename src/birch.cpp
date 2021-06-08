#include <iostream>

#include "birch.h"
#include "Integer.h"

int main()
{
  Integer<Z64> a = Z64(240);
  Integer<Z64> b = Z64(46);

  EuclideanDomain< Integer<Z64> >::XGcdRes d_s_t = a.xgcd(b);

  std::cout << std::get<0>(d_s_t) << ", ";
  std::cout << std::get<1>(d_s_t) << ", ";
  std::cout << std::get<2>(d_s_t);
  std::cout << std::endl;
  
  return 0;
}
