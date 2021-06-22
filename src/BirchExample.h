#ifndef __BIRCH_EXAMPLE_H_
#define __BIRCH_EXAMPLE_H_

#include <map>
#include <vector>

#include "birch.h"
#include "NumberFieldElement.h"
#include "QuadFormInt.h"

template<typename R, size_t n>
class BirchExample
{
public:
  BirchExample(const typename QuadFormZZ<R,n>::SymVec &, const R &, size_t,
	       const std::vector< std::vector< std::vector<R> > > &);

  // n = 3
  static BirchExample<Z64,3> getExample_GV_7_2(void);
  static BirchExample<Z64,3> getExample_CMF_49_2_a_a(void);
  // n = 4
  static BirchExample<Z64,4> getExample_GV_7_3(void);
 

  typename QuadFormZZ<R,n>::SymVec coeffs;
  R spinor_prime;
  size_t dim;
  std::vector< std::map< R, std::vector< NumberFieldElement<Z> > > > evs;
  
};

#include "BirchExample.inl"

#endif // __BIRCH_EXAMPLE_H_
