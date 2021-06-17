#ifndef __TEST_BIRCH_H_
#define __TEST_BIRCH_H_

#include <map>
#include <vector>

#include "birch.h"
#include "BirchExample.h"
#include "Genus.h"
#include "QuadFormInt.h"

template<typename R, size_t n>
class TestBirch
{
public:
  TestBirch(const typename QuadFormZZ<R,n>::SymVec & coeffs);

  TestBirch(const BirchExample<R,n> &);

  void testDim(const R & spinor_prime, size_t dim) const;

  void testEigenvalues(const R & spinor_prime,
		       const std::map< R, std::vector< NumberFieldElement<Z> > > & evs) const;
  
protected:
  std::shared_ptr< Genus<R,n> > _p_genus;

  inline void _init(const typename QuadFormZZ<R,n>::SymVec & coeffs);
  
};

void runBirchTests(void);

#include "TestBirch.inl"

#endif // __TEST_BIRCH_H_
