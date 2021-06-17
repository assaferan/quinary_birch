#include <vector>

#include "birch.h"
#include "Integer.h"
#include "QuadFormInt.h"

template<typename R, size_t n>
BirchExample<R,n>::BirchExample(const typename QuadFormZZ<R,n>::SymVec & q,
				const R & spinor,
				size_t d,
				const std::vector< std::vector<R> > & aps)
{
  for (size_t i = 0; i < n*(n+1)/2; i++)
    coeffs[i] = q[i];
  spinor_prime = spinor;
  dim = d;

  std::shared_ptr<const RationalField<Z> > QQ = std::make_shared<const RationalField<Z> >();
  UnivariatePolyRat<Z> f = UnivariatePolyRat<Z>::x(QQ) - Rational<Z>::one();
  std::shared_ptr< const NumberField<Z> > QNF = std::make_shared< const NumberField<Z> >(f);
  
  Integer<R> p = R(2);
  for (size_t i = 0; i < aps.size(); i++) {
    std::vector< NumberFieldElement<Z> > vec;
    for (size_t j = 0; j < aps[i].size(); j++) {
      Rational<Z> ap_rat = birch_util::convertInteger<R,Z>(aps[i][j]);
      NumberFieldElement<Z> ap_nf(QNF, ap_rat);
      vec.push_back(ap_nf);
    }
    evs[p.num()] = vec;
    p = p.nextPrime();
  }
}

template<typename R, size_t n>
inline BirchExample<Z64,3> BirchExample<R,n>::getExample_GV_7_2(void)
{
  std::vector< std::vector<Z64> > aps;
  std::vector<Z64> eis;
  std::vector<Z64> cusp;

  // initialize the eisenstein form
  Integer<Z64> p = 2;
  while (p.num() < 100) {
    if (p.num() == 11)
      eis.push_back(0);
    else
      eis.push_back(p.num()+1);
    p = p.nextPrime();
  }

  // initialize the cusp form
  cusp = {-2, -1, 1, -2, 0, 4, -2, 0, -1, 0, 7, 3, -8, -6, 8, -6, 5, 12, -7, -3, 4, -10, -6, 15, -7};

  aps.push_back(eis);
  aps.push_back(cusp);

  QuadFormZZ<Z64,3>::SymVec coeffs = {2,0,2,1,0,6};
  
  BirchExample<Z64,3> example(coeffs, 1, 2, aps);

  return example;
}
