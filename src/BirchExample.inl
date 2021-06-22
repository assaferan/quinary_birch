#include <vector>

#include "birch.h"
#include "Integer.h"
#include "QuadFormInt.h"

template<typename R, size_t n>
BirchExample<R,n>::BirchExample(const typename QuadFormZZ<R,n>::SymVec & q,
				const R & spinor,
				size_t d,
				const std::vector< std::vector< std::vector<R> > > & aps)
{
  for (size_t i = 0; i < n*(n+1)/2; i++)
    coeffs[i] = q[i];
  spinor_prime = spinor;
  dim = d;

  std::shared_ptr<const RationalField<Z> > QQ = std::make_shared<const RationalField<Z> >();
  UnivariatePolyRat<Z> f = UnivariatePolyRat<Z>::x(QQ) - Rational<Z>::one();
  std::shared_ptr< const NumberField<Z> > QNF = std::make_shared< const NumberField<Z> >(f);

  evs.resize(aps.size());
  for (size_t k = 0; k < aps.size(); k++) {
    Integer<R> p = R(2);
    for (size_t j = 0; j < aps[k][0].size(); j++) {
      std::vector< NumberFieldElement<Z> > vec;
      for (size_t i = 0; i < aps[k].size(); i++) {
	Rational<Z> ap_rat = birch_util::convertInteger<R,Z>(aps[k][i][j]);
	NumberFieldElement<Z> ap_nf(QNF, ap_rat);
	vec.push_back(ap_nf);
      }
      evs[k][p.num()] = vec;
      p = p.nextPrime();
    }
  }
}

template<typename R, size_t n>
inline BirchExample<Z64,3> BirchExample<R,n>::getExample_GV_7_2(void)
{
  std::vector< std::vector< std::vector<Z64> > > aps(1);
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

  aps[0].push_back(eis);
  aps[0].push_back(cusp);

  QuadFormZZ<Z64,3>::SymVec coeffs = {2,0,2,1,0,6};
  
  BirchExample<Z64,3> example(coeffs, 1, 2, aps);

  return example;
}

template<typename R, size_t n>
inline BirchExample<Z64,4> BirchExample<R,n>::getExample_GV_7_3(void)
{
  // aps[k][j] are the T_p^(k+1) eigenvalues of the form f_j
  std::vector< std::vector< std::vector<Z64> > > > aps(2);
  
  std::vector<Z64> eis,cusp,eis2,a,b;

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

  eis2.resize(cusp.size());
  a.resize(cusp.size());
  b.resize(cusp.size());
  for (size_t i = 0; i < cusp.size(); i++)
    eis2[i] = eis[i]*eis[i];
  for (size_t i = 0; i < cusp.size(); i++)
    a[i] = cusp[i]*cusp[i];
  for (size_t i = 0; i < cusp.size(); i++)
    b[i] = eis[i]*cusp[i];
  
  aps[0].push_back(eis2);
  aps[0].push_back(a);
  aps[0].push_back(b);

  std::vector<Z64> eis2_2,a_2,b_2;
  eis2_2.resize(cusp.size());
  a_2.resize(cusp.size());
  b_2.resize(cusp.size());
  for (size_t i = 0; i < cusp.size(); i++)
    eis2_2[i] = 2*eis[i]*(eis[i]-1);
  for (size_t i = 0; i < cusp.size(); i++)
    a_2[i] = 2*(cusp[i]*cusp[i] - eis[i]);
  for (size_t i = 0; i < cusp.size(); i++)
    b_2[i] = cusp[i]*cusp[i] + eis[i]*(eis[i]-2);

  aps[1].push_back(eis2);
  aps[1].push_back(a);
  aps[1].push_back(b);
  
  QuadFormZZ<Z64,4>::SymVec coeffs = {2,0,2,0,1,6,1,0,0,6};
  
  BirchExample<Z64,4> example(coeffs, 1, 3, aps);

  return example;
}
