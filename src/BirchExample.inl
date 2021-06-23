#include <vector>

#include "birch.h"
#include "Integer.h"
#include "QuadFormInt.h"

// In this constructor we assume the aps are the actual eigenvalues and use them to populate both evs and traces
template<typename R, size_t n>
BirchExample<R,n>::BirchExample(const QuadFormZZ<R,n> & q,
				const R & spinor,
				size_t d,
				const std::vector< std::vector< std::vector<R> > > & aps)
{
  qf = q;
  spinor_prime = spinor;
  dim = d;

  std::shared_ptr<const RationalField<Z> > QQ = std::make_shared<const RationalField<Z> >();
  UnivariatePolyRat<Z> f = UnivariatePolyRat<Z>::x(QQ) - Rational<Z>::one();
  std::shared_ptr< const NumberField<Z> > QNF = std::make_shared< const NumberField<Z> >(f);

  evs.resize(aps.size());
  traces.resize(aps.size());
  for (size_t k = 0; k < aps.size(); k++) {
    Integer<R> p = R(2);
    for (size_t j = 0; j < aps[k][0].size(); j++) {
      std::vector< NumberFieldElement<Z> > vec_nf;
      std::vector<R> vec;
      for (size_t i = 0; i < aps[k].size(); i++) {
	Rational<Z> ap_rat = birch_util::convertInteger<R,Z>(aps[k][i][j]);
	NumberFieldElement<Z> ap_nf(QNF, ap_rat);
	vec_nf.push_back(ap_nf);
	vec.push_back(aps[k][i][j]);
      }
      evs[k][p.num()] = vec_nf;
      traces[k][p.num()] = vec;
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
  QuadFormZZ<Z64,3> qf(coeffs);
    
  BirchExample<Z64,3> example(qf, 1, 2, aps);

  return example;
}

template<typename R, size_t n>
inline BirchExample<Z64,3> BirchExample<R,n>::getExample_CMF_49_2_a_a(void)
{
  std::vector< std::vector< std::vector<Z64> > > aps(1);
  std::vector<Z64> eis;
  std::vector<Z64> cusp;

  // initialize the eisenstein form
  Integer<Z64> p = 2;
  while (p.num() < 100) {
    if (p.num() == 7)
      eis.push_back(0);
    else
      eis.push_back(p.num()+1);
    p = p.nextPrime();
  }

  // initialize the cusp form
  cusp = { 1, 0, 0, 0, 4, 0, 0, 0, 8, 2, 0, -6, 0, -12, 0, -10, 0, 0, 4, 16, 0, 8, 0, 0, 0};

  aps[0].push_back(eis);
  aps[0].push_back(cusp);

  QuadFormZZ<Z64,3>::SymVec coeffs = {6,1,6,1,-1,20};
  QuadFormZZ<Z64,3> qf(coeffs);
  
  BirchExample<Z64,3> example(qf, 1, 2, aps);

  return example;
}

template<typename R, size_t n>
inline BirchExample<Z64,4> BirchExample<R,n>::getExample_GV_7_3(void)
{
  // aps[k][j] are the T_p^(k+1) eigenvalues of the form f_j
  std::vector< std::vector< std::vector<Z64> > > aps(2);
  
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

  aps[1].push_back(eis2_2);
  aps[1].push_back(a_2);
  aps[1].push_back(b_2);
  
  QuadFormZZ<Z64,4>::SymVec coeffs = {2,0,2,0,1,6,1,0,0,6};
  QuadFormZZ<Z64,4> qf(coeffs);
  
  BirchExample<Z64,4> example(qf, 1, 3, aps);

  return example;
}

template<typename R, size_t n>
inline BirchExample<Z64,5> BirchExample<R,n>::getExample_RT_Table1(void)
{
  std::vector< std::vector< std::vector<Z64> > > aps(2);

  // initialize the cusp form
  std::vector<Z64> traces1 = {-8,-10,-4,-14,-22,-4,-47,-12,41,50,-504,-102,174,30,42,156,-252,472,106,-481,-744,927,-632,-297,2,-992,-1222,1436,-954,19,516,-258, 1080, 1030, -974, -1119, 1152, 108, -2707, -182, 2568, -2804, -3035, 583, 2276, 6754, 360, 3569, -3346, 2220, -2780, -3878, -819, 6112, -5343, -808, 3592, 2954, -8334, -2942, 6360, -856, 3548, -6322, -9443, 108, 1596, -2129, 1856, 480, 1704, 4601, 6298, -4998, 7706, -18293, 5316, 4324, -4679, -3476, -910, 3552, -4878, 15213, -6909, -7130, 12908, -4005, -7334, -77, 12248, 6447, -14197, 1960, 3288};

  std::vector<Z64> traces2 = {10,11,-44,-9,-67,-158,260,41,-198,-187,2744,-730,800,442,-5052};

  aps[0].push_back(traces1);
  aps[1].push_back(traces2);

  QuadFormZZ<Z64,5> qf = QuadFormZZ<Z64,5>::getQuinaryForms(167);
  
  BirchExample<Z64,5> example(qf, 167, 1, aps);

  return example;
}
