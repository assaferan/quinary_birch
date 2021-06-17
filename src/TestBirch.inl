#include <map>
#include <vector>

#include "birch.h"
#include "BirchExample.h"
#include "Eigenvector.h"
#include "Genus.h"
#include "Integer.h"
#include "QuadFormInt.h"

template<typename R, size_t n>
TestBirch<R,n>::TestBirch(const typename QuadFormZZ<R,n>::SymVec & coeffs)
{
  this->_init(coeffs);
}

template<typename R, size_t n>
inline void TestBirch<R,n>::testDim(const R & spinor_prime, size_t dim) const
{
  std::map<R,size_t> dims = this->_p_genus->dimensionMap();

  assert(dims[spinor_prime] == dim);

  return;
}

template<typename R, size_t n>
inline void TestBirch<R,n>::testEigenvalues(const R & spinor_prime,
					    const  std::map< R, std::vector< NumberFieldElement<Z> > > & evs) const
{
  EigenvectorManager<R,n> manager;
  std::vector< std::vector< NumberFieldElement<Z> > > evecs = _p_genus->eigenvectors()[spinor_prime];
  for (std::vector< NumberFieldElement<Z> > evec : evecs)
      manager.addEigenvector(_p_genus->eigenvector(evec, spinor_prime));
  manager.finalize();

  for (std::pair< R, std::vector< NumberFieldElement<Z> > > ev : evs) {
    Integer<R> p = ev.first;
    assert(ev.second == _p_genus->eigenvalues(manager, p.num()));
  }

  return;
}

template<typename R, size_t n>
inline TestBirch<R,n>::TestBirch(const BirchExample<R,n> & example)
{
  this->_init(example.coeffs);
  this->testDim(example.spinor_prime, example.dim);
  this->testEigenvalues(example.spinor_prime, example.evs);
}

template<typename R, size_t n>
inline void TestBirch<R,n>::_init(const typename QuadFormZZ<R,n>::SymVec & coeffs)
{
   QuadFormInt<R,n> q(coeffs);

  Integer<R> disc = q.discriminant();
  typename Integer<R>::FactorData facs = disc.factorization();

  std::vector<PrimeSymbol<R> > symbols;
  PrimeSymbol<R> symb;

  for (std::pair<Integer<R>, size_t> fa : facs) {
    symb.p = fa.first.num();
    symb.power = fa.second;
    symb.ramified = true;
    symbols.push_back(symb);
  }

  this->_p_genus = std::make_shared< Genus<R,n> >(q, symbols);
}

inline void runBirchTests(void) {
  TestBirch<Z64,3> test(BirchExample<Z64,3>::getExample_GV_7_2());
}
