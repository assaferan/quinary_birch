#include <map>
#include <unordered_set>
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

  std::vector< EigenvalueVector > evalues(evecs.size());
  std::vector< EigenvalueVector > computed_evalues(evecs.size());
  
  for (std::pair< R, std::vector< NumberFieldElement<Z> > > ev : evs) {
    Integer<R> p = ev.first;
    std::vector< NumberFieldElement<Z> > computed = _p_genus->eigenvalues(manager, p.num());
    for (size_t i = 0; i < evecs.size(); i++) {
      evalues[i].vec.push_back(ev.second[i]);
      computed_evalues[i].vec.push_back(computed[i]);
    }
#ifdef DEBUG
    std::cerr << "ev.second = " << ev.second << std::endl;
    std::cerr << "computed eigenvalues: " << computed << std::endl;
#endif
    
    //    assert(ev.second == _p_genus->eigenvalues(manager, p.num()));
  }

  std::unordered_set< EigenvalueVector > evalue_set(evalues.begin(), evalues.end());
  std::unordered_set< EigenvalueVector > computed_set(computed_evalues.begin(), computed_evalues.end());

  assert( evalue_set == computed_set);

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
   QuadFormZZ<R,n> q(coeffs);

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
  TestBirch<Z64,3> test_7_2(BirchExample<Z64,3>::getExample_GV_7_2());
  TestBirch<Z64,4> test_7_3(BirchExample<Z64,4>::getExample_GV_7_3());
}