#include <chrono>
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
TestBirch<R,n>::TestBirch(const QuadFormZZ<R,n> & q)
{
  this->_init(q);
  std::cerr << "Testing orthogonal modular forms for " << std::endl  << q << std::endl;
}

template<typename R, size_t n>
inline bool TestBirch<R,n>::testDim(const R & spinor_prime, size_t dim) const
{
  std::map<R,size_t> dims = this->_p_genus->dimensionMap();

  return (dims[spinor_prime] == dim);
}

template<typename R, size_t n>
inline bool TestBirch<R,n>::testEigenvalues(const R & spinor_prime,
					    const  std::vector< std::map< R, std::vector< NumberFieldElement<Z> > > > & evs,
					    size_t num_evs) const
{
  EigenvectorManager<R,n> manager;
  std::vector< std::vector< NumberFieldElement<Z> > > evecs = _p_genus->eigenvectors()[spinor_prime];
  for (std::vector< NumberFieldElement<Z> > evec : evecs)
      manager.addEigenvector(_p_genus->eigenvector(evec, spinor_prime));
  manager.finalize();

  std::vector< EigenvalueVector > evalues(evecs.size());
  std::vector< EigenvalueVector > computed_evalues(evecs.size());

  std::cerr << "Testing Hecke eigensystem with spinor = " << spinor_prime << std::endl;
  
  for (size_t k = 0; k < evs.size(); k++) {
    size_t num_processed = 0;
    for (std::pair< R, std::vector< NumberFieldElement<Z> > > ev : evs[k]) {
      Integer<R> p = ev.first;
      std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
      std::vector< NumberFieldElement<Z> > computed = _p_genus->eigenvalues(manager, p.num(), k+1);
      std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
      std::cerr << "computing eigenvalues took " << std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count() << " ms\n";
      for (size_t i = 0; i < evecs.size(); i++) {
	evalues[i].vec.push_back(ev.second[i]);
	computed_evalues[i].vec.push_back(computed[i]);
      }
      // #ifdef DEBUG
      std::cerr << "Testing eigenvalues of T_" << p << "^" << k << "..." << std::endl;
      std::cerr << "ev.second = " << ev.second << std::endl;
      std::cerr << "computed eigenvalues: " << computed << std::endl;
      // #endif
      num_processed++;
      if (num_processed == num_evs) break;
      //    assert(ev.second == _p_genus->eigenvalues(manager, p.num()));
    }

    std::unordered_set< EigenvalueVector > evalue_set(evalues.begin(), evalues.end());
    std::unordered_set< EigenvalueVector > computed_set(computed_evalues.begin(), computed_evalues.end());

    if ( evalue_set != computed_set)
      return false;
  }
  
  return true;
}

template<typename R, size_t n>
inline bool TestBirch<R,n>::testEigenvalueTraces(const R & spinor_prime,
						 const  std::vector< std::map< R, std::vector<R> > > & traces,
						 size_t num_evs) const
{
  EigenvectorManager<R,n> manager;
  std::vector< std::vector< NumberFieldElement<Z> > > evecs = _p_genus->eigenvectors()[spinor_prime];
  for (std::vector< NumberFieldElement<Z> > evec : evecs)
      manager.addEigenvector(_p_genus->eigenvector(evec, spinor_prime));
  manager.finalize();

  std::vector< std::vector<R> > evalues(evecs.size());
  std::vector< std::vector<R> > computed_evalues(evecs.size());

  for (size_t k = 0; k < traces.size(); k++) {
    size_t num_processed = 0;
    for (std::pair< R, std::vector<R> > ev_tr : traces[k]) {
      Integer<R> p = ev_tr.first;
      std::vector< NumberFieldElement<Z> > computed = _p_genus->eigenvalues(manager, p.num(), k+1);
      for (size_t i = 0; i < evecs.size(); i++) {
	evalues[i].vec.push_back(ev_tr.second[i]);
	computed_evalues[i].vec.push_back(birch_util::convertInteger<R,Z>(computed[i].trace()));
      }
      // #ifdef DEBUG
      std::cerr << "ev_tr.second = " << ev_tr.second << std::endl;
      std::cerr << "computed eigenvalues: " << computed << std::endl;
      // #endif
      num_processed++;
      if (num_processed == num_evs) break;
      //    assert(ev.second == _p_genus->eigenvalues(manager, p.num()));
    }

    std::unordered_set< std::vector<R> > evalue_set(evalues.begin(), evalues.end());
    std::unordered_set< std::vector<R> > computed_set(computed_evalues.begin(), computed_evalues.end());

    if ( evalue_set != computed_set)
      return false;
  }
  
  return true;
}

template<typename R, size_t n>
inline TestBirch<R,n>::TestBirch(const BirchExample<R,n> & example, size_t num_evs)
{
  this->_init(example.qf);
  bool pass_dim = this->testDim(example.spinor_prime, example.dim);
  if (!pass_dim)
    throw std::runtime_error("Dimension test failed.\n");
  else {
    bool pass_ev = this->testEigenvalues(example.spinor_prime, example.evs, num_evs);
    if (!pass_ev)
      throw std::runtime_error("Eigenvalue test failed.\n");
  }
  
}

template<typename R, size_t n>
inline void TestBirch<R,n>::_init(const QuadFormZZ<R,n> & q)
{
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

inline void runBirchTests(size_t num_evs) {
  // tetss reading from Nipp's tables
  std::vector< std::vector< QuadFormZZ<Z64,5> > > qfs;
  qfs = QuadFormZZ<Z64,5>::getQuinaryForms(256);
  qfs = QuadFormZZ<Z64,5>::getQuinaryForms(300);

  TestBirch<Z64,3> test_7_2(BirchExample<Z64,3>::getExample_GV_7_2(), num_evs);
  TestBirch<Z64,3> test_cmf_49_2_a_a(BirchExample<Z64,3>::getExample_CMF_49_2_a_a(), num_evs);
  TestBirch<Z64,4> test_7_3(BirchExample<Z64,4>::getExample_GV_7_3(), num_evs);
  TestBirch<Z64,5> test_RT_table1(BirchExample<Z64,5>::getExample_RT_Table1(), num_evs);
  
}
