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
  TestBirch(const QuadFormZZ<R,n> &, ReductionMethod alg=GREEDY);

  TestBirch(const BirchExample<R,n> &, size_t, ReductionMethod alg=GREEDY);

  bool testDim(const R &, size_t) const;

  bool testEigenvalues(const R &,
		       const std::vector< std::map< R, std::vector< NumberFieldElement<Z> > > > &,
		       size_t) const;

  bool testEigenvalueTraces(const R &,
			    const std::vector< std::map< R, std::vector<R> > > &,
			    size_t) const;

  inline Genus<R,n> & genus(void) {return *_p_genus;}
  
protected:
  std::shared_ptr< Genus<R,n> > _p_genus;

  inline void _init(const QuadFormZZ<R,n> &, ReductionMethod alg);
  
};

struct EigenvalueVector
{
  inline const NumberFieldElement<Z> & operator[](size_t i) const
  {return vec[i]; }
  
  inline NumberFieldElement<Z> & operator[](size_t i)
  {return vec[i]; }

  inline size_t size(void) const
  { return vec.size(); }

  inline bool operator==(const EigenvalueVector & other) const
  {
    if (vec.size() != other.size()) return false;
    for (size_t i = 0; i < vec.size(); i++) {
      // could be isomorphic number fields. Right now just compare minimal polyomials
      // The right thing to do would be to loop over isomorphisms between the fields,
      // but it requires some more infrastructure.
      //      if (vec[i].parent() != other[i].parent()) return false;
      //if (vec[i] != other[i]) return false;
      if (vec[i].minimalPolynomial() != other[i].minimalPolynomial()) return false;
    }
    return true;
  }
  
  std::vector< NumberFieldElement<Z> > vec;
};

namespace std {
  
  template<>
  struct hash< EigenvalueVector >
  {
    Z64 operator()(const EigenvalueVector & vec) const
    {
      Z64 fnv = FNV_OFFSET;

      for (size_t i = 0; i < vec.size(); i++)
	fnv = (fnv ^ std::hash< NumberFieldElement<Z> >{}(vec[i])) * FNV_PRIME;
  
      return fnv;
    }
  };
}

void runQuinaryBirch(const Z64 & disc, size_t num_evs, ReductionMethod alg = GREEDY);
void runBirchTests(size_t num_evs = 0, ReductionMethod alg = GREEDY);

#include "TestBirch.inl"

#endif // __TEST_BIRCH_H_
