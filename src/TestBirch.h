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
  TestBirch(const typename QuadFormZZ<R,n>::SymVec &);

  TestBirch(const BirchExample<R,n> &, size_t);

  void testDim(const R &, size_t) const;

  void testEigenvalues(const R &,
		       const std::vector< std::map< R, std::vector< NumberFieldElement<Z> > > > &,
		       size_t) const;
  
protected:
  std::shared_ptr< Genus<R,n> > _p_genus;

  inline void _init(const typename QuadFormZZ<R,n>::SymVec &);
  
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
    for (size_t i = 0; i < vec.size(); i++)
      if (vec[i] != other[i]) return false;
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

void runBirchTests(size_t num_evs = 0);

#include "TestBirch.inl"

#endif // __TEST_BIRCH_H_
