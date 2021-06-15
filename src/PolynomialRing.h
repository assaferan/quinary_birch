#ifndef __POLYNOMIAL_RING_H_
#define __POLYNOMIAL_RING_H_

#include "Ring.h"

template <typename R, typename S>
class PolynomialFp<R,S>;

template <typename R, typename S>
class PolynomialRingFp : public virtual Ring< PolynomialRingFp<R,S>, PolynomialFp<R,S> >
{
  public:

  IntegerRing() {}
  
  inline Integer<R> zero() const override
  {return Integer<R>::zero(); }
  
  inline Integer<R> one() const override
  {return Integer<R>::one(); }

  inline std::shared_ptr<const IntegerRing<R> > getPtr() const override
  {return std::enable_shared_from_this< const IntegerRing<R> >::shared_from_this(); }
};

#endif // __POLYNOMIAL_RING_H_
