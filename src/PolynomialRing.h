#ifndef __POLYNOMIAL_RING_H_
#define __POLYNOMIAL_RING_H_

#include "Ring.h"

template <typename R, typename S>
class PolynomialFp<R,S>;

template <typename R, typename S>
class PolynomialRingFp : public virtual Ring< PolynomialRingFp<R,S>, PolynomialFp<R,S> >
{
public:

  PolynomialRingFp(std::shared_ptr<const Fp<R,S> > GF) : _base(GF) {}
  
  inline PolynomialFp<R,S> zero(void) const override
  {return PolynomialFp<R,S>::zero(_base); }
  
  inline PolynomialFp<R,S> one(void) const override
  {return PolynomialFp<R,S>::one(_base); }

  inline std::shared_ptr<const PolynomialRingFp<R,S> > getPtr(void) const override
  {return std::enable_shared_from_this< const PolynomialRingFp<R,S> >::shared_from_this(); }

protected:
  std::shared_ptr<const Fp<R,S> > _base;
};

#endif // __POLYNOMIAL_RING_H_
