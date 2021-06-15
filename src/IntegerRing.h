#ifndef __INTEGER_RING_H_
#define __INTEGER_RING_H_

#include "Ring.h"

template<typename R>
class Integer;

template<typename R>
class IntegerRing : public virtual Ring< IntegerRing<R>, Integer<R> >
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

#endif // __INTEGER_RING_H_
