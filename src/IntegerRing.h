#ifndef __INTEGER_RING_H
#define __INTEGER_RING_H

#include "Ring.h"

template<typename R>
class Integer;

template<typename R>
class IntegerRing : public virtual Ring< IntegerRing<R>, Integer<R> >
{
public:

  inline static IntegerRing<R>& getInstance()
  { static IntegerRing<R> instance; return instance; }

  IntegerRing(IntegerRing<R> const&) = delete;
  void operator=(IntegerRing<R> const &) = delete;
  
  inline Integer<R> zero() const override
  {return Integer<R>::zero(); }
  
  inline Integer<R> one() const override
  {return Integer<R>::one(); }

  inline std::shared_ptr<const IntegerRing<R> > getPtr() const override
  {return std::enable_shared_from_this< const IntegerRing<R> >::shared_from_this(); }

private:
  IntegerRing() {}
};

#endif // __INTEGER_RING_H
