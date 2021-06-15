#ifndef __RATIONAL_FIELD_H_
#define __RATIONAL_FIELD_H_

#include "Ring.h"

template<typename R>
class Rational;

template<typename R>
class RationalField : public virtual Ring< RationalField<R>, Rational<R> >
{
public:

  RationalField() {}
  
  inline Rational<R> zero(void) const override
  {return Rational<R>::zero(); }
  
  inline Rational<R> one(void) const override
  {return Rational<R>::one(); }

  inline std::shared_ptr<const RationalField<R> > getPtr(void) const override
  {return std::enable_shared_from_this< const RationalField<R> >::shared_from_this(); }

  
};

#endif // __RATIONAL_FIELD_H_
