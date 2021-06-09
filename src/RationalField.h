#ifndef __RATIONAL_FIELD_H_
#define __RATIONAL_FIELD_H_

#include "Ring.h"

template<typename R>
class Rational;

template<typename R>
class RationalField : public virtual Ring< RationalField<R>, Rational<R> >
{
public:

  static RationalField<R>& getInstance()
  { static RationalField<R> instance; return instance; }

  RationalField(RationalField<R> const&) = delete;
  void operator=(RationalField<R> const &) = delete;
  
  inline Rational<R> zero() const override
  {Rational<R> z = 0; return z;}
  
  inline Rational<R> one() const override
  {Rational<R> z = 1; return z; }

  inline std::shared_ptr<const RationalField<R> > getPtr() const override
  {return std::enable_shared_from_this< const RationalField<R> >::shared_from_this(); }

private:
  RationalField() {}
};

#endif // __RATIONAL_FIELD_H_
