#ifndef __NUMBER_FIELD_H_
#define __NUMBER_FIELD_H_

#include "birch.h"

template<typename R>
class NumberField : public virtual Ring< NumberField<R>, NumberFieldElement<R,S> >
{
public:
  NumberField(const UnivariatePolyRat<R> & mod) : _f(mod) {}
  NumberField(const UnivariatePolyInt<R> &);

  inline const UnivariatePolyRat<R> & modulus() const
  {return this->_f; }

  inline std::shared_ptr< const NumberField<R> > getPtr() const override
  { return std::enable_shared_from_this< const NumberField<R> >::shared_from_this(); }

  inline NumberField<R> zero(void) const override
  {return NumberFieldElement<R>::zero(this->getPtr()); }
  
  inline NumberField<R> one(void) const override
  {return NumberFieldElement<R>::one(this->getPtr()); }
  
protected:
  UnivariatePolyRat<R> _f;
};

#include "NumberField.inl"

#endif // _NUMBER_FIELD_H
