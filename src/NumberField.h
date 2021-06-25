#ifndef __NUMBER_FIELD_H_
#define __NUMBER_FIELD_H_

#include "antic/nf.h"

#include "birch.h"

template<typename R>
class NumberField : public virtual Ring< NumberField<R>, NumberFieldElement<R> >
{
public:
  NumberField(const UnivariatePolyRat<R> & mod) : _f(mod)
  {_initAntic(); }
  
  NumberField(const UnivariatePolyInt<R> &);

  inline const UnivariatePolyRat<R> & modulus(void) const
  {return this->_f; }

  inline std::shared_ptr< const NumberField<R> > getPtr(void) const override
  { return std::enable_shared_from_this< const NumberField<R> >::shared_from_this(); }

  inline NumberFieldElement<R> zero(void) const override
  {return NumberFieldElement<R>::zero(this->getPtr()); }
  
  inline NumberFieldElement<R> one(void) const override
  {return NumberFieldElement<R>::one(this->getPtr()); }

  ~NumberField()
  { nf_clear(_nf_antic); }
  
protected:
  UnivariatePolyRat<R> _f;
  nf_t _nf_antic;

  void _initAntic(void);
};

#include "NumberField.inl"

#endif // _NUMBER_FIELD_H_
