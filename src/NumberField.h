#ifndef __NUMBER_FIELD_H_
#define __NUMBER_FIELD_H_

#include "birch.h"
#include "nf.h"

template<typename R>
class NumberField : public virtual Ring< NumberField<R>, NumberFieldElement<R> >
{
public:
  NumberField(const UnivariatePolyRat<R> & mod) : _f(mod)
  {
    fmpq_poly_t poly;
    fmpq_poly_init(poly);
    assert(!mod.isZero());
      
    mpq_t* c_mpq = new mpq_t(mod.degree()+1);
    
    for (int i = 0; i <= mod.degree(); i++) {
      Rational<R> c = mod.coefficient(i);
      fmpq_set_si(c_mpq[i], c.num().num(), c.denom().num());
    }
    
    fmpq_poly_set_array_mpq(poly, c_mpq, mod.degree()+1);
    nf_init(_nf_antic, poly);

    fmpq_poly_clear(poly);
    delete[] c_mpq;
  }
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
  {
    nf_clear(_nf_antic);
  }
  
protected:
  UnivariatePolyRat<R> _f;
  nf_t _nf_antic;

  void _initAntic(void);
};

#include "NumberField.inl"

#endif // _NUMBER_FIELD_H_
