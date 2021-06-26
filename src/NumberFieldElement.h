#ifndef __NUMBER_FIELD_ELEMENT_H_
#define __NUMBER_FIELD_ELEMENT_H_

#include "antic/nf_elem.h"

#include "birch.h"
#include "FieldElement.h"
#include "UnivariatePoly.h"

template<typename R>
class NumberFieldElement : public virtual FieldElement< NumberFieldElement<R>, NumberField<R> >
{
public:
  NumberFieldElement() = default;
  NumberFieldElement(std::shared_ptr<const NumberField<R> > fld) : _K(fld) { _initAntic(); }
  NumberFieldElement(std::shared_ptr<const NumberField<R> > fld,
		     const UnivariatePolyRat<R> & poly)
    : _K(fld) {_initAntic(poly); }
  NumberFieldElement(std::shared_ptr<const NumberField<R> > fld,
		     const R & a)
    : _K(fld) { _initAntic(a); }
  NumberFieldElement(std::shared_ptr<const NumberField<R> > fld,
		     const Integer<R> & a)
    : _K(fld) { _initAntic(a); } 
  NumberFieldElement(std::shared_ptr<const NumberField<R> > fld,
		     const Rational<R> & a) : _K(fld) { _initAntic(a); }

  NumberFieldElement(const NumberFieldElement<R> &);

  // assignment operator
  NumberFieldElement<R> & operator=(const NumberFieldElement<R> &) override;
  NumberFieldElement<R> & operator=(const R &);
  
  // arithmetic
  NumberFieldElement<R> operator-() const override;
  NumberFieldElement<R> operator+(const NumberFieldElement<R> & ) const override;
  NumberFieldElement<R> operator-(const NumberFieldElement<R> & ) const override;
  NumberFieldElement<R> operator*(const NumberFieldElement<R> & ) const override;
  NumberFieldElement<R> operator/(const NumberFieldElement<R> & ) const override;

  NumberFieldElement<R> operator*(const R & ) const;
  NumberFieldElement<R> operator/(const R & ) const;

  NumberFieldElement<R> & operator+=(const NumberFieldElement<R> & ) override;
  NumberFieldElement<R> & operator-=(const NumberFieldElement<R> & ) override;
  NumberFieldElement<R> & operator*=(const NumberFieldElement<R> & ) override;
  NumberFieldElement<R> & operator/=(const NumberFieldElement<R> & ) override;

  NumberFieldElement<R>& operator*=(const R & );
  NumberFieldElement<R>& operator/=(const R & );

  NumberFieldElement<R> inverse(void) const override;

  inline std::shared_ptr< const NumberField<R> > parent(void) const override
  {return this->_K;}

  UnivariatePolyInt<R> minimalPolynomial(void) const;

  bool isZero(void) const override;
  bool isOne(void) const override;

  static NumberFieldElement<R> zero(std::shared_ptr<const NumberField<R> > fld);
  static NumberFieldElement<R> one(std::shared_ptr<const NumberField<R> > fld);
  
  inline NumberFieldElement<R>& makeZero(void) override
  { nf_elem_zero(_nf_elt_antic, _K->antic()); return *this; }

  inline NumberFieldElement<R>& makeOne(void) override
  { nf_elem_one(_nf_elt_antic, _K->antic()); return *this; }

  inline NumberFieldElement<R>* getPtr(void) override {return this;}
  inline const NumberFieldElement<R>* getPtr(void) const override {return this;}

  inline void print(std::ostream& os) const override
  { char* elt_str = nf_elem_get_str_pretty(_nf_elt_antic, "x", _K->antic());
    os << elt_str; return; }

  inline const nf_elem_t & getPoly(void) const
  { return _nf_elt_antic; }

  inline R trace(void) const
  { R trace;
    Z num, denom;
    fmpq_t ftrace;
    fmpq_init(ftrace);
    nf_elem_trace(ftrace, _nf_elt_antic, _K->antic());
    fmpq_get_mpz_frac(num.get_mpz_t(), denom.get_mpz_t(), ftrace);
    fmpq_clear(ftrace);
    
    return trace;
  }

  inline R norm(void) const
  { R norm;
    Z num, denom;
    fmpq_t fnorm;
    fmpq_init(fnorm);
    nf_elem_norm(fnorm, _nf_elt_antic, _K->antic());
    fmpq_get_mpz_frac(num.get_mpz_t(), denom.get_mpz_t(), fnorm);
    fmpq_clear(fnorm);
    
    return norm;
  }

  ~NumberFieldElement()
  {
    // apparently when the element is zero, antic tries to free the polynomial
    // that was not allocated
    if (!nf_elem_is_zero(_nf_elt_antic, _K->antic()))
      nf_elem_clear(_nf_elt_antic, _K->antic());
  }
  
protected:
  std::shared_ptr<const NumberField<R> > _K;
  
  nf_elem_t _nf_elt_antic;

  MatrixRat<R> _multByMatrix(void) const;
  void _initAntic(void);
};

namespace std
{
  template<typename R>
  struct hash< NumberFieldElement<R> >
  {
    Z64 operator()(const NumberFieldElement<R>& elt) const
    {
      Z64 fnv = std::hash< UnivariatePolyInt<R> >{}(elt.minimalPolynomial());
            
      return fnv;
    }
  };
  
}

#include "NumberFieldElement.inl"

#endif // __NUMBER_FIELD_ELEMENT_H_
