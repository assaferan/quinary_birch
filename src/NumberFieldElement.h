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
  NumberFieldElement(std::shared_ptr<const NumberField<R> > fld) : _K(fld), _elt(std::make_shared< const RationalField<R> >()) { _initAntic(); }
  NumberFieldElement(std::shared_ptr<const NumberField<R> > fld,
		     const UnivariatePolyRat<R> & poly)
    : _K(fld), _elt(poly) {_initAntic(); }
  NumberFieldElement(std::shared_ptr<const NumberField<R> > fld,
		     const R & a)
    : _K(fld), _elt(a) { _initAntic(); }
  NumberFieldElement(std::shared_ptr<const NumberField<R> > fld,
		     const Integer<R> & a)
    : _K(fld), _elt(a) { _initAntic(); } 
  NumberFieldElement(std::shared_ptr<const NumberField<R> > fld,
		     const Rational<R> & a) : _K(fld), _elt(a) { _initAntic(); }

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
  { _elt.makeZero(); return *this; }

  inline NumberFieldElement<R>& makeOne(void) override
  { _elt.makeOne(); return *this; }

  inline NumberFieldElement<R>* getPtr(void) override {return this;}
  inline const NumberFieldElement<R>* getPtr(void) const override {return this;}

  inline void print(std::ostream& os) const override
  { _elt.print(os); return; }

  inline const UnivariatePolyRat<R> & getPoly(void) const
  { return _elt; }

  inline R trace(void) const
  { UnivariatePolyInt<R> f = this->minimalPolynomial(); return f.coefficient(f.degree()-1); }

  inline R norm(void) const
  { UnivariatePolyInt<R> f = this->minimalPolynomial(); R sign = (f.degree() % 2 == 0) ? 1 : -1; return sign*f.coefficient(0);}

  ~NumberFieldElement()
  { nf_elem_clear(_nf_elt_antic); }
  
protected:
  std::shared_ptr<const NumberField<R> > _K;
  UnivariatePolyRat<R> _elt;
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
