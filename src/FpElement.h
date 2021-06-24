#ifndef __FP_ELEMENT_H_
#define __FP_ELEMENT_H_

#include "FieldElement.h"
#include "Fp.h"

template<typename R, typename S>
class FpElement : public virtual FieldElement< FpElement<R,S>, Fp<R,S> >
{
public:
  
  //c-tors
  // We allow a default constructor for static memory allocation
  // This might have consequence, so we should be extra careful.
  FpElement() = default;
  
  FpElement(std::shared_ptr<const Fp<R,S> > fld) : _GF(fld) {}
  
  FpElement(std::shared_ptr<const Fp<R,S> > fld, const R & val)
    : _GF(fld), _val(val) {}

  // access = get methods
  inline const R & lift(void) const {return _val; }
  inline R reducedLift(void) const {return _GF->mod(_val).lift(); }
  inline std::shared_ptr< const Fp<R,S> > parent(void) const override
  {return _GF; }

  // arithmetic
  inline FpElement<R,S> operator+() const {return FpElement(_GF, _val); }
  inline FpElement<R,S> operator++(int)
  {return ((*this) += _GF->one()); }
  FpElement<R,S> operator-() const override;
  FpElement<R,S> operator+(const FpElement<R,S> &other) const override;
  FpElement<R,S> operator-(const FpElement<R,S> &other) const override;
  FpElement<R,S> operator*(const FpElement<R,S> &other) const override;
  FpElement<R,S> operator/(const FpElement<R,S> &other) const override;
  FpElement<R,S> inverse(void) const override;
  
  FpElement<R,S> & operator+=(const FpElement<R,S> &other) override;
  FpElement<R,S> & operator-=(const FpElement<R,S> &other) override;
  FpElement<R,S> & operator*=(const FpElement<R,S> &other) override;
  FpElement<R,S> & operator/=(const FpElement<R,S> &other) override;

  int legendre(void) const;
  FpElement<R,S> sqrt(void) const;
  // assignment and conversion
  FpElement<R,S> & operator=(const FpElement<R,S> &other) override;
  inline FpElement<R,S> & operator=(const R &other)
  { this->_val = other; return (*this); }
  
  //boolean
  inline bool isZero(void) const override
  {return (this->_val % this->_GF->prime() == 0);}

  inline bool isOne(void) const override
  {return (this->_val % this->_GF->prime() == 1);}
  
  inline bool operator==(const FpElement<R,S> &other) const override;
  using RingElement< FpElement<R,S>, Fp<R,S> >::operator!=;
  
  // This is for sorting, we use the lift for that
  // !! TODO - Is it useful in anyway? shoul we get rid fof that
  inline bool operator<(const FpElement<R,S> &other) const
  { return (this->_GF->mod(this->_val).lift() < this->_GF->mod(other._val).lift()); }
  inline bool operator>(const FpElement<R,S> &other) const
  { return (other < (*this)); }
  inline bool operator>=(const FpElement<R,S> &other) const
  { return ((other < (*this)) || ((*this) == other)); }
  
  inline bool operator==(const R &other) const;
  inline bool operator!=(const R &other) const;
  inline bool isSquare(void) const;

  inline void setField(std::shared_ptr<const Fp<R,S>> fld) {this->_GF = fld;}

  inline FpElement<R,S>& makeZero() override {this->_val = 0; return (*this); }
  inline FpElement<R,S>& makeOne() override {this->_val = 1; return (*this); }

  inline static FpElement<R,S> zero(std::shared_ptr< const Fp<R,S> > GF)
  {FpElement<R,S> z(GF); z.makeZero(); return z;}

  inline static FpElement<R,S> one(std::shared_ptr< const Fp<R,S> > GF)
  {FpElement<R,S> z(GF); z.makeOne(); return z;}
  
  inline FpElement<R,S>* getPtr(void) override {return this;}

  inline const FpElement<R,S>* getPtr(void) const override {return this;}

  inline void print(std::ostream& os) const override
  {
    // making sure this is in the range [0,p) before printing
    os << (this->parent()->mod(this->lift())).lift();
    return;
  }
  
protected:
  std::shared_ptr< const Fp<R,S> > _GF;
  R _val;
  
};

#include "FpElement.inl"

#endif //__FP_ELEMENT_H
