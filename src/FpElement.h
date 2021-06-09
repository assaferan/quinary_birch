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
  inline const R & lift() const {return _val; }
  inline std::shared_ptr< const Fp<R,S> > parent() const override
  {return _GF; }

  // arithmetic
  inline FpElement<R,S> operator+() const {return FpElement(_GF, _val); }
  inline FpElement<R,S> operator++(int)
  {(this->_val)++; return (*this); }
  inline FpElement<R,S> operator-() const override;
  inline FpElement<R,S> operator+(const FpElement<R,S> &other) const override;
  inline FpElement<R,S> operator-(const FpElement<R,S> &other) const override;
  inline FpElement<R,S> operator*(const FpElement<R,S> &other) const override;
  inline FpElement<R,S> operator/(const FpElement<R,S> &other) const override;
  inline FpElement<R,S> inverse(void) const override;
  
  inline FpElement<R,S> & operator+=(const FpElement<R, S> &other) override;
  inline FpElement<R,S> & operator-=(const FpElement<R, S> &other) override;
  inline FpElement<R,S> & operator*=(const FpElement<R, S> &other) override;
  inline FpElement<R,S> & operator/=(const FpElement<R, S> &other) override;

  inline int legendre() const;
  inline FpElement<R, S> sqrt() const;
  // assignment and conversion
  inline FpElement<R,S> & operator=(const FpElement<R,S> &other) override;
  inline FpElement<R,S> & operator=(const R &other)
  { this->_val = other; return (*this); }
  
  //boolean
  inline bool isZero(void) const override
  {return (this->_val % this->_GF->prime() == 0);}

  inline bool isOne(void) const override
  {return (this->_val % this->_GF->prime() == 1);}
  
  inline bool operator==(const FpElement<R, S> &other) const override;
  
  // This is for sorting, we use the lift for that
  // !! TODO - Is it useful in anyway? shoul we get rid fof that
  inline bool operator<(const FpElement<R, S> &other) const
  { return (this->_val < other._val); }
  inline bool operator>(const FpElement<R, S> &other) const
  { return (this->_val > other._val); }
  inline bool operator>=(const FpElement<R, S> &other) const
  { return (this->_val >= other._val); }
  
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
  
  inline FpElement<R,S>* getPtr() override {return this;}

  inline const FpElement<R,S>* getPtr() const override {return this;}

  inline void print(std::ostream& os) const
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
