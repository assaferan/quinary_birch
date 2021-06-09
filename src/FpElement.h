#ifndef __FP_ELEMENT_H_
#define __FP_ELEMENT_H_

#include "Field.h"
#include "Fp.h"

template<typename R, typename S>
class FpElement : public Field< FpElement<R,S> >
{
public:
  //c-tors
  FpElement(std::shared_ptr<const Fp<R,S> > fld) : _GF(fld) {}
  
  FpElement(std::shared_ptr<const Fp<R,S> > fld, const R & val)
    : _GF(fld), _val(val) {}

  // access = get methods
  const R & lift() const {return _val; }
  const std::shared_ptr< const Fp<R,S> > & field() const {return _GF; }

  // arithmetic
  FpElement<R,S> operator+() const {return FpElement(_GF, _val); }
  FpElement<R,S> operator++(int)
  {(this->_val)++; return (*this); }
  FpElement<R,S> operator-() const override;
  FpElement<R,S> operator+(const FpElement<R,S> &other) const override;
  FpElement<R,S> operator-(const FpElement<R,S> &other) const override;
  FpElement<R,S> operator*(const FpElement<R,S> &other) const override;
  FpElement<R,S> operator/(const FpElement<R,S> &other) const override;
  FpElement<R,S> inverse(void) const override;
  
  FpElement<R,S> & operator+=(const FpElement<R, S> &other) override;
  FpElement<R,S> & operator-=(const FpElement<R, S> &other) override;
  FpElement<R,S> & operator*=(const FpElement<R, S> &other) override;
  FpElement<R,S> & operator/=(const FpElement<R, S> &other) override;
  
  FpElement<R, S> sqrt() const;
  // assignment and conversion
  FpElement<R, S> & operator=(const FpElement<R,S> &other) override;
  FpElement<R, S> & operator=(const R &other)
  { this->_val = other; return (*this); }
  
  //boolean
  bool isZero(void) const override
  {return (this->_val % this->_GF->prime() == 0);}

  bool isOne(void) const override
  {return (this->_val % this->_GF->prime() == 1);}
  
  bool operator==(const FpElement<R, S> &other) const override;
  
  // This is for sorting, we use the lift for that
  // !! TODO - Is it useful in anyway? shoul we get rid fof that
  bool operator<(const FpElement<R, S> &other) const
  { return (this->_val < other._val); }
  bool operator>(const FpElement<R, S> &other) const
  { return (this->_val > other._val); }
  bool operator>=(const FpElement<R, S> &other) const
  { return (this->_val >= other._val); }
  
  bool operator==(const R &other) const;
  bool operator!=(const R &other) const;
  bool isSquare(void) const;

  void setField(std::shared_ptr<const Fp<R,S>> fld) {this->_GF = fld;}

  FpElement<R,S>& makeZero() override {this->_val = 0;}
  FpElement<R,S>& makeOne() override {this->_val = 1;}
  
  FpElement<R,S>* getPtr() override {return this;}

  const FpElement<R,S>* getPtr() const override {return this;}
  
protected:
  std::shared_ptr< const Fp<R,S> > _GF;
  R _val;
  
};

#include "FpElement.inl"

#endif //__FP_ELEMENT_H
