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
  const std::shared_ptr< const Fp<R, S> > & field() const {return _GF; }

  // arithmetic
  FpElement<R, S> operator+() const {return FpElement(_GF, _val); }
  FpElement<R, S> operator++(int)
  {(this->_val)++; return (*this); }
  FpElement<R, S> operator-() const;
  FpElement<R, S> operator+(const FpElement<R, S> &other) const;
  FpElement<R, S> operator-(const FpElement<R, S> &other) const;
  FpElement<R, S> operator*(const FpElement<R, S> &other) const;
  FpElement<R, S> operator/(const FpElement<R, S> &other) const;
  FpElement<R, S> inverse(void) const;
  
  FpElement<R, S> & operator+=(const FpElement<R, S> &other);
  FpElement<R, S> & operator-=(const FpElement<R, S> &other);
  FpElement<R, S> & operator*=(const FpElement<R, S> &other);
  FpElement<R, S> & operator/=(const FpElement<R, S> &other);
  
  FpElement<R, S> sqrt() const;
  // assignment and conversion
  FpElement<R, S> & operator=(const FpElement<R, S> &other);
  FpElement<R, S> & operator=(const R &other)
  { this->_val = other; return (*this); }
  
  //boolean
  bool isZero(void) const
  {return (this->_val % this->_GF->prime() == 0);}
  bool operator==(const FpElement<R, S> &other) const {
    assert( (_GF != 0) && (other._GF != 0) ); 
    if (this->_GF->prime() != other._GF->prime()) return false;
    return ((this->_val - other._val) % (this->_GF->prime()) == 0);
  }
  bool operator!=(const FpElement<R, S> &other) const
  { return !((*this)==other); }
  
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
  
protected:
  std::shared_ptr< const Fp<R, S> > _GF;
  R _val;
  
};

#include "FpElement.inl"

#endif //__FP_ELEMENT_H
