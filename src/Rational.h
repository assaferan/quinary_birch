#ifndef __RATIONAL_H_
#define __RATIONAL_H_

#include "Field.h"
#include "Integer.h"

/**
 * A rational number based on the integer class R.
 * This allows for both arbitrary-precision and finite-precision rationals
 */

template<typename R>
class Rational : public Field< Rational<R> >
{
public:
  
  // c-tors
  Rational(const Integer<R>& num, const Integer<R>& denom)
    : _num(num), _denom(denom)
  { reduce(); }

  Rational(const R & num, const R & denom = 1)
    : _num(num), _denom(denom)
  { reduce(); }

  Rational(const Integer<R> & num) : _num(num), _denom(1) {}
  
  // default c-tor
  Rational() : _num(0), _denom(1) {}

  // copy c-tor
  Rational(const Rational<R> & other)
    : _num(other._num), _denom(other._denom) {}

  // access
  const Integer<R> & num() const {return this->_num; }
  const Integer<R> & denom() const {return this->_denom; }

  // arithmetic
  Rational<R> operator+() const {return Rational(_num, _denom); }
  Rational<R> operator-() const {return Rational(-_num, _denom); }
  Rational<R> operator+(const Rational<R> &) const override;
  Rational<R> operator-(const Rational<R> &b) const override
  {return (*this)+(-b); }
  Rational<R> operator*(const Rational<R> &) const override;
  Rational<R> operator*(const Integer<R> & b) const {
    Rational<R> b_rat(b);
    return (*this)*b_rat;
  }

  Rational<R> operator/(const Rational<R> &) const override;
  Rational<R> operator/(const Integer<R> & b) const {
    Rational<R> b_rat(b);
    return (*this)/b_rat;
  }
    
  // assignment
  Rational<R> & operator=(const Rational<R> & b)
  {
    if (this != &b) {
      _num = b._num;
      _denom = b._denom;
    }
    return (*this);
  }
  
  Rational<R> & operator+=(const Rational<R> &b)
  {return ((*this) = (*this) + b);}
  Rational<R> & operator-=(const Rational<R> &b)
  {return ((*this) = (*this) - b);}
  Rational<R> & operator*=(const Rational<R> &b)
  {return ((*this) = (*this) * b);}
  Rational<R> & operator*=(const R &b)
  {return ((*this) = (*this) * b);}
  Rational<R> & operator/=(const Rational<R> &b) override
  {return ((*this) = (*this) / b);}
  Rational<R> & operator/=(const R &b)
  {return ((*this) = (*this) / b);}

  Rational<R> inverse() const
  { Rational<R> inv(_denom, _num); return inv;}
  
  // comparison
  bool operator==(const Rational<R> &) const;
  
  bool operator<(const Rational<R> &) const;
  bool operator>(const Rational<R> &b) const {return b < (*this); }
  bool operator<=(const Rational<R> &b) const
  {return ((*this) == b) || ((*this) < b); }
  bool operator>=(const Rational<R> &b) const
  {return ((*this) == b) || ((*this) > b); }

  // other
  Integer<R> floor() const    
  { return _num / _denom; }

  Integer<R> ceiling() const
  { return (_num + _denom - 1)/_denom; }

  bool isIntegral() const
  {return ((_denom.isOne()) || ((-_denom).isOne())); }

  // other
  friend Rational<R> operator*(Integer<R> b, const Rational<R> & r) {
    return r*b;
  }

  friend Rational<R> operator-(Integer<R> b, const Rational<R> & r) {
    Rational<R> b_rat(b);
    return b_rat-r;
  }

  friend Rational<R> operator+(Integer<R> b, const Rational<R> & r) {
    Rational<R> b_rat(b);
    return b_rat+r;
  }

  friend Rational<R> operator/(Integer<R> b, const Rational<R> & r) {
    Rational<R> b_rat(b);
    return b_rat/r;
  }

  void print(std::ostream & os) const
  {
    if (_denom.isOne())
      os << _num;
    else if ((-_denom).isOne())
      os << -_num;
    else
      os << _num << "/" << _denom;
    return;
  }

  // zero and one global constants

  inline bool isZero() const { return _num.isZero(); }

  // assign zero
  inline Rational<R> & makeZero() { _num.makeZero(); return (*this); }

  inline bool isOne() const { return (_num == _denom); }

  // assign to one
  inline Rational<R> & makeOne()
  { _num.makeOne(); _denom.makeOne(); return (*this); }

  static Rational<R> zero() { Rational<R> a; return a.makeZero(); }
  static Rational<R> one() { Rational<R> a; return a.makeOne(); }
  
  Rational<R>* getPtr() { return this; }

  const Rational<R>* getPtr() const { return this; }
  
protected:
  Integer<R> _num;
  Integer<R> _denom;

  void reduce(void);
};

#include "Rational.inl"

#endif // __RATIONAL_H_
