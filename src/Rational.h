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
  Rational<R> operator+(const Rational<R> &) const;
  Rational<R> operator-(const Rational<R> &b) const {return (*this)+(-b); }
  Rational<R> operator*(const Rational<R> &) const;
  Rational<R> operator*(const Integer<R> & b) const {
    Rational<R> b_rat(b);
    return (*this)*b_rat;
  }

  Rational<R> operator/(const Rational<R> &) const;
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
  Rational<R> & operator/=(const Rational<R> &b)
  {return ((*this) = (*this) / b);}
  Rational<R> & operator/=(const R &b)
  {return ((*this) = (*this) / b);}

  // comparison
  bool operator==(const Rational<R> &) const;
  bool operator!=(const Rational<R> &b) const {return !((*this)==b); }
  bool operator<(const Rational<R> &) const;
  bool operator>(const Rational<R> &b) const {return b < (*this); }
  bool operator<=(const Rational<R> &b) const
  {return ((*this) == b) || ((*this) < b); }
  bool operator>=(const Rational<R> &b) const
  {return ((*this) == b) || ((*this) > b); }

  // other
  Integer<R> floor() const    
  {
    Integer<R> num = _num;
    Integer<R> denom = _denom;
    
    if (denom < 0) {
      denom = -denom;
      num = - num;
    }
      
    return ((num >= Integer<R>::zero()) ? num : (num - denom + 1)) / denom;
  }

  Integer<R> ceiling() const
  {
    Integer<R> num = _num;
    Integer<R> denom = _denom;
    
    if (denom < 0) {
      denom = -denom;
      num = - num;
    }
      
    return ((num >= Integer<R>::zero()) ? (num + denom - 1) : num) / denom;

  }

  bool is_integral() const
  {R one = 1; return ((_denom == one) || (_denom == -one)); }

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

  friend std::ostream& operator<<(std::ostream & os, const Rational<R> & r)
  {
    R one = 1;
    if (r.denom() == one) return os << r.num();
    if (r.denom() == -one) return os << -r.num();
    os << r.num() << "/" << r.denom();
    return os;
  }
  
protected:
  Integer<R> _num;
  Integer<R> _denom;

  void reduce(void);
};

#include "Rational.inl"

#endif // __RATIONAL_H_
