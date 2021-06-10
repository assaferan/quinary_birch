#ifndef __RATIONAL_H_
#define __RATIONAL_H_

#include "FieldElement.h"
#include "Integer.h"
#include "RationalField.h"

/**
 * A rational number based on the integer class R.
 * This allows for both arbitrary-precision and finite-precision rationals
 */

template<typename R>
class Rational : public virtual FieldElement< Rational<R>, RationalField<R> >
{
public:
  
  // c-tors
  Rational(const Integer<R>& num, const Integer<R>& denom = Integer<R>::one())
    : _num(num), _denom(denom)
  { reduce(); }

  Rational(const R & num, const R & denom = 1)
    : _num(num), _denom(denom)
  { reduce(); }

  // Rational(const Integer<R> & num) : _num(num), _denom(1) {}
  
  // default c-tor
  Rational() : _num(0), _denom(1) {}

  // copy c-tor
  Rational(const Rational<R> & other)
    : _num(other._num), _denom(other._denom) {}

  // access
  inline const Integer<R> & num() const {return this->_num; }
  inline const Integer<R> & denom() const {return this->_denom; }

  // arithmetic
  inline Rational<R> operator+() const {return Rational(_num, _denom); }
  inline Rational<R> operator-() const {return Rational(-_num, _denom); }
  Rational<R> operator+(const Rational<R> &) const override;
  inline Rational<R> operator-(const Rational<R> &b) const override
  {return (*this)+(-b); }
  Rational<R> operator*(const Rational<R> &) const override;
  inline Rational<R> operator*(const Integer<R> & b) const {
    Rational<R> b_rat(b);
    return (*this)*b_rat;
  }

  Rational<R> operator/(const Rational<R> &) const override;
  inline Rational<R> operator/(const Integer<R> & b) const {
    Rational<R> b_rat(b);
    return (*this)/b_rat;
  }
    
  // assignment
  inline Rational<R> & operator=(const Rational<R> & b)
  {
    if (this != &b) {
      _num = b._num;
      _denom = b._denom;
    }
    return (*this);
  }
  
  inline Rational<R> & operator+=(const Rational<R> &b)
  {return ((*this) = (*this) + b);}
  inline Rational<R> & operator-=(const Rational<R> &b)
  {return ((*this) = (*this) - b);}
  inline Rational<R> & operator*=(const Rational<R> &b)
  {return ((*this) = (*this) * b);}
  inline Rational<R> & operator/=(const Rational<R> &b) override
  {return ((*this) = (*this) / b);}

  inline Rational<R> inverse() const
  { Rational<R> inv(_denom, _num); return inv;}
  
  // comparison
  bool operator==(const Rational<R> &) const;
  
  bool operator<(const Rational<R> &) const;
  inline bool operator>(const Rational<R> &b) const {return b < (*this); }
  inline bool operator<=(const Rational<R> &b) const
  {return ((*this) == b) || ((*this) < b); }
  inline bool operator>=(const Rational<R> &b) const
  {return ((*this) == b) || ((*this) > b); }

  // other
  inline Integer<R> floor() const    
  { return _num / _denom; }

  inline Integer<R> ceiling() const
  { return (_num + _denom - 1)/_denom; }

  inline bool isIntegral() const
  {return ((_denom.isOne()) || ((-_denom).isOne())); }

  // other
  
  inline friend Rational<R> operator*(const Integer<R> & b, const Rational<R> & r)
  { return r*b; }
  /*
  inline friend Rational<R> operator-(const Integer<R> & b, const Rational<R> & r) {
    Rational<R> b_rat(b);
    return b_rat-r;
  }

  inline friend Rational<R> operator+(const Integer<R> & b, const Rational<R> & r) {
    Rational<R> b_rat(b);
    return b_rat+r;
  }

  inline friend Rational<R> operator/(const Integer<R> & b, const Rational<R> & r) {
    Rational<R> b_rat(b);
    return b_rat/r;
  }
  */
  
  inline void print(std::ostream & os) const
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

  inline static Rational<R> zero()
  { return Rational<R>::zero(); }
  
  inline static Rational<R> one()
  { return Rational<R>::one(); }
  
  inline Rational<R>* getPtr() { return this; }

  inline const Rational<R>* getPtr() const { return this; }

  inline std::shared_ptr<const RationalField<R> > parent() const override
  {return std::make_shared< RationalField<R> >(); }

  inline size_t valuation(const Integer<R>& p) const
  {return _num.valuation(p) - _denom.valuation(p);}

  inline Rational<R> abs() const
  {return (_num*_denom < Integer<R>::zero()) ? -(*this) : *this; }
  
protected:
  Integer<R> _num;
  Integer<R> _denom;

  void reduce(void);
};

#include "Rational.inl"

#endif // __RATIONAL_H_
