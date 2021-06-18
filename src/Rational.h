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
  { _reduce(); }

  Rational(const R & num, const R & denom = 1)
    : _num(num), _denom(denom)
  { _reduce(); }

  // Rational(const Integer<R> & num) : _num(num), _denom(1) {}
  
  // default c-tor
  Rational() : _num(0), _denom(1) {}

  // copy c-tor
  Rational(const Rational<R> & other)
    : _num(other._num), _denom(other._denom) {}

  // access
  inline const Integer<R> & num(void) const {return this->_num; }
  inline const Integer<R> & denom(void) const {return this->_denom; }

  // arithmetic
  inline Rational<R> operator+(void) const {return Rational(_num, _denom); }
  inline Rational<R> operator-(void) const override {return Rational(-_num, _denom); }
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
  inline Rational<R> & operator=(const Rational<R> & b) override
  {
    if (this != &b) {
      _num = b._num;
      _denom = b._denom;
    }
    return (*this);
  }
  
  inline Rational<R> & operator+=(const Rational<R> &b) override
  {return ((*this) = (*this) + b);}
  inline Rational<R> & operator-=(const Rational<R> &b) override
  {return ((*this) = (*this) - b);}
  inline Rational<R> & operator*=(const Rational<R> &b) override
  {return ((*this) = (*this) * b);}
  inline Rational<R> & operator/=(const Rational<R> &b) override
  {return ((*this) = (*this) / b);}

  inline Rational<R> inverse(void) const override
  { Rational<R> inv(_denom, _num); return inv;}
  
  // comparison
  bool operator==(const Rational<R> &) const override;
  
  bool operator<(const Rational<R> &) const;
  inline bool operator>(const Rational<R> &b) const {return b < (*this); }
  inline bool operator<=(const Rational<R> &b) const
  {return ((*this) == b) || ((*this) < b); }
  inline bool operator>=(const Rational<R> &b) const
  {return ((*this) == b) || ((*this) > b); }

  // other
  inline Integer<R> floor(void) const    
  { return _num / _denom; }

  inline Integer<R> ceiling(void) const
  { return (_num + _denom - Integer<R>::one())/_denom; }

  inline bool isIntegral(void) const
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
  
  inline void print(std::ostream & os) const override
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

  inline bool isZero(void) const override { return _num.isZero(); }

  // assign zero
  inline Rational<R> & makeZero(void) override { _num.makeZero(); return (*this); }

  inline bool isOne(void) const override { return (_num == _denom); }

  // assign to one
  inline Rational<R> & makeOne(void) override
  { _num.makeOne(); _denom.makeOne(); return (*this); }

  inline static Rational<R> zero(void) 
  { Rational<R> z; return z.makeZero(); }
  
  inline static Rational<R> one(void)
  { Rational<R> one; return one.makeOne(); }
  
  inline Rational<R>* getPtr(void) override { return this; }

  inline const Rational<R>* getPtr(void) const override { return this; }

  inline std::shared_ptr<const RationalField<R> > parent(void) const override
  {return std::make_shared< RationalField<R> >(); }

  inline size_t valuation(const Integer<R>& p) const
  {return _num.valuation(p) - _denom.valuation(p);}

  inline Rational<R> abs(void) const
  {return (_num*_denom < Integer<R>::zero()) ? -(*this) : *this; }

  inline bool isLocalSquare(const Integer<R>& p) const
  {return _num.isLocalSquare(p) == _denom.isLocalSquare(p); }

  static Rational<R> bernoulliNumber(size_t n);

  static Rational<R> bernoulliNumber(size_t n, const Integer<R> & d);

  void reduce(void);
  
protected:
  Integer<R> _num;
  Integer<R> _denom;

  
  static std::vector< Rational<R> > _bernoulliUpTo(size_t n);
  static std::vector< Rational<R> > _bernoulliPoly(size_t n);
};

namespace std
{
  template<typename R>
  struct hash< Rational<R> >
  {
    Z64 operator()(const Rational<R> & a) const
    {
      Z64 fnv = FNV_OFFSET;
      fnv = (fnv ^ std::hash< Integer<R> >{}(a.num())) * FNV_PRIME;
      fnv = (fnv ^ std::hash< Integer<R> >{}(a.denom())) * FNV_PRIME;
      return fnv;
    }
  };
}

#include "Rational.inl"

#endif // __RATIONAL_H_
