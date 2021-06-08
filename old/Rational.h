#ifndef __RATIONAL_H_
#define __RATIONAL_H_

#include "Field.h"

template<typename R>
class Rational<R>;

template<typename R>
class Rationals : public Field
{
public:

  static Rationals<R>& getInstance();

  Rationals(Rationals<R> const&) = delete;
  void operator=(Rationals<R> const &) = delete;
  
  Rational<R> zero() const;
  Rational<R> one() const;

private:
  Rationals() {}
};

template<typename R>
class Rational : public FieldElement
{
public:
  
  // c-tors
  Rational(const R& num, const R& denom);

  Rational(const R & num);
  
  // default c-tor
  Rational();
  
  // copy c-tor
  Rational(const Rational<R> & other);

  // assignment
  Rational<R> & operator=(const Rational<R> & b);
  
  // access
  const R & num() const {return this->num_; }
  const R & denom() const {return this->denom_; }

  // arithmetic
  Rational<R> operator+() const {return Rational(num_, denom_); }
  Rational<R> operator-() const {return Rational(-num_, denom_); }
  Rational<R> operator+(const Rational<R> &) const;
  Rational<R> operator-(const Rational<R> &b) const {return (*this)+(-b); }
  Rational<R> operator*(const Rational<R> &) const;
  Rational<R> operator*(const R &) const;
  Rational<R> operator/(const Rational<R> &) const;
  Rational<R> operator/(const R &) const;
  
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
  R floor() const;

  R ceiling() const;
  
  bool is_integral() const
  {R one = 1; return ((denom_ == one) || (denom_ == -one)); }

protected:
  R num_;
  R denom_;

  void reduce(void);
};

// arithmetic
template <typename R>
Rational<R> operator*(R b, const Rational<R> & r);

template <typename R>
Rational<R> operator-(R b, const Rational<R> & r);

template <typename R>
Rational<R> operator+(R b, const Rational<R> & r);

template <typename R>
Rational<R> operator/(R b, const Rational<R> & r);

template<typename R>
Rational<R> abs(const Rational<R> & r)
{ R zero = 0; return (r > zero) ? r : -r;}

template <typename R>
std::ostream& operator<<(std::ostream &, const Rational<R> & );

#include "Rational.inl"

#endif // __RATIONAL_H_
