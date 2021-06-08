#ifndef __INTEGER_H
#define __INTEGER_H

#include "Ring.h"

template<typename R>
class Integer<R>;

template<typename R>
class Integers : public Ring
{
public:

  static Integers<R>& getInstance();

  Integers(Integers<R> const&) = delete;
  void operator=(Integers<R> const &) = delete;
  
  Integer<R> zero() const;
  Integer<R> one() const;

private:
  Integers() {}
};

template<typename R>
class Integer : public RingElement
{
public:
  
  // c-tors
  Integer(const R & );
  
  // default c-tor
  Integer();
  
  // copy c-tor
  Integer(const Integer<R> & other);

  // assignment
  Integer<R> & operator=(const Integer<R> & b);

  // access
  const R & num(void) const { return this->num_; }
  
  // arithmetic
  Integer<R> operator+() const {return Integer(num_); }
  Integer<R> operator-() const {return Integer(-num_); }
  Integer<R> operator+(const Integer<R> &) const;
  Integer<R> operator-(const Integer<R> &b) const {return (*this)+(-b); }
  Integer<R> operator*(const Integer<R> &) const;
  Integer<R> operator*(const R &) const;
  Integer<R> operator/(const Integer<R> &) const;
  Integer<R> operator/(const R &) const;
  Integer<R> operator%(const Integer<R> &) const;
  Integer<R> operator%(const R &) const;
  
  Integer<R> & operator+=(const Integer<R> &b)
  {return ((*this) = (*this) + b);}
  Integer<R> & operator-=(const Integer<R> &b)
  {return ((*this) = (*this) - b);}
  Integer<R> & operator*=(const Integer<R> &b)
  {return ((*this) = (*this) * b);}
  Integer<R> & operator*=(const R &b)
  {return ((*this) = (*this) * b);}
  Integer<R> & operator/=(const Integer<R> &b)
  {return ((*this) = (*this) / b);}
  Integer<R> & operator/=(const R &b)
  {return ((*this) = (*this) / b);}
  Integer<R> & operator%=(const Integer<R> &b)
  {return ((*this) = (*this) % b);}
  Integer<R> & operator%=(const R &b)
  {return ((*this) = (*this) % b);}

  // comparison
  bool operator==(const Integer<R> &) const;
  bool operator!=(const Integer<R> &b) const {return !((*this)==b); }
  bool operator<(const Integer<R> &) const;
  bool operator>(const Integer<R> &b) const {return b < (*this); }
  bool operator<=(const Integer<R> &b) const
  {return ((*this) == b) || ((*this) < b); }
  bool operator>=(const Integer<R> &b) const
  {return ((*this) == b) || ((*this) > b); }

protected:
  R num_;

};

// arithmetic
template <typename R>
Integer<R> operator*(R b, const Integer<R> & r);

template <typename R>
Integer<R> operator-(R b, const Integer<R> & r);

template <typename R>
Integer<R> operator+(R b, const Integer<R> & r);

template <typename R>
Integer<R> operator/(R b, const Integer<R> & r);

template <typename R>
Integer<R> operator%(R b, const Integer<R> & r);

template<typename R>
Integer<R> abs(const Integer<R> & r)
{ R zero = 0; return (r > zero) ? r : -r;}

template <typename R>
std::ostream& operator<<(std::ostream &, const Integer<R> & );

#include "Integer.inl"

#endif //__INTEGER_H
