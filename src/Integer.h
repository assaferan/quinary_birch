#ifndef __INTEGER_H_
#define __INTEGER_H_

#include "IntegerRing.h"
#include "EuclideanDomainElement.h"

/**
 * An integer based on the integer class R.
 * This allows for both arbitrary-precision and finite-precision integers
 */

template<typename R>
class Integer : public virtual EuclideanDomainElement< Integer<R> >
{
public:
  
  // c-tors
  Integer(const R & num)
    : _num(num) {}
  
  /**
   * Get a zero integer.
   */
  Integer()
    : _num(0) {}
  
  // copy c-tor
  Integer(const Integer<R> & other)
    : _num(other._num) {}
  
  // assignment
  Integer<R> & operator=(const Integer<R> & b);

  // access
  inline const R & num(void) const { return this->_num; }
  
  // arithmetic
  inline Integer<R> & operator+=(const Integer<R> &b);
  inline Integer<R> & operator-=(const Integer<R> &b);
  inline Integer<R> & operator*=(const Integer<R> &b);

  // !! - TODO - move these to a monoid / ordered / valued ring
  // comparison
  inline bool operator==(const Integer<R> &) const;
  inline bool operator!=(const Integer<R> &b) const {return !((*this)==b); }
  inline bool operator<(const Integer<R> &) const;
  inline bool operator>(const Integer<R> &b) const {return b < (*this); }
  inline bool operator<=(const Integer<R> &b) const
  {return ((*this) == b) || ((*this) < b); }
  inline bool operator>=(const Integer<R> &b) const
  {return ((*this) == b) || ((*this) > b); }

  // global constants

  inline bool isZero() const { return (_num == 0); }

  // assign zero
  inline Integer<R> & makeZero() { _num = 0; return (*this); }

  inline bool isOne() const { return (_num == 1); }

  // assign to one
  inline Integer<R> & makeOne() { _num = 1; return (*this); }

  inline static Integer<R> zero() { Integer<R> a; return a.makeZero(); }
  inline static Integer<R> one() { Integer<R> a; return a.makeOne(); }
  
  // euclidean division

  /**
   * Perform the eucldiean division of *this and b. Returns the
   * quotient and the remainder.
   */
  inline typename EuclideanDomainElement<Integer<R> >::DivRes
  euclideanDivision(const Integer<R>& b) const;

  inline void print(std::ostream&) const;

  inline Integer<R>* getPtr() { return this; }

  inline const Integer<R>* getPtr() const { return this; }

  inline std::shared_ptr<const Ring< Integer<R> > > parent() const
  {return std::make_shared(IntegerRing<R>::getInstance()); }
  
protected:
  R _num;

};

#include "Integer.inl"

#endif //__INTEGER_H_
