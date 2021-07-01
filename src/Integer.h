#ifndef __INTEGER_H_
#define __INTEGER_H_

#include <vector>

#include "polyhedral_common/src_number/NumberTheory.h"

#include "birch.h"
#include "IntegerRing.h"
#include "EuclideanDomainElement.h"

/**
 * An integer based on the integer class R.
 * This allows for both arbitrary-precision and finite-precision integers
 */

template<typename R>
class Integer : public virtual EuclideanDomainElement< Integer<R>, IntegerRing<R> >
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
  Integer<R> & operator=(const Integer<R> & b) override;

  // conversion
  // This might be dangerous, so we don't do it for now
  // inline operator R() const { return this->_num; }
  
  // access
  inline const R & num(void) const { return this->_num; }
  
  // arithmetic
  Integer<R> & operator+=(const Integer<R> &b) override;
  Integer<R> & operator-=(const Integer<R> &b) override;
  Integer<R> & operator*=(const Integer<R> &b) override;

  // !! - TODO - move these to a monoid / ordered / valued ring
  // comparison
  bool operator==(const Integer<R> &) const override;
  inline bool operator!=(const Integer<R> &b) const override {return !((*this)==b); }
  bool operator<(const Integer<R> &) const;
  inline bool operator>(const Integer<R> &b) const {return b < (*this); }
  inline bool operator<=(const Integer<R> &b) const
  {return ((*this) == b) || ((*this) < b); }
  inline bool operator>=(const Integer<R> &b) const
  {return ((*this) == b) || ((*this) > b); }

  // global constants

  inline bool isZero(void) const override { return (_num == 0); }

  // assign zero
  inline Integer<R> & makeZero(void) override { _num = 0; return (*this); }

  inline bool isOne(void) const override { return (_num == 1); }

  // assign to one
  inline Integer<R> & makeOne(void) override { _num = 1; return (*this); }

  inline static Integer<R> zero(void) { Integer<R> a; return a.makeZero(); }
  inline static Integer<R> one(void) { Integer<R> a; return a.makeOne(); }
  
  // euclidean division

  /**
   * Perform the eucldiean division of *this and b. Returns the
   * quotient and the remainder.
   */
  typename EuclideanDomainElement<Integer<R>, IntegerRing<R> >::DivRes
  euclideanDivision(const Integer<R>& b) const override;

  void print(std::ostream&) const override;

  inline Integer<R>* getPtr(void) override { return this; }

  inline const Integer<R>* getPtr(void) const override { return this; }

  inline std::shared_ptr<const IntegerRing<R> > parent(void) const override
  {return std::make_shared< IntegerRing<R> >(); }

  int hilbertSymbol(const Integer<R> & b, const Integer<R>& p) const;
  
  // temporary - to have it exported for where it is used 
  static int hilbertSymbolZ(const Z & x, const Z & y, const Z & p);

  typedef std::vector< std::pair<Integer<R>, size_t> > FactorData;

  typename Integer<R>::FactorData factorization(void) const;
  
  size_t valuation(const Integer<R>& p) const;

  inline Integer<R> abs() const {return _num < 0 ? -(*this) : *this; }

  int kroneckerSymbol(const Integer<R> & n) const;
  
  bool isLocalSquare(const Integer<R>& p) const;

  Integer<R> binomialCoefficient(const Integer<R> & k) const;

  Integer<R> nextPrime(void) const;

  bool isSquare(void) const;
  
protected:
  R _num;

};

template<typename R>
inline Integer<R>  operator*(const R& a,const Integer<R> &b)
{return b*a;}

namespace std
{
  template<typename R>
  struct hash< Integer<R> >
  {
    Z64 operator()(const Integer<R> & a) const
    {
      Z64 fnv = std::hash<R>{}(a.num());
      return fnv;
    }
  };
}

#include "Integer.inl"

#endif //__INTEGER_H_
