#ifndef __EUCLIDEAN_DOMAIN_ELEMENT_H_
#define __EUCLIDEAN_DOMAIN_ELEMENT_H_

#include <tuple>

#include "RingElement.h"

/**
 * An abstract class defining the interface of a Euclidean domain.
 */

template<class Derived>
class EuclideanDomainElement : public virtual RingElement<Derived>
{
  public:

  typedef std::pair<Derived, Derived> DivRes;
  typedef std::tuple<Derived, Derived, Derived> XGcdRes;
  
  /**
   * Perform the eucldiean division of *this and b. Returns the
   * quotient and the remainder.
   *
   * @param b: the divisor.
   *
   * @return the quotient and the remainder.
   */
  virtual DivRes euclideanDivision(const Derived& b) const = 0;

  /**
   * Perform the extended euclidean division on *this and b.
   * Returns the GCD.
   * @param b: the divisor.
   * @return the gcd and the bezout coefficients.
   */
  inline virtual XGcdRes xgcd(const Derived& b) const;

  /**
   * Get GCD of *this and other.
   *
   * @param other: the other element to get a gcd with.
   * @return the gcd.
   */
  inline virtual Derived gcd(const Derived& other) const;
  
  /**
   * Get the quotient of *this and b.
   *
   * @param b: the divisor
   * @return the quotient
   */
  inline virtual Derived quotient(const Derived& b) const;

  /**
   * Get the remainder of *this and b.
   * @param b: the divisor
   * @return the remainder
   */
  inline virtual Derived remainder(const Derived& b) const;

  /**
   * Exact division.
   *
   * @param d: the divisor.
   * @return the equotient.
   */
  inline virtual Derived operator/ (const Derived& d) const;

  /**
   * Exact division assignment.
   *
   * @param d: the divisor.
   * @return a reference to this after assignment.
   */
  inline virtual Derived& operator/= (const Derived& d);

  /**
   * Get the remainder of *this and b;
   * @param b: the divisor
   * @return the remainder
   */
  inline virtual Derived operator%(const Derived& b) const;

  /**
   * Assign *this to be the remainder of *this and b.
   * @param b: the divisor
   * @return this after assignment.
   */
  inline virtual Derived& operator%=(const Derived& b);
};

#include "EuclideanDomainElement.inl"

#endif // __EUCLIDEAN_DOMAIN_ELEMENT_H_
