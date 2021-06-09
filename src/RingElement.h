#ifndef __RING_ELEMENT_H_
#define __RING_ELEMENT_H_

#include <ostream>
#include <string>

#include "Ring.h"

// We follow the same template pattern as in BPAS
// Derived is a concrete class derived from RingElement

// on the other hand, we want to be able to implement as little
// as possible in the derived class,
// so we define as many methods as possible already here
// i.e. RingElement is not a completely virtual ADT

// Forward declaration of friend operators

template<class Derived>
class RingElement;

template<class Derived>
std::ostream& operator<< (std::ostream& ostream,
			  const RingElement<Derived> & d);

template<class Derived>
std::ostream& operator<< (std::ostream& ostream,
			  RingElement<Derived> && d);

template <class Derived>
class RingElement
{
public:

  virtual Derived* getPtr() = 0;

  virtual const Derived* getPtr() const = 0;
  
  /**
   * Determine if *this ring element is zero, that is the additive identity.
   *
   * returns true iff *this is zero.
   */
   
  virtual bool isZero() const = 0;
  
  /**
   * Make *this ring element zero.
   */
  virtual Derived& makeZero() = 0;

  /**
   * Determine if *this ring element is one, that is the multiplication identity.
   *
   * returns true iff *this is one.
   */
  virtual bool isOne() const = 0;

  /**
   * Make *this ring element one.
   */
  virtual Derived& makeOne() = 0;

  /**
   * Copy assignment.
   */
  virtual Derived& operator= (const Derived&) = 0;

  /**
   * Addition.
   */
  inline virtual Derived operator+ (const Derived&) const;

  /**
   * Addition assignment.
   */
  virtual Derived& operator+= (const Derived&) = 0;

  /**
   * Subtraction.
   */
  inline virtual Derived operator- (const Derived&) const;

  /**
   * Subtraction assignment.
   */
  virtual Derived& operator-= (const Derived&) = 0;

  /**
   * Negation.
   */
  inline virtual Derived operator- () const;

  /**
   * Multiplication.
   */
  inline virtual Derived operator* (const Derived&) const;

  /**
   * Multiplication assignment.
   */
  virtual Derived& operator*= (const Derived&) = 0;

  /**
   * Exponentiation.
   */
  inline virtual Derived operator^ (unsigned long long int e) const;

  /**
   * Exponentiation assignment.
   */
  inline virtual Derived& operator^= (unsigned long long int e);

  /**
   * Equality test,
   *
   * returns true iff equal
   */
  virtual bool operator== (const Derived&) const = 0;

  /**
   * Inequality test,
   *
   * returns true iff not equal.
   */
  inline virtual bool operator!= (const Derived& other) const
  { return !((*this) == other); }

  /**
   * Print the Ring element.
   */
  virtual void print(std::ostream&) const = 0;

  /**
   * Convert the Ring element to a string.
   *
   * Simple delegation of printing to a stringstream to obtain a string.
   * Overriding the print method is sufficient for sub-classes
   * to make use of this method.
   *
   * returns the string representation of the Ring element.
   */
  inline virtual std::string toString() const;

  /**
   * Output operator.
   *
   * Defines a to string conversion.
   */
  friend std::ostream& operator<< <Derived>(std::ostream& ostream,
					    const RingElement<Derived>& d);

  friend std::ostream& operator<< <Derived>(std::ostream& ostream,
					    RingElement<Derived>&& d);

  virtual std::shared_ptr<const Ring<Derived> > parent() const = 0;

};

#include "RingElement.inl"

#endif // __RING_ELEMENT_H_
