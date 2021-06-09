#ifndef __FIELD_ELEMENT_H
#define __FIELD_ELEMENT_H

#include "EuclideanDomainElement.h"

template<class Derived, class DerivedParent>
class FieldElement : public virtual EuclideanDomainElement<Derived, DerivedParent>
{
  public:

  /**
   * Get the inverse of *this.
   *
   * @return the inverse
   */
  virtual Derived inverse() const = 0;

  inline virtual typename EuclideanDomainElement<Derived,DerivedParent>::DivRes
  euclideanDivision(const Derived & b) const override;

  inline virtual Derived gcd(const Derived &) const override;

  inline virtual Derived operator^ (long long int e) const;

  inline virtual Derived& operator^= (long long int e);
};

#include "FieldElement.inl"

#endif // __FIELD_ELEMENT_H
