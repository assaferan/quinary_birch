#ifndef __FIELD_H
#define __FIELD_H

#include "EuclideanDomain.h"

template<class Derived>
class Field : public virtual EuclideanDomain<Derived>
{
  public:

  /**
   * Get the inverse of *this.
   *
   * @return the inverse
   */
  virtual Derived inverse() const = 0;

  virtual typename EuclideanDomain<Derived>::DivRes
  euclideanDivision(const Derived & b) const override;

  virtual Derived gcd(const Derived &) const override;

  virtual Derived operator^ (long long int e) const;

  virtual Derived& operator^= (long long int e);
};

#include "Field.inl"

#endif // __FIELD_H
