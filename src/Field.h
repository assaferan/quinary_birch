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

  typename EuclideanDomain<Derived>::DivRes
  euclideanDivision(const Derived & b) const;
  
};

#include "Field.inl"

#endif // __FIELD_H
