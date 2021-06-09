#include <cassert>

template<class Derived>
typename EuclideanDomain<Derived>::DivRes
Field<Derived>::euclideanDivision(const Derived & b) const
{
  assert(!b.isZero());
  Derived a = *(this->getPtr());
  Derived b_inv = b.inverse();
  Derived z = a;
  z.makeZero();
  
  return std::make_pair(a*b_inv, z);
}

template<class Derived>
Derived Field<Derived>::gcd(const Derived & b) const
{
  if (b.isZero() && this->isZero()) {
    Derived z = b;
    z.makeZero();
    return z;
  }
  Derived id = b;
  id.makeOne();
  
  return id;
}
