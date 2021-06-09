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

template<class Derived>
Derived Field<Derived>::operator^ (long long int e) const
{
  unsigned long long int abs_e = e < 0 ? -e : e;
  Derived ret = this->operator^(abs_e);
  return e < 0 ? ret.inverse() : ret;
}

template<class Derived>
Derived& Field<Derived>::operator^= (long long int e)
{
  unsigned long long int abs_e = e < 0 ? -e : e;
  this->operator^=(abs_e);
  if (e < 0)
    *(this->getPtr()) = this->inverse();
  return (*this->getPtr());
}
