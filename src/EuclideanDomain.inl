#include <cassert>

template<class Derived>
typename EuclideanDomain<Derived>::XGcdRes
EuclideanDomain<Derived>::xgcd(const Derived& b) const
{
  Derived a = *(this->getPtr());
  Derived old_r = a;
  Derived r = b;
  // we don't know if there is a default constructor
  // so we use the copy constructor instead.
  Derived s = r;
  Derived t = r;
  Derived old_s = r;
  Derived old_t = r;
  Derived temp = r;
  
  old_s.makeOne();
  old_t.makeZero();
  s.makeZero();
  t.makeOne();
  
  while (!r.isZero()) {
    typename EuclideanDomain<Derived>::DivRes q_r = old_r.euclideanDivision(r);
    old_r = r;
    r = q_r.second;
    
    temp = t;
    t = old_t - q_r.first * t;
    old_t = temp;

    temp = s;
    s = old_s - q_r.first * s;
    old_s = temp;

    assert(old_r == old_s*a+old_t*b);
    assert(r == s*a+t*b);
  }
  
  return std::make_tuple(old_r,old_s,old_t);
}

// This could use xgcd, but this implementation is lightly quicker
template<class Derived>
Derived EuclideanDomain<Derived>::gcd(const Derived& b) const
{
  Derived old_r = *(this->getPtr());
  Derived r = b;
  while (!r.isZero()) {
#ifdef DEBUG_LEVEL_FULL
    std::cerr << "r = " << r << std::endl;
    std::cerr << "old_r = " << old_r << std::endl;
#endif
    typename EuclideanDomain<Derived>::DivRes q_r = old_r.euclideanDivision(r);
    old_r = r;
    r = q_r.second;
  }
  
  return old_r;
}
  
/**
 * Get the quotient of *this and b.
 *
 * @param b: the divisor
 * @return the quotient
 */
template<class Derived>
Derived EuclideanDomain<Derived>::quotient(const Derived& b) const
{
  typename EuclideanDomain<Derived>::DivRes q_r = this->euclideanDivision(b);
  return q_r.first;
}

/**
 * Get the remainder of *this and b.
 * @param b: the divisor
 * @return the remainder
 */
template<class Derived>
Derived EuclideanDomain<Derived>::remainder(const Derived& b) const
{
  typename EuclideanDomain<Derived>::DivRes q_r = this->euclideanDivision(b);
  return q_r.second;
}

/**
 * Exact division.
 *
 * @param d: the divisor.
 * @return the equotient.
 */
template<class Derived>
Derived EuclideanDomain<Derived>::operator/ (const Derived& d) const
{
  return this->quotient(d);
}

/**
 * Exact division assignment.
 *
 * @param d: the divisor.
 * @return a reference to this after assignment.
 */
template<class Derived>
Derived& EuclideanDomain<Derived>::operator/= (const Derived& d)
{
 ((*this) = (*this) / d);
 return *(this->getPtr());
}

/**
 * Get the remainder of *this and b;
 * @param b: the divisor
 * @return the remainder
 */
template<class Derived>
Derived EuclideanDomain<Derived>::operator%(const Derived& b) const
{
  return this->remainder(b);
}

/**
 * Assign *this to be the remainder of *this and b.
 * @param b: the divisor
 * @return this after assignment.
 */
template<class Derived>
Derived& EuclideanDomain<Derived>::operator%=(const Derived& b)
{
  ((*this) = (*this) % b);
  return *(this->getPtr());
}
