// implementation of Ring methods

#include <sstream>

template <class Derived>
Derived Ring<Derived>::operator+ (const Derived& other) const
{
  Derived sum = (*this);
  sum += other;
  return sum;
}

template <class Derived>
Derived Ring<Derived>::operator- (const Derived& other) const
{
  Derived diff = (*this);
  diff -= other;
  return diff;
}

template <class Derived>
Derived Ring<Derived>::operator* (const Derived& other) const
{
  Derived prod = (*dynamic_cast<const Derived*>(this));
  prod *= other;
  return prod;
}

template <class Derived>
Derived Ring<Derived>::operator- () const {
  Derived neg = (*this);
  return neg.zero() - (*this);
}

template <class Derived>
Derived& Ring<Derived>::operator^= (unsigned long long int e)
{
  if (e == 0) return this->one();
  if (e == 1) return (*this);

  Derived a = (*this);

  // !! - TODO - eliminate recursion here
  (*this) ^= (e>>1);
  (*this) *= (*this);

  if (e%2) (*this) *= a;

  return (*this);
}

template <class Derived>
Derived Ring<Derived>::operator^ (unsigned long long int e) const
{
  Derived pow = (*this);
  pow ^= e;
  return pow;
}

template <class Derived>
std::string Ring<Derived>::toString() const {
  std::stringstream ss;
  print(ss);
  return ss.str();
}

/**
 * Output operator.
 *
 * Defines a to string conversion.
 */
template <class Derived>
std::ostream& operator<< (std::ostream& ostream, const Ring<Derived>& d) {
  d.print(ostream);
  return ostream;
}

template <class Derived>
std::ostream& operator<< (std::ostream& ostream, Ring<Derived>&& d) {
  d.print(ostream);
  return ostream;
}
