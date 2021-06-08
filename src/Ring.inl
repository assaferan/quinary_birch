// implementation of Ring methods

#include <sstream>

template <class Derived>
Derived Ring<Derived>::operator+ (const Derived& other) const
{
  Derived sum = *(this->getPtr());
  sum += other;
  return sum;
}

template <class Derived>
Derived Ring<Derived>::operator- (const Derived& other) const
{
  Derived diff = *(this->getPtr());
  diff -= other;
  return diff;
}

template <class Derived>
Derived Ring<Derived>::operator* (const Derived& other) const
{
  Derived prod = *(this->getPtr());
  prod *= other;
  return prod;
}

template <class Derived>
Derived Ring<Derived>::operator- () const {
  Derived neg = *(this->getPtr());
  return neg.zero() - (*this->getPtr());
}

template <class Derived>
Derived& Ring<Derived>::operator^= (unsigned long long int e)
{
  if (e == 0) return this->one();
  if (e == 1) return *(this->getPtr());

  Derived a = *(this->getPtr());;

  // !! - TODO - eliminate recursion here
  (*this) ^= (e>>1);
  (*this) *= (*(this->getPtr()));

  if (e%2) (*this) *= a;

  return (*(this->getPtr()));
}

template <class Derived>
Derived Ring<Derived>::operator^ (unsigned long long int e) const
{
  Derived pow = *(this->getPtr());
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
