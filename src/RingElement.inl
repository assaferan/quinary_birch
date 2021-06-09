// implementation of Ring methods

#include <sstream>

template <class Derived, class DerivedParent>
Derived RingElement<Derived,DerivedParent>::operator+ (const Derived& other) const
{
  Derived sum = *(this->getPtr());
  sum += other;
  return sum;
}

template <class Derived, class DerivedParent>
Derived RingElement<Derived,DerivedParent>::operator- (const Derived& other) const
{
  Derived diff = *(this->getPtr());
  diff -= other;
  return diff;
}

template <class Derived, class DerivedParent>
Derived RingElement<Derived,DerivedParent>::operator* (const Derived& other) const
{
  Derived prod = *(this->getPtr());
  prod *= other;
  return prod;
}

template <class Derived, class DerivedParent>
Derived RingElement<Derived,DerivedParent>::operator- () const {
  Derived neg = *(this->getPtr());
  return neg.makeZero() - (*this->getPtr());
}

template <class Derived, class DerivedParent>
Derived& RingElement<Derived,DerivedParent>::operator^= (unsigned long long int exp)
{
  if (exp == 1) return *(this->getPtr());
  
  this->makeOne();
  Derived base = *(this->getPtr());;
  unsigned long long int e = exp;
  
  while (e) {
    if (e & 1)
      (*this) *= base;
    e >>= 1;
    base *= base;
  }

  return (*(this->getPtr()));
}

template <class Derived, class DerivedParent>
Derived RingElement<Derived,DerivedParent>::operator^ (unsigned long long int e) const
{
  Derived pow = *(this->getPtr());
  pow ^= e;
  return pow;
}

template <class Derived, class DerivedParent>
std::string RingElement<Derived,DerivedParent>::toString() const {
  std::stringstream ss;
  print(ss);
  return ss.str();
}

/**
 * Output operator.
 *
 * Defines a to string conversion.
 */
template <class Derived, class DerivedParent>
std::ostream& operator<< (std::ostream& ostream,
			  const RingElement<Derived,DerivedParent>& d) {
  d.print(ostream);
  return ostream;
}

template <class Derived, class DerivedParent>
std::ostream& operator<< (std::ostream& ostream, RingElement<Derived,DerivedParent>&& d) {
  d.print(ostream);
  return ostream;
}
