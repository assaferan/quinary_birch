#ifndef __RING_H_
#define __RING_H_

#include <memory>

// This is an ADT for a ring - basically a factory for proucing ring elements
// which should be a parent of a ring element
// The DerivedElement class is a class derived from RingElement

// Pondering - we could decide that all the actions belong to this
// class (as Jeff did with Fp), and then let the element wrap it ?

template <class Derived, class DerivedElement>
class Ring : public virtual std::enable_shared_from_this<const Derived>
{
public:

  virtual std::shared_ptr<const Derived> getPtr() const = 0;
  
  // producing the global constants of the ring
  virtual DerivedElement zero() const = 0;
  virtual DerivedElement one() const = 0;
};

#endif // __RING_H_
