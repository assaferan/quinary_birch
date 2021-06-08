#ifndef __RING_H
#define __RING_H

class RingElement;

class Ring {
 public:
  virtual RingElement zero() const = 0;
  virtual RingElement one() const = 0;

  virtual ~Ring() {}
};

class RingElement {
 public:

  RingElement(const Ring & R) : parent_(R) {}

  const Ring & parent(void) const
  {return this->parent_; }
  
  // assignment operator
  virtual RingElement & operator=(const RingElement &) = 0;
  
  // arithmetic
  virtual RingElement operator-() const = 0;
  virtual RingElement operator+(const RingElement & ) const = 0;
  virtual RingElement operator-(const RingElement & ) const = 0;
  virtual RingElement operator*(const RingElement & ) const = 0;

  virtual RingElement & operator+=(const RingElement & ) = 0;
  virtual RingElement & operator-=(const RingElement & ) = 0;
  virtual RingElement & operator*=(const RingElement & ) = 0;

  virtual ~RingElement() {}
  
 protected:
  const Ring & parent_;

};

#endif // __RING_H
