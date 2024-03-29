#ifndef __VECTOR_H_
#define __VECTOR_H_

#include "birch.h"

template<class R, class Parent, size_t n>
class Vector {
  
  static_assert(std::is_base_of<RingElement<R,Parent>, R>::value, "R template parameter has to inherit from RingElement.");
  static_assert(std::is_base_of<Ring<Parent,R>, Parent>::value, "Parent template parameter has to inherit from Ring.");
  
public:

  // c-tor - default constructor constructs the zero vector
  Vector(std::shared_ptr<const Parent> base_ring) : _base(base_ring)
  {
    for (size_t i = 0; i < n; i++) this->_v[i] = base_ring->zero();
  }
  // access
  inline const R& operator[](size_t i) const {return this->_v[i]; }
  inline R& operator[](size_t i) {return this->_v[i];}
  
  inline std::shared_ptr<const Parent> baseRing(void) const { return _base; };

  // arithmetic
  inline Vector<R,Parent,n> operator-(void) const;
  inline Vector<R,Parent,n> operator+(const Vector<R,Parent,n> &) const;
  inline Vector<R,Parent,n> operator-(const Vector<R,Parent,n> &) const;
  inline Vector<R,Parent,n> operator*(const R & a) const;

  inline Vector<R,Parent,n> & operator+=(const Vector<R,Parent,n> &);
  inline Vector<R,Parent,n> & operator-=(const Vector<R,Parent,n> &);
  
  // considering the vector as a row vector
  inline Vector<R,Parent,n> operator*(const SquareMatrix<R,Parent,n>& mat) const;

  // inner product
  inline static R innerProduct(const Vector<R,Parent,n> & , const Vector<R,Parent,n> &);

  // assignment
  inline Vector<R,Parent,n> & operator=(const Vector<R,Parent,n>& other)
  {
    if (this != (&other)) {
      for (size_t i = 0; i < n; i++) this->_v[i] = other[i];
    }
    return (*this);
  }
  
  // booleans
  inline bool operator==(const Vector<R,Parent,n> &) const;

  inline bool operator!=(const Vector<R,Parent,n> & other) const
  {return !((*this)==other);}

  inline bool operator<(const Vector<R,Parent,n> &) const;

  inline bool isZero(void) const
  {
    for (size_t i = 0; i < n; i++)
      if (!this->_v[i].isZero())
	return false;
    return true;
  }

  inline std::ostream & prettyPrint(std::ostream &, size_t upTo = n) const;
  
protected:
  std::shared_ptr<const Parent> _base;
  R _v[n];
  
};

//scalar multiple
template<class R, class Parent, size_t n>
inline Vector<R,Parent,n> operator*(const R & a, const Vector<R,Parent,n> & v)
{ return v*a; }

// we put it outsisde the class to avoid partial specialization
// !! TODO - figure out a better way

template<typename R, typename S, typename T, size_t n>
std::shared_ptr< VectorFp<S,T,n> >
mod(const VectorInt<R,n> & a, std::shared_ptr< Fp<S,T> > GF);

namespace std
{
  template<class R, class Parent, size_t n>
  struct hash<Vector<R,Parent,n> >
  {
    Z64 operator()(const Vector<R,Parent,n>& vec) const
    {
      Z64 fnv = FNV_OFFSET;
      for (size_t i = 0; i < n; i++)
	fnv = (fnv ^ std::hash<R>{}(vec[i])) * FNV_PRIME;
            
      return fnv;
    }
  };
}

#include "Vector.inl"

#endif // __VECTOR_H_
