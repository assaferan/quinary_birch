#ifndef __VECTOR_INT_H_
#define __VECTOR_INT_H_

#include "birch.h"

template<typename R, size_t n>
class VectorInt
{  
public:

  // c-tor - default constructor constructs the zero vector
  VectorInt()
  {
    for (size_t i = 0; i < n; i++) this->_v[i] = 0;
  }
  // access
  inline const R& operator[](size_t i) const {return this->_v[i]; }
  inline R& operator[](size_t i) {return this->_v[i];}

  // arithmetic
  inline VectorInt<R,n> operator-(void) const;
  inline VectorInt<R,n> operator+(const VectorInt<R,n> &) const;
  inline VectorInt<R,n> operator-(const VectorInt<R,n> &) const;
  inline VectorInt<R,n> operator*(const R & a) const;

  inline VectorInt<R,n> & operator+=(const VectorInt<R,n> &);
  inline VectorInt<R,n> & operator-=(const VectorInt<R,n> &);
  
  // considering the vector as a row vector
  inline VectorInt<R,n> operator*(const SquareMatrixInt<R,n>& mat) const;

  // inner product
  inline static R innerProduct(const VectorInt<R,n> & , const VectorInt<R,n> &);

  // assignment
  inline VectorInt<R,n> & operator=(const VectorInt<R,n>& other)
  {
    if (this != (&other)) {
      for (size_t i = 0; i < n; i++) this->_v[i] = other[i];
    }
    return (*this);
  }
  
  // booleans
  inline bool operator==(const VectorInt<R,n> &) const;

  inline bool operator!=(const VectorInt<R,n> & other) const
  {return !((*this)==other);}

  inline bool operator<(const VectorInt<R,n> &) const;

  inline bool isZero(void) const
  {
    for (size_t i = 0; i < n; i++)
      if (!this->_v[i].isZero())
	return false;
    return true;
  }

  inline std::ostream & prettyPrint(std::ostream &, size_t upTo = n) const;
  
protected:
  R _v[n];
  
};

//scalar multiple
template<typename R, size_t n>
inline VectorInt<R,n> operator*(const R & a, const VectorInt<R,n> & v)
{ return v*a; }

namespace std
{
  template<typename R, size_t n>
  struct hash<VectorInt<R,n> >
  {
    Z64 operator()(const VectorInt<R,n>& vec) const
    {
      Z64 fnv = FNV_OFFSET;
      for (size_t i = 0; i < n; i++)
	fnv = (fnv ^ vec[i]) * FNV_PRIME;
            
      return fnv;
    }
  };
}

#include "VectorInt.inl"

#endif // __VECTOR_INT_H_
