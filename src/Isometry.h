#ifndef __ISOMETRY_H_
#define __ISOMETRY_H_

#include "birch.h"
#include "Integer.h"
#include "QuadForm.h"
#include "SquareMatrix.h"
#include "Vector.h"

// !! TODO - think if we really need a different class for Isometry

// we only use isometries with integral types

template<typename R, size_t n>
class Isometry
{
public:
  // c-tors
  Isometry() :
    _a(SquareMatrixInt<R,n>::identity(std::make_shared< IntegerRing<R> >())),
    _scale(Integer<R>::one())
  {}

  Isometry(const SquareMatrixInt<R,n> & mat) : _a(mat), _scale(Integer<R>::one()) {}

  Isometry(const SquareMatrixInt<R,n> & mat, const Integer<R> & scale) :
    _a(mat), _scale(scale) { this->rescale(); }

  // access - set/get
  inline const Integer<R> & getScale(void) const
  { return this->_scale; }
  
  inline void setValues(const SquareMatrixInt<R,n> & mat)
  { this->_a = mat; }

  inline void setIdentity(void)
  { this->_a.setIdentity(); this->_scale = Integer<R>::one(); }

  inline void setScale(const Integer<R> & scale)
  { this->_scale = scale; }

  void rescale(void);

  inline const Integer<R> & operator()(size_t i, size_t j) const
  { return this->a(i, j); }

  inline Integer<R> & operator()(size_t i, size_t j)
  { return this->a(i, j); }

  // basic operations
  
  Isometry<R,n> inverse(void) const;

  inline Isometry<R,n> transpose(void) const
  { return Isometry(this->_a.transpose(), this->_scale); }

  // arithmetic
  
  inline Isometry<R,n> operator*(const Isometry<R,n>& s) const
  { return Isometry((this->_a)*s._a, (this->_scale)*s._scale); }

  inline VectorInt<R,n> operator*(const VectorInt<R,n>& vec) const
  { return (this->_a)*vec; }

  // assignment
  inline Isometry<R,n> & operator=(const Isometry<R,n>& other)
  { if (this != &other) {
      this->_a = other._a;
      this->_scale = other._scale;
    }
    return *this;
  }
  
  SquareMatrixInt<R,n> transform(const SquareMatrixInt<R,n>& from) const;

  // we save some clocks by returning once a single coordinate is mismatched.
  bool isIsometry(const QuadFormInt<R,n>& from, const QuadFormInt<R,n>& to) const;
  
  void updatePerm(const VectorInt<size_t,n> & perm);

  inline friend std::ostream& operator<<(std::ostream& os, const Isometry<R,n>& s)
  { os << s._a; return os; }

  inline bool operator==(const Isometry<R,n> & other) const
  {return (this->_a == other._a);}

  inline bool operator<(const Isometry<R,n> & other) const
  {return (this->_a < other._a);}

  SquareMatrixInt<R,n> _a;
  
protected:
 
  Integer<R> _scale;
};

#include "Isometry.inl"

#endif // __ISOMETRY_H_
