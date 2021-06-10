#ifndef __ISOMETRY_H_
#define __ISOMETRY_H_

#include "Integer.h"
#include "QuadForm.h"
#include "SquareMatrix.h"
#include "Vector.h"

// !! TODO - think if we really need a different class for Isometry

// we only use isometries with integral types

template<typename R, size_t n>
using Mat = SquareMatrix<Integer<R>, IntegerRing<R>, n>;

template<typename R, size_t n>
using Vec = Vector<Integer<R>, IntegerRing<R>, n>;

template<typename R, size_t n>
using QF = QuadForm<Integer<R>, IntegerRing<R>, n>;

template<typename R, size_t n>
class Isometry
{
public:
  // c-tors
  Isometry() :
    _a(Mat<R,n>::identity(IntegerRing<R>::getInstance().getPtr())),
    _scale(Integer<R>::one())
  {}

  Isometry(const Mat<R,n> & mat) : _a(mat), _scale(Integer<R>::one()) {}

  Isometry(const Mat<R,n> & mat, const Integer<R> & scale) :
    _a(mat), _scale(scale) { this->rescale(); }

  // access - set/get
  inline const Integer<R> & getScale(void) const
  { return this->_scale; }
  
  inline void setValues(const Mat<R,n> & mat)
  { this->_a = mat; }

  inline void setIdentity(void)
  { this->_a.set_identity(); this->_scale = Integer<R>::one(); }

  inline void setScale(const Integer<R> & scale)
  { this->_scale = scale; }

  inline void rescale(void);

  inline const Integer<R> & operator()(size_t i, size_t j) const
  { return this->a(i, j); }

  inline Integer<R> & operator()(size_t i, size_t j)
  { return this->a(i, j); }

  // basic operations
  
  inline Isometry<R,n> inverse(void) const;

  inline Isometry<R,n> transpose(void) const
  { return Isometry(this->_a.transpose(), this->_scale); }

  // arithmetic
  
  inline Isometry<R, n> operator*(const Isometry<R,n>& s) const
  { return Isometry((this->_a)*s._a, (this->_scale)*s._scale); }

  Vec<R,n> operator*(const Vec<R,n>& vec) const
  { return (this->_a)*vec; }

  // assignment
  inline Isometry<R,n> & operator=(const Isometry<R,n>& other)
  { if (this != &other) {
      this->_a = other._a;
      this->_scale = other._scale;
    }
    return *this;
  }
  
  inline Mat<R,n> transform(const Mat<R,n>& from) const;

  // we save some clocks by returning once a single coordinate is mismatched.
  inline bool isIsometry(const QF<R,n>& from, const QF<R,n>& to) const;
  
  inline void updatePerm(const Vec<size_t,n> & perm);

  friend std::ostream& operator<<(std::ostream& os, const Isometry<R,n>& s)
  { os << s._a; return os; }

  inline bool operator==(const Isometry<R,n> & other) const
  {return (this->_a == other._a);}

  inline bool operator<(const Isometry<R,n> & other) const
  {return (this->_a < other._a);}

  Mat<R,n> _a;
  
protected:
 
  Integer<R> _scale;
};

#include "Isometry.inl"

#endif // __ISOMETRY_H_
