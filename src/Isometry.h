#ifndef __ISOMETRY_H_
#define __ISOMETRY_H_

#include "birch.h"
#include "Integer.h"
#include "QuadForm.h"
#include "SquareMatrixInt.h"
#include "Vector.h"

// !! TODO - think if we really need a different class for Isometry

// we only use isometries with integral types

template<typename R, size_t n>
class Isometry
{
public:
  // c-tors
  Isometry() :
    _a(SquareMatrixInt<R,n>::identity()),
    _scale(1)
  {}

  Isometry(const SquareMatrixInt<R,n> & mat) : _a(mat), _scale(1) {}

  Isometry(const SquareMatrixInt<R,n> & mat, const R & scale) :
    _a(mat), _scale(scale) { this->rescale(); }

  // access - set/get
  inline const SquareMatrixInt<R,n> & integralMatrix(void) const
  {return this->_a;}
  
  inline const R & getScale(void) const
  { return this->_scale; }
  
  inline void setValues(const SquareMatrixInt<R,n> & mat)
  { this->_a = mat; }

  inline void setIdentity(void)
  { this->_a.setIdentity(); this->_scale = 1; }

  inline void setScale(const R & scale)
  { this->_scale = scale; }

  void rescale(void);

  inline const R & operator()(size_t i, size_t j) const
  { return this->_a(i,j); }

  inline R & operator()(size_t i, size_t j)
  {return this->_a(i,j); }

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
  { os << s._a << " scaled by " << s._scale; return os; }

  inline bool operator==(const Isometry<R,n> & other) const
  {return (other._scale * this->_a == this->_scale * other._a);}

  inline bool operator!=(const Isometry<R,n> & other) const
  {return !((*this) == other); }   
  /*
  inline bool operator<(const Isometry<R,n> & other) const
  {return (other._scale * this->_a < this->_scale * other._a);}
  */
  inline Rational<R> determinant() const
  {return Rational<R>(_a.determinant(), Integer<R>(_scale)^n); }

  inline bool isOne(void) const
  {return _a == _scale * SquareMatrixInt<R,n>::identity(); }

  inline bool isIdentity(void) const
  {return this->isOne();}

  inline SquareMatrixInt<R,n> hermiteForm(const R & d) const
  {return _a.hermiteForm(d); }

  inline VectorRat<R,n> operator*(const VectorRat<R,n> & v) const
  {
    std::shared_ptr< const RationalField<R> > QQ = std::make_shared< const RationalField<R> >();
    SquareMatrixRat<R,n> mat(QQ);
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < n; j++) {
	Rational<R> elt(_a(i,j), _scale);
	mat(i,j) = elt;
      }
    return mat * v;
  }
  
protected:
  SquareMatrixInt<R,n> _a;
  R _scale;
};

namespace std {

  template <typename R, size_t n>
  struct hash< Isometry<R,n> >
  {
    Z64 operator()(const Isometry<R,n> & s) const
    {
      Z64 fnv = FNV_OFFSET;
      for (size_t i = 0; i < n; i++)
	for (size_t j = 0; j < n; j++) {
	  Rational<R> elt(s(i,j), s.getScale());
	  fnv = (fnv ^ std::hash< Rational<R> >{}(elt)) * FNV_PRIME;
	}
      return fnv;
    }
  };
}

#include "Isometry.inl"

#endif // __ISOMETRY_H_
