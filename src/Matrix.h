#ifndef __MATRIX_H_
#define __MATRIX_H_

#include <cassert>
#include <type_traits>
#include <vector>

#include "birch.h"
#include "UnivariatePolyInt.h"

// R should be a ring element,
// Parent the parent ring class

// This is a general matrix class (matrices of mxn)

template<class R, class Parent>
class Matrix
{
  static_assert(std::is_base_of<RingElement<R,Parent>,R>::value);
  static_assert(std::is_base_of<Ring<Parent,R>,Parent>::value);
  
public:
  Matrix(const std::vector<R> & data, size_t nrows, size_t ncols)
    : _nrows(nrows), _ncols(ncols), _data(data)
  {
    assert(data.size() == nrows*ncols);
    if (data.size() > 0) {
      _base = data[0].parent();
    }
  }
  
  template <size_t n>
  Matrix(const R data[n][n]);
  
  template <size_t n>
  Matrix(const SquareMatrix<R,Parent,n> &);
  
  Matrix(std::shared_ptr<const Parent> base_ring, size_t nrows, size_t ncols)
    : _nrows(nrows),
      _ncols(ncols),
      _data(nrows*ncols, base_ring->zero()),
      _base(base_ring)
  {}
  
  inline const R & operator()(size_t row, size_t col) const
  {
    assert( (_ncols*row+col < _nrows*_ncols) && (0 <= _ncols*row+col) );
    return _data[_ncols*row+col];
  }
  
  inline R & operator()(size_t row, size_t col)
  {
    assert( (_ncols*row+col < _nrows*_ncols) && (0 <= _ncols*row+col) );
    return _data[_ncols*row+col];
  }
  
  inline size_t nrows(void) const {return _nrows;}
  inline size_t ncols(void) const {return _ncols;}

  // return the i-th row
  std::vector<R> operator[](size_t i) const;
  
  R determinant(void) const;

  size_t rank(void) const;

  Matrix<R,Parent> kernel(void) const;

  Matrix<R,Parent> leftKernel(void) const;

  // restrict matrix to the subspace specified by the argument
  Matrix<R,Parent> restrict(const Matrix<R,Parent> &) const;

  R trace(void) const;
  
  UnivariatePolyInt<Z> charPoly(void) const;
  
  static Matrix<R,Parent> diagonalJoin(const std::vector< Matrix<R,Parent> > &);

  static Matrix<R,Parent> identity(std::shared_ptr< const Parent >, size_t);
  
  // TODO - just change access resolution to the same vector instead
  Matrix<R,Parent> transpose() const;

  // arithmetic
  Matrix<R,Parent> operator*(const Matrix<R,Parent> &) const;

  Matrix<R,Parent>& operator+=(const Matrix<R,Parent> &);
  Matrix<R,Parent>& operator-=(const Matrix<R,Parent> &);
  Matrix<R,Parent>& operator*=(const Matrix<R,Parent> &);

  Matrix<R,Parent> operator*(const R & a) const;

  // algorithms
  
  void swapRows(size_t, size_t);
  
  // in-place row-echelon form for the matrix echelon,
  // returns the rank and the transformation matrix trans
  static size_t rowEchelon(Matrix<R,Parent> & , Matrix<R,Parent>& );

  // printing
  
  inline void print(std::ostream & os) const
  {
    for (size_t i = 0; i < _nrows; i++) {
      for (size_t j = 0; j < _ncols; j++)
	os << (*this)(i,j) << " ";
      os << std::endl;
    }
    return;
  }

  bool isZero(void) const;
  bool isOne(void) const;

  Matrix<R,Parent>& makeZero(void);
  
protected:
  size_t _nrows;
  size_t _ncols;
  std::vector<R> _data;

  std::shared_ptr< const Parent > _base;
};

template<class R, class Parent>
inline Matrix<R, Parent> operator*(const R & a, const Matrix<R,Parent> & mat)
{ return mat*a; }

#include "Matrix.inl"

#endif // __MATRIX_H_
