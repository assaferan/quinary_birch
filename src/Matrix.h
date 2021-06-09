#ifndef __MATRIX_H_
#define __MATRIX_H_

#include <cassert>
#include <type_traits>
#include <vector>

#include "RingElement.h"

// #include "Polynomial.h"

// R should be a ring element,
// Parent the parent ring class

// This is a general matrix class (matrices of mxn)

template<class R, class Parent>
class Matrix
{
  static_assert(std::is_base_of<RingElement<R, Parent>, R>::value);
  static_assert(std::is_base_of<Ring<Parent, R>, Parent>::value);
  
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
  
  /*
  template <size_t n>
  Matrix(const SquareMatrix<R, n> &);
  */
  
  Matrix(std::shared_ptr<const Parent> base_ring, size_t nrows, size_t ncols)
    : _nrows(nrows), _ncols(ncols), _data(nrows*ncols, base_ring->zero())),
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
  
  inline size_t nrows() const {return _nrows;}
  inline size_t ncols() const {return _ncols;}

  // return the i-th row
  inline std::vector<R> operator[](size_t i) const;
  
  inline R determinant() const;

  inline size_t rank() const;

  inline Matrix<R,Parent> kernel() const;

  inline Matrix<R,Parent> leftKernel() const;

  // restrict matrix to the subspace specified by the argument
  inline Matrix<R,Parent> restrict(const Matrix<R,Parent> & ) const;

  inline R trace() const;
  
  // UnivariatePoly<Z> char_poly() const;
  
  inline static Matrix<R,Parent> diagonalJoin(const std::vector< Matrix<R,Parent> > &);

  inline static Matrix<R,Parent> identity(size_t);
  
  // TODO - just change access resolution to the same vector instead
  inline Matrix<R,Parent> transpose() const;

  // arithmetic
  inline Matrix<R,Parent> operator*(const Matrix<R,Parent> &) const override;

  inline Matrix<R,Parent>& operator+=(const Matrix<R,Parent> &) override;
  inline Matrix<R,Parent>& operator-=(const Matrix<R,Parent> &) override;
  inline Matrix<R,Parent>& operator*=(const Matrix<R,Parent> &) override;

  inline Matrix<R,Parent> operator*(const R & a) const;

  // algorithms
  
  inline void swapRows(size_t, size_t);
  
  // in-place row-echelon form for the matrix echelon,
  // returns the rank and the transformation matrix trans
  inline static size_t rowEchelon(Matrix<R,Parent> & , Matrix<R,Parent>& );

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

  inline Matrix<R,Parent>* getPtr() override {return this;}

  inline const Matrix<R,Parent>* getPtr() const override {return this;}

  inline bool isZero() const override;
  inline bool isOne() const override;

  inline Matrix<R,Parent>& makeZero() override;
  inline Matrix<R,Parent>& makeOne() override;
  
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
