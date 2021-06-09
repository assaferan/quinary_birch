#ifndef __MATRIX_H_
#define __MATRIX_H_

#include <cassert>
#include <type_traits>
#include <vector>

// #include "Polynomial.h"

// R should be a ring, and the matrix class forms a ring by itself

template<class R>
class Matrix : public virtual Ring< Matrix<R> >
{
  static_assert(std::is_base_of<Ring, R>::value);
  
public:
  Matrix(const std::vector<R> & data, size_t nrows, size_t ncols)
    : _nrows(nrows), _ncols(ncols), _data(data) {}
  
  template <size_t n>
  Matrix(const R data[n][n]);
  
  /*
  template <size_t n>
  Matrix(const SquareMatrix<R, n> &);
  */
  
  Matrix(size_t nrows, size_t ncols)
    : _nrows(nrows), _ncols(ncols), _data(nrows*ncols) {}
  
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

  inline Matrix<R> kernel() const;

  inline Matrix<R> leftKernel() const;

  // restrict matrix to the subspace specified by the argument
  inline Matrix<R> restrict(const Matrix<R> & ) const;

  inline R trace() const;
  
  // UnivariatePoly<Z> char_poly() const;
  
  inline static Matrix<R> diagonalJoin(const std::vector< Matrix<R> > &);

  inline static Matrix<R> identity(size_t);
  
  // TODO - just change access resolution to the same vector instead
  inline Matrix<R> transpose() const;

  // arithmetic
  inline Matrix<R> operator*(const Matrix<R> &) const override;

  inline Matrix<R>& operator+=(const Matrix<R> &) override;
  inline Matrix<R>& operator-=(const Matrix<R> &) override;
  inline Matrix<R>& operator*=(const Matrix<R> &) override;

  inline Matrix<R> operator*(const R & a) const;

  // algorithms
  
  inline void swapRows(size_t, size_t);
  
  // in-place row-echelon form for the matrix echelon,
  // returns the rank and the transformation matrix trans
  inline static size_t rowEchelon(Matrix<R> & , Matrix<R>& );

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

  inline Matrix<R>* getPtr() override {return this;}

  inline const Matrix<R>* getPtr() const override {return this;}

  inline bool isZero() const override;
  inline bool isOne() const override;

  inline Matrix<R>& makeZero() override;
  inline Matrix<R>& makeOne() override;
  
protected:
  size_t _nrows;
  size_t _ncols;
  std::vector<R> _data;
};

template<typename R>
inline Matrix<R> operator*(const R & a, const Matrix<R> & mat)
{ return mat*a; }

#include "Matrix.inl"

#endif // __MATRIX_H_
