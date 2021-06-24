#ifndef __SQUARE_MATRIX_INT_H_
#define __SQUARE_MATRIX_INT_H_

#include <stdalign.h>

// We create this one to mkae matrix operations more efficient when the underlying class is integral
// At first tries to let it inherit from SquareMatrix<Integer<R>, IntegerRing<R>, n>
// but just access to elements was inefficient

template <typename R, size_t n>
class SquareMatrixInt
{
public:
  // c-tor
  SquareMatrixInt() = default;

  SquareMatrixInt(const R mat[n][n]);
  SquareMatrixInt(const SquareMatrixInt<R,n> & mat);
  SquareMatrixInt(const SquareMatrix<Integer<R>,IntegerRing<R>,n> & );

  // assignment
  SquareMatrixInt<R,n> & operator=(const SquareMatrixInt<R,n> &);
  
  // access
  inline const R& operator()(size_t i, size_t j) const {return this->_mat[i][j]; }
  inline R& operator()(size_t i, size_t j) {return this->_mat[i][j];}

  // return the i-th row
  VectorInt<R,n> operator[](size_t i) const;
  
  // arithmetic
  SquareMatrixInt<R,n> operator*(const SquareMatrixInt<R,n>&) const;
  VectorInt<R,n> operator*(const VectorInt<R,n>& vec) const;
  SquareMatrixInt<R,n> operator*(const R &) const;
  SquareMatrixInt<R,n> operator/(const R &) const;

  // booleans
  bool operator==(const SquareMatrixInt<R,n>&) const;
  inline bool operator!=(const SquareMatrixInt<R,n>& other) const
  { return !((*this) == other);}
  
  // ordering of the matrices for Minkowski reduction
  bool operator<(const SquareMatrixInt<R,n>&) const;
  
  bool isUpperTriangular(void) const;
  bool isLowerTriangular(void) const;
  bool isSymmetric(void) const;
  bool isPositiveDefinite(void) const;
  
  // basic operations
  void setIdentity(void);
  
  SquareMatrixInt<R,n> transpose(void) const;
  // !! TODO - save the inverse and track it
  // to save computation
  SquareMatrixInt<R,n> inverse(void) const;
  SquareMatrixInt<R,n> adjugate(size_t dim) const;
  
  template<size_t m>
  SquareMatrixInt<R,m> submatrix(size_t idxs[m]) const;
  R determinant(void) const;

  VectorInt<R,n> solve(const VectorInt<R,n> & vec) const;

  // more complex operations that might be useful outside the class
  bool cholesky(SquareMatrixInt<R,n>& L, VectorInt<R,n> & D) const;
  // should have sizeof(S) > 3*sizeof(R)+log(n)
  template <typename S>
  bool ldl(SquareMatrixInt<S,n>& L, VectorInt<S,n> & D) const;
  
  // elementary operations
  void swapRows(size_t row1, size_t row2);
  void swapCols(size_t col1, size_t col2);
  void multiplyRow(size_t row, const R & val);
  void multiplyCol(size_t col, const R & val);
  void addRow(size_t row_to, size_t row_from, const R & val);
  void addCol(size_t col_to, size_t col_from, const R & val);

  SquareMatrixInt<R,n> hermiteForm(const R & d) const;
  
  // global constants
  static SquareMatrixInt<R,n> identity();

  std::ostream & prettyPrint(std::ostream &, size_t upTo = n) const;
  
protected:
  alignas(32) R _mat[n][n];
  
  // helper functions
  void _deepCopy(const R mat[n][n]);
  
  VectorInt<R,n> _forwardSubstitution(const VectorInt<R,n> & vec) const;
  VectorInt<R,n> _backwardSubstitution(const VectorInt<R,n> & vec) const;
  SquareMatrixInt<R,n> _inverseLowerTriangular(void) const;
  SquareMatrixInt<R,n> _inverseUpperTriangular(void) const;
  
};

// left multiplication
template<typename R, size_t n>
inline SquareMatrixInt<R,n> operator*(const R & a, const SquareMatrixInt<R,n> & M)
{ return M*a; }

// printing
template<typename R, size_t n>
std::ostream& operator<<(std::ostream&, const SquareMatrixInt<R,n>&);

#include "SquareMatrixInt.inl"

#endif // __SQUARE_MATRIX_INT_H_
