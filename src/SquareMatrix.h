#ifndef __SQUARE_MATRIX_H_
#define __SQUARE_MATRIX_H_

#include "birch.h"

// !! TODO - It might be useful to let SquareMatrix inherit
// from Ring (as it is a ring), and implement the overrides ?

template<class R, class Parent, size_t n>
class SquareMatrix
{
  
  static_assert(std::is_base_of<RingElement<R,Parent>, R>::value, "R template parameter must inherit from RingElement.");
  static_assert(std::is_base_of<Ring<Parent,R>, Parent>::value, "Parent template parameter must inherit from Ring.");
  
public:
  // c-tor
  SquareMatrix(std::shared_ptr<const Parent> base_ring) : _base(base_ring) {}

  SquareMatrix(const R mat[n][n]);
  SquareMatrix(const SquareMatrix<R,Parent,n> & mat);

  // assignment
  SquareMatrix<R,Parent,n> & operator=(const SquareMatrix<R,Parent,n> &);
  
  // access
  inline const R& operator()(size_t i, size_t j) const {return this->_mat[i][j]; }
  inline R& operator()(size_t i, size_t j) {return this->_mat[i][j];}

  inline std::shared_ptr<const Parent> baseRing(void) const { return this->_base; };

  // return the i-th row
  Vector<R,Parent,n> operator[](size_t i) const;
  
  // arithmetic
  SquareMatrix<R,Parent,n> operator*(const SquareMatrix<R,Parent,n>&) const;
  Vector<R,Parent,n> operator*(const Vector<R,Parent,n>& vec) const;
  SquareMatrix<R,Parent,n> operator*(const R &) const;
  SquareMatrix<R,Parent,n> operator/(const R &) const;

  // booleans
  bool operator==(const SquareMatrix<R,Parent,n>&) const;
  inline bool operator!=(const SquareMatrix<R,Parent,n>& other) const
  { return !((*this) == other);}
  
  // ordering of the matrices for Minkowski reduction
  bool operator<(const SquareMatrix<R,Parent,n>&) const;
  
  bool isUpperTriangular(void) const;
  bool isLowerTriangular(void) const;
  bool isSymmetric(void) const;
  bool isPositiveDefinite(void) const;
  
  // basic operations
  void setIdentity(void);
  
  SquareMatrix<R,Parent,n> transpose(void) const;
  // !! TODO - save the inverse and track it
  // to save computation
  SquareMatrix<R,Parent,n> inverse(void) const;
  SquareMatrix<R,Parent,n> adjugate(size_t dim) const;
  
  template<size_t m>
  SquareMatrix<R,Parent,m> submatrix(size_t idxs[m]) const;
  R determinant(void) const;

  Vector<R,Parent,n> solve(const Vector<R,Parent,n> & vec) const;

  // more complex operations that might be useful outside the class
  bool cholesky(SquareMatrix<R,Parent,n>& L, Vector<R,Parent,n> & D) const;
  bool ldl(SquareMatrix<R,Parent,n>& L, Vector<R,Parent,n> & D) const;
  
  // elementary operations
  void swapRows(size_t row1, size_t row2);
  void swapCols(size_t col1, size_t col2);
  void multiplyRow(size_t row, const R & val);
  void multiplyCol(size_t col, const R & val);
  void addRow(size_t row_to, size_t row_from, const R & val);
  void addCol(size_t col_to, size_t col_from, const R & val);

  // !! TODO - this might be only relevant for SquareMatrixInt
  SquareMatrix<R,Parent,n> hermiteForm(const R & d) const;
  // static functions
  // compute S[idx1]*F*S[idx2]^t
  // this is needed in this form for jordan decomposition
  // F and FParent will be the field of fractions of R, RParent
  template<class F, class FParent>
  static F innerProduct(const SquareMatrix<R,Parent,n> &,
			const SquareMatrix<F,FParent,n> &,
			size_t, size_t);
  
  // global constants
  static SquareMatrix<R,Parent,n> identity(std::shared_ptr<const Parent>);

  std::ostream & prettyPrint(std::ostream &, size_t upTo = n) const;
  
protected:
  std::shared_ptr<const Parent> _base;
  R _mat[n][n];
  
  // helper functions
  void _deepCopy(const R mat[n][n]);
  
  Vector<R,Parent,n> _forwardSubstitution(const Vector<R,Parent,n> & vec) const;
  Vector<R,Parent,n> _backwardSubstitution(const Vector<R,Parent,n> & vec) const;
  SquareMatrix<R,Parent,n> _inverseLowerTriangular(void) const;
  SquareMatrix<R,Parent,n> _inverseUpperTriangular(void) const;
};

// left multiplication
template<class R, class Parent, size_t n>
inline SquareMatrix<R,Parent,n> operator*(const R & a, const SquareMatrix<R,Parent,n> & M)
{ return M*a; }

// printing
template<class R, class Parent, size_t n>
std::ostream& operator<<(std::ostream&, const SquareMatrix<R,Parent,n>&);

// we put it outside the class to avoid partial specialization
// !! TODO - figure out a better way

template<typename R, typename S, typename T, size_t n>
std::shared_ptr< SquareMatrixFp<S,T,n> >
mod(const SquareMatrixInt<R,n> & a, std::shared_ptr< Fp<S,T> > GF);

#include "SquareMatrix.inl"

#endif // __SQUARE_MATRIX_H_
