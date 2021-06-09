#ifndef __SQUARE_MATRIX_H_
#define __SQUARE_MATRIX_H_

#include "birch.h"
#include "Rational.h"
#include "RingElement.h"
#include "Vector.h"

// !! TODO - It might be useful to let SquareMatrix inherit
// from Ring (as it is a ring), and implement the overrides ?

template<class R, class Parent, size_t n>
class SquareMatrix
{
  
  static_assert(std::is_base_of<RingElement<R,Parent>, R>::value);
  static_assert(std::is_base_of<Ring<Parent,R>, R>::value);
  
public:
  // c-tor
  SquareMatrix(std::shared_ptr<const Parent> base_ring) : _base(base_ring) {}

  SquareMatrix(const R mat[n][n]);
  SquareMatrix(const SquareMatrix<R,Parent,n> & mat);

  // assignment
  inline SquareMatrix<R,Parent,n> & operator=(const SquareMatrix<R,Parent,n> &);
  
  // access
  inline const R& operator()(size_t i, size_t j) const {return mat[i][j]; }
  inline R& operator()(size_t i, size_t j) {return mat[i][j];}

  // return the i-th row
  inline Vector<R,Parent,n> operator[](size_t i) const;
  
  // arithmetic
  inline SquareMatrix<R,Parent,n> operator*(const SquareMatrix<R,Parent,n>&) const;
  inline Vector<R,Parent,n> operator*(const Vector<R,Parent,n>& vec) const;
  inline SquareMatrix<R,Parent,n> operator*(const R &) const;
  inline SquareMatrix<R,Parent,n> operator/(const R &) const;

  // booleans
  inline bool operator==(const SquareMatrix<R,Parent,n>&) const;
  inline bool operator!=(const SquareMatrix<R,Parent,n>& other) const
  { return !((*this) == other);}
  
  // ordering of the matrices for Minkowski reduction
  inline bool operator<(const SquareMatrix<R,Parent,n>&) const;
  
  inline bool isUpperTriangular() const;
  inline bool isLowerTriangular() const;
  inline bool isSymmetric() const;
  inline bool isPositiveDefinite() const;
  
  // basic operations
  inline void setIdentity(void);
  
  inline SquareMatrix<R,Parent,n> transpose(void) const;
  // !! TODO - save the inverse and track it
  // to save computation
  inline SquareMatrix<R,Parent,n> inverse(void) const;
  inline SquareMatrix<R,Parent,n> adjugate(size_t dim) const;
  
  template<size_t m>
  inline SquareMatrix<R,Parent,m> submatrix(size_t idxs[m]) const;
  inline R determinant(void) const;

  inline Vector<R,Parent,n> solve(const Vector<R,Parent,n> & vec) const;

  // more complex operations that might be useful outside the class
  inline bool cholesky(SquareMatrix<R,Parent,n>& L, Vector<R,Parent,n> & D) const;
  inline bool ldl(SquareMatrix<R,Parent,n>& L, Vector<R,Parent,n> & D) const;
  
  // elementary operations
  inline void swapRows(size_t row1, size_t row2);
  inline void swapCols(size_t col1, size_t col2);
  inline void multiplyRow(size_t row, const R & val);
  inline void multiplyCol(size_t col, const R & val);
  inline void addRow(size_t row_to, size_t row_from, const R & val);
  inline void addCol(size_t col_to, size_t col_from, const R & val);

  inline SquareMatrix<R,Parent,n> hermiteForm(const R & d) const;
  // static functions
  // compute S[idx1]*F*S[idx2]^t
  // this is needed in this form for jordan decomposition
  // F and FParent will be the field of fractions of R, RParent
  template<class F, class FParent>
  inline static F innerProduct(const SquareMatrix<R,Parent,n> &,
			       const SquareMatrix<F,FParent,n> &,
			       size_t, size_t);
  
  // global constants
  inline static SquareMatrix<R,Parent,n> identity(std::shared_ptr<const Parent>);

  inline std::ostream & prettyPrint(std::ostream &, size_t upTo = n) const;
  
protected:
  std::shared_ptr<const Parent> _base;
  R mat[n][n];
  
  // helper functions
  inline void deepCopy(const R mat[n][n]);
  
  inline Vector<R,Parent,n> forwardSubstitution(const Vector<R,Parent,n> & vec) const;
  inline Vector<R,Parent,n> backwardSubstitution(const Vector<R,Parent,n> & vec) const;
  inline SquareMatrix<R,Parent,n> inverseLowerTriangular(void) const;
  inline SquareMatrix<R,Parent,n> inverseUpperTriangular(void) const;
};

// left multiplication
template<class R, class Parent, size_t n>
inline SquareMatrix<R,Parent,n> operator*(const R & a, const SquareMatrix<R,Parent,n> & M)
{ return M*a; }

// printing
template<class R, class Parent, size_t n>
inline std::ostream& operator<<(std::ostream&, const SquareMatrix<R,Parent,n>&);

namespace std
{
  template<class R, class Parent, size_t n>
  struct hash<Vector<R,Parent,n> >
  {
    Z64 operator()(const Vector<R,Parent,n>& vec) const
    {
      Z64 fnv = FNV_OFFSET;
      for (size_t i = 0; i < n; i++)
	fnv = (fnv ^ vec[i]) * FNV_PRIME;
            
      return fnv;
    }
  };
}

#include "SquareMatrix.inl"

#endif // __SQUARE_MATRIX_H_
