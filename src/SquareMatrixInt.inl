// implementation file for SquareMatrixInl.h
#include <cassert>

// SquareMatrixInt

template<typename R, size_t n>
inline void SquareMatrixInt<R,n>::_deepCopy(const R mat[n][n])
{
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      this->_mat[i][j] = mat[i][j];
}

// c-tors
template<typename R, size_t n>
SquareMatrixInt<R,n>::SquareMatrixInt(const R mat[n][n])
{
  _deepCopy(mat);
}

template<typename R, size_t n>
SquareMatrixInt<R,n>::SquareMatrixInt(const SquareMatrixInt<R,n> & other)
{
  _deepCopy(other._mat);
}

template<typename R, size_t n>
SquareMatrixInt<R,n>::SquareMatrixInt(const SquareMatrix<Integer<R>,IntegerRing<R>,n> & other)
{
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      this->_mat[i][j] = other(i,j).num();
}

// access
template<typename R, size_t n>
inline VectorInt<R,n> SquareMatrixInt<R,n>::operator[](size_t i) const
{
  assert(i < n);
  VectorInt<R,n> v;
  for (size_t j = 0; j < n; j++)
    v[j] = (*this)(i,j);
  return v;
}

// assignment
template<typename R, size_t n>
inline SquareMatrixInt<R,n> &
SquareMatrixInt<R,n>::operator=(const SquareMatrixInt<R,n> & other)
{
  if (this !=  &other) {
    _deepCopy(other._mat);
  }
  return (*this);
}

// arithmetic

// matrix multiplication is a major bottleneck, hence we attempt to optimize it here

template<typename R, size_t n>
inline SquareMatrixInt<R,n>
SquareMatrixInt<R,n>::operator*(const SquareMatrixInt<R,n>& other) const
{
  SquareMatrixInt<R,n> prod;

  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++) {
      prod._mat[i][j] = 0;
      for (size_t k = 0; k < n; k++)
	prod._mat[i][j] += this->_mat[i][k]*other._mat[k][j];
    }
  
  return prod;
}

template<>
SquareMatrixInt<Z64, 4>
SquareMatrixInt<Z64,4>::operator*(const SquareMatrixInt<Z64,4>& B) const
{
  SquareMatrixInt<Z64,4> C;
  const SquareMatrixInt<Z64,4> & A = (*this);
  
  const __m256 BCx = _mm_load_ps((Z64*)&B._mat[0]);
  const __m256 BCy = _mm_load_ps((Z64*)&B._mat[1]);
  const __m256 BCz = _mm_load_ps((Z64*)&B._mat[2]);
  const __m256 BCw = _mm_load_ps((Z64*)&B._mat[3]);

  Z64* leftRowPointer = &A._mat[0];
  Z64* resultRowPointer = &C._mat[0];

  for (size_t i = 0; i < 4; ++i, leftRowPointer += 4, resultRowPointer += 4) {
    __m256 ARx = _mm_set1_ps(leftRowPointer[0]);
    __m256 ARy = _mm_set1_ps(leftRowPointer[1]);
    __m256 ARz = _mm_set1_ps(leftRowPointer[2]);
    __m256 ARw = _mm_set1_ps(leftRowPointer[3]);

    __mm256 X = ARx * BCx;
    __mm256 Y = ARy * BCy;
    __mm256 Z = ARz * BCz;
    __mm256 W = ARw * BCw;

    __mm256 R = X+Y+Z+W;
    _mm_store_ps(resultRowPointer, R);
  }
  
  return C;
}

template<typename R, size_t n>
inline VectorInt<R,n> SquareMatrixInt<R,n>::operator*(const VectorInt<R,n>& vec) const
{
  VectorInt<R,n> prod;
  for (size_t i = 0; i < n; i++) {
    prod[i] = 0;
    for (size_t j = 0; j < n; j++)
      prod[i] += this->_mat[i][j] * vec[j];
  }
  return prod;
}

template<typename R, size_t n>
inline SquareMatrixInt<R,n> SquareMatrixInt<R,n>::operator*(const R & scalar) const {
  SquareMatrixInt<R,n> prod;
  for (size_t row = 0; row < n; row++)
    for (size_t col = 0; col < n; col++)
      prod(row,col) = scalar*this->_mat[row][col];
  return prod;
}

template<typename R, size_t n>
inline SquareMatrixInt<R,n>  SquareMatrixInt<R,n>::operator/(const R & scalar) const
{
  assert (scalar != 0);

  SquareMatrixInt<R,n> quo;
  for (size_t row = 0; row < n; row++)
    for (size_t col = 0; col < n; col++)
      quo(row,col) = this->_mat[row][col] / scalar;
  return quo;
}

// booleans
template<typename R, size_t n>
bool SquareMatrixInt<R,n>::operator==(const SquareMatrixInt<R,n>& other) const
{
  for (size_t row = 0; row < n; row++)
    for (size_t col = 0; col < n; col++)
      if (this->_mat[row][col] != other(row, col)) return false;
  return true;
}

// !! - TODO - replace == and < by compare to save time
// Also, this is only relevant for integral matrices
template<typename R, size_t n>
inline bool SquareMatrixInt<R,n>::operator<(const SquareMatrixInt<R,n>& other) const
{
  for (size_t i = 0; i < n; i++) {
    if (this->_mat[i][i] < other(i,i)) return true;
    if (this->_mat[i][i] > other(i,i)) return false;
  }
  for (size_t col = 0; col < n-1; col++)
    for (size_t row = 0; row < n-1-col; row++) {
      if (abs(this->_mat[row][row+col+1]) > abs(other(row, row+col+1))) return true;
      if (abs(this->_mat[row][row+col+1]) < abs(other(row, row+col+1))) return false;
    }
  for (size_t col = 0; col < n-1; col++)
    for (size_t row = 0; row < n-1-col; row++) {
      if (this->_mat[row][row+col+1] > other(row, row+col+1)) return true;
      if (this->_mat[row][row+col+1] < other(row, row+col+1)) return false;
    }
  if (!isSymmetric()) {
    for (size_t col = 0; col < n-1; col++)
      for (size_t row = 0; row < n-1-col; row++) {
	if (this->_mat[row+col+1][col] > other(row+col+1, col)) return true;
	if (this->_mat[row+col+1][col] < other(row+col+1, col)) return false;
      }
  }
  return false;
}

template<typename R, size_t n>
inline bool SquareMatrixInt<R,n>::isUpperTriangular(void) const
{
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < i; j++)
      if (this->_mat[i][j] != 0) return false;
  return true;
}

template<typename R, size_t n>
inline bool SquareMatrixInt<R,n>::isLowerTriangular(void) const
{
  for (size_t i = 0; i < n; i++)
    for (size_t j = i+1; j < n; j++)
      if (this->_mat[i][j] != 0) return false;
  return true;
}

template<typename R, size_t n>
inline bool SquareMatrixInt<R,n>::isSymmetric(void) const
{
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < i; j++)
      if (this->_mat[i][j] != this->_mat[j][i]) return false;
  return true;
}

template<typename R, size_t n>
inline bool SquareMatrixInt<R,n>::isPositiveDefinite(void) const
{
  SquareMatrixInt<R,n> L;
  VectorInt<R,n> D;
  return ldl(L,D);
}
  
// basic operations
template<typename R, size_t n>
inline SquareMatrixInt<R,n> SquareMatrixInt<R,n>::transpose(void) const
{
  SquareMatrixInt<R,n> trans;
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      trans(i,j) = this->_mat[j][i];
  return trans;
}

template<typename R, size_t n>
template<size_t m>
inline SquareMatrixInt<R,m> SquareMatrixInt<R,n>::submatrix(size_t idxs[m]) const
{
  SquareMatrixInt<R,m> sub;
  for (size_t i = 0; i < m; i++)
    for (size_t j = 0; j < m; j++)
      sub(i,j) = (*this)(idxs[i], idxs[j]);
  
  return sub;
}

template<typename R, size_t n>
inline R SquareMatrixInt<R,n>::determinant(void) const
{
  // Instead of the previous ad-hoc method, we use Bareiss algorithm
  // to compute the determinant.
  // TODO - can do Cholesky, will be faster
  // !! TODO - Leibniz should be better when n <= 5 !?
  R sign = 1;
  SquareMatrixInt<R,n+1> M;
  M(0,0) = 1;
  // init
  for (size_t row = 0; row < n; row++)
    for (size_t col = 0; col < n; col++)
      M(row+1, col+1) = this->_mat[row][col];
  // we initialize the first row and column to 0, even though we don't use them
  // otherwise the compiler gets mad
  // !! TODO - use M00 as a single variable outside and have M as nXn matrix
  for (size_t row = 0; row < n; row++) {
    M(row+1,0) = 0;
    M(0,row+1) = 0;
  }
  
  for (size_t k = 1; k < n; k++) {
    if (M(k,k) == 0) {
      bool found = false;
      for (size_t i = k+1; i <= n; i++) {
	if (M(i,k) != 0) {
	  M.swapRows(k,i);
	  sign = -sign;
	  found = true;
	  break;
	}
      }
      if (!found)
	return 0;
    }
    for (size_t i = k+1; i <= n; i++)
      for (size_t j = k+1; j <= n; j++)
	M(i,j) = (M(i,j)*M(k,k) - M(i,k)*M(k,j))/M(k-1,k-1);
  }
  return sign*M(n,n);
}

// this solves using forward substitution in O(n^2)
// for lower triangular matrices
template<typename R, size_t n>
inline VectorInt<R,n>
SquareMatrixInt<R,n>::_forwardSubstitution(const VectorInt<R,n> & vec) const
{
  VectorInt<R,n> sol;
  R sum = 0;
  assert(isLowerTriangular());
 
  for (size_t i = 0; i < n; i++) {
    sum = 0;
    for (size_t j = 0; j < i; j++)
      sum += this->_mat[i][j] * sol[j];
    sol[i] = (vec[i] - sum) / this->_mat[i][i];
  } 
 
  return sol;
}

// this solves using forward substitution in O(n^2)
// for lower triangular matrices
template<typename R, size_t n>
inline SquareMatrixInt<R,n>
SquareMatrixInt<R,n>::_inverseLowerTriangular(void) const
{
  SquareMatrixInt<R,n> inv;
  inv.setIdentity();
  R sum = 0;
  assert(isLowerTriangular());

  for (size_t col = 0; col < n; col++) {
    for (size_t i = col; i < n; i++) {
      sum = 0;
      for (size_t j = 0; j < i; j++)
	sum += this->_mat[i][j] * inv(j, col);
      R delta = (i == col) ? 1 : 0;
      inv(i,col) = (delta - sum) / this->_mat[i][i];
    } 
  }
  return inv;
}

template<typename R, size_t n>
inline SquareMatrixInt<R,n>
SquareMatrixInt<R,n>::_inverseUpperTriangular(void) const
{
  SquareMatrixInt<R,n> inv;
  inv.setIdentity();
  R sum = 0;
  assert(isUpperTriangular());

  for (size_t col = 0; col < n; col++) {
    for (size_t i = col+1; i > 0; i--) {
      sum = 0;
      for (size_t j = i; j < n; j++)
	sum += this->_mat[i-1][j] * inv(j,col);
      R delta = (i-1 == col) ? 1 : 0;
      inv(i-1,col) = (delta - sum) / this->_mat[i-1][i-1];
    } 
  }
  return inv;
}

// this solves using backward substitution in O(n^2)
// for upper triangular matrices
template<typename R, size_t n>
inline VectorInt<R,n>
SquareMatrixInt<R,n>::_backwardSubstitution(const VectorInt<R,n> & vec) const
{
  VectorInt<R,n> sol;
  R sum = 0;
  assert(isUpperTriangular());
  
  for (size_t i = n; i > 0; i--) {
    sum = 0;
    for (size_t j = i; j < n; j++)
      sum += this->_mat[i-1][j] * sol[j];
    sol[i-1] = (vec[i-1] - sum) / this->_mat[i-1][i-1];
  } 
  
  return sol;
}

// returns false if the matrix is not positive definite
template<typename R, size_t n>
inline bool
SquareMatrixInt<R,n>::cholesky(SquareMatrixInt<R,n>& L,
			       VectorInt<R,n> & D) const
{
  assert(isSymmetric());
  L.setIdentity();
  R sum = 0;
  for (size_t j = 0; j < n; j++) {
    sum.makeZero();
    for (size_t k = 0; k < j; k++)
      sum += L(j,k)*L(j,k)*D[k];
    D[j] = this->_mat[j][j] - sum;
    if (D[j] == 0) return false;
    for (size_t i = j+1; i < n; i++) {
      sum = 0;
      for (size_t k = 0; k < j; k++)
	sum += L(i,k)*L(j,k)*D[k];
      L(i,j) = (this->_mat[i][j] - sum) / D[j];
    } 
  }
  /*
    #ifdef DEBUG
    // verify that L*Q*L^t = D
    SquareMatrixInt<R,n> diag(_base);
    for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
    diag(i,j) = (i == j) ? D[i] : Math<R>::zero();
    assert(L*(*this)*L.transpose() == diag);
    #endif
  */
  return true;
}

// !! TODO - This is integral Cholesky - should be able to
// determine which to call by the class R
template<typename R, size_t n>
inline bool
SquareMatrixInt<R,n>::ldl(SquareMatrixInt<R,n>& L,  VectorInt<R,n> & D) const
{
  assert(isSymmetric());
  R prod_diag = 1;
  R d = 0;
  R inner_sum = 0;
  // This works but inefficiently - for some reason we get O(n^4) operations.
  // !! TODO - check it out later
  // Oh I see - we should do the L update in two passes...
  for (size_t i = 0; i < n; i++)
    {
      L(i, i) = prod_diag;
      d = prod_diag;
      for (size_t j = 0; j < i; j++)
	{
	  L(i,j) = 0;
	  for (size_t k = j; k < i; k++)
	    {
	      inner_sum = 0;
	      for (size_t r = 0; r <= k; r++)
		inner_sum += L(k, r)*((*this)(i,r))*L(k,j);
	      inner_sum *= -L(i, i) / D[k];
	      L(i,j) += inner_sum;
	    }
	  //d = d.gcd(L(i, j));
	  d = birch_util::gcd(d, L(i,j));
	}
      for (size_t j = 0; j <= i; j++)
	L(i,j) /= d;
      D[i] = 0;
      for (size_t j = 0; j <= i; j++)
	for (size_t k = 0; k <= i; k++)
	  D[i] += L(i, j)*((*this)(j,k))*L(i, k);
      if (D[i] == 0) return false;
      // prod_diag = prod_diag.lcm(D[i]);
      prod_diag = birch_util::lcm(prod_diag, D[i]);
      for (size_t j = i+1; j < n; j++)
	L(i,j) = 0;
    }
  
#ifdef DEBUG
  // verify that L*Q*L^t = D
  SquareMatrixInt<R,n> diag;
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      diag(i,j) = (i == j) ? D[i] :0;
  assert(L*(*this)*L.transpose() == diag);
#endif
  return true;
}

// This solves only for symmetric positive definite matrices
// using LDL
template<typename R, size_t n>
inline VectorInt<R,n>
SquareMatrixInt<R,n>::solve(const VectorInt<R,n> & vec) const
{
  assert(isSymmetric());
  VectorInt<R,n> sol;
  SquareMatrixInt<R,n> L;
  VectorInt<R,n> D;
  bool is_positive_definite = cholesky(L, D);
  assert(is_positive_definite);
  sol = L._forwardSubstitution(vec);
  for (size_t i = 0; i < n; i++)
    sol[i] /= D[i];
  sol = L.transpose()._backwardSubstitution(sol);
  return sol;
}

template<typename R, size_t n>
inline void SquareMatrixInt<R,n>::swapRows(size_t row1, size_t row2)
{
  assert((row1 < n) && (row2 < n));
  
  VectorInt<R,n> temp_row;

  for (size_t col = 0; col < n; col++)
    temp_row[col] = this->_mat[row1][col];

  for (size_t col = 0; col < n; col++)
    this->_mat[row1][col] = this->_mat[row2][col];

  for (size_t col = 0; col < n; col++)
    this->_mat[row2][col] = temp_row[col];

  return;
}

template<typename R, size_t n>
inline void SquareMatrixInt<R,n>::swapCols(size_t col1, size_t col2)
{
  assert((col1 < n) && (col2 < n));

  VectorInt<R,n> temp_col;

  for (size_t row = 0; row < n; row++)
    temp_col[row] = this->_mat[row][col1];

  for (size_t row = 0; row < n; row++)
    this->_mat[row][col1] = this->_mat[row][col2];

  for (size_t row = 0; row < n; row++)
    this->_mat[row][col2] = temp_col[row];

  return;
}

template<typename R, size_t n>
inline void SquareMatrixInt<R,n>::multiplyRow(size_t row, const R & val)
{
  assert(row < n);

  for (size_t col = 0; col < n; col++)
    this->_mat[row][col] *= val;
  return;
}

template<typename R, size_t n>
inline void SquareMatrixInt<R,n>::multiplyCol(size_t col, const R & val)
{
  assert(col < n);

  for (size_t row = 0; row < n; row++)
    this->_mat[row][col] *= val;
  return;
}

template<typename R, size_t n>
inline void SquareMatrixInt<R,n>::addRow(size_t row_to, size_t row_from, const R & val)
{
  assert((row_to < n) && (row_from < n));

  for (size_t col = 0; col < n; col++) {
    this->_mat[row_to][col] += val * this->_mat[row_from][col];
  }
  return;
}

template<typename R, size_t n>
inline void SquareMatrixInt<R,n>::addCol(size_t col_to, size_t col_from, const R & val)
{
  assert((col_to < n) && (col_from < n));

  for (size_t row = 0; row < n; row++) {
    this->_mat[row][col_to] += val * this->_mat[row][col_from];
  }
  return;
}

template<typename R, size_t n>
inline SquareMatrixInt<R,n> SquareMatrixInt<R,n>::adjugate(size_t dim) const
{
  SquareMatrixInt<R,n> adj;
#ifdef DEBUG
  adj = SquareMatrixInt<R,n>::identity();
#endif
  // We will use this only for dim <= 4
  // and we write it down explicitly for each case
  // !! TODO - This is not the most effective way to do this
  const SquareMatrixInt<R,n> & a = (*this);
  if (dim <= n) {
    if (dim == 1) {
      adj(0,0) = 1;
    }
    if (dim == 2) {
      adj(0,0) = a(1,1);
      adj(1,1) = a(0,0);
      adj(0,1) = -a(0,1);
      adj(1,0) = -a(1,0);
    }
    if (dim == 3) {
      adj(0,0) = a(1,1)*a(2,2) - a(1,2)*a(2,1);
      adj(1,0) = - a(1,0)*a(2,2) + a(1,2)*a(2,0);
      adj(2,0) = a(1,0)*a(2,1) - a(1,1)*a(2,0);
      adj(0,1) = - a(0,1)*a(2,2) + a(2,1)*a(0,2);
      adj(1,1) = a(0,0)*a(2,2) - a(0,2)*a(2,0);
      adj(2,1) = - a(0,0)*a(2,1) + a(0,1)*a(2,0);
      adj(0,2) = a(0,1)*a(1,2) - a(1,1)*a(0,2);
      adj(1,2) = - a(0,0)*a(1,2) + a(1,0)*a(0,2);
      adj(2,2) = a(1,1)*a(0,0) - a(1,0)*a(0,1);
    }
    if (dim == 4) {
      adj(0,0) = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2));
      adj(0,0) -= a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1));
      adj(0,0) += a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1));
      adj(0,1) = -a(1,0)*(a(2,2)*a(3,3)-a(2,3)*a(3,2));
      adj(0,1) += a(1,2)*(a(2,0)*a(3,3)-a(2,3)*a(3,0));
      adj(0,1) -= a(1,3)*(a(2,0)*a(3,2)-a(2,2)*a(3,0));
      adj(0,2) = a(1,0)*(a(2,1)*a(3,3)-a(2,3)*a(3,1));
      adj(0,2) -= a(1,1)*(a(2,0)*a(3,3)-a(2,3)*a(3,0));
      adj(0,2) += a(1,3)*(a(2,0)*a(3,1)-a(2,1)*a(3,0));
      adj(0,3) = -a(1,0)*(a(2,1)*a(3,2)-a(2,2)*a(3,1));
      adj(0,3) += a(1,1)*(a(2,0)*a(3,2)-a(2,2)*a(3,0));
      adj(0,3) -= a(1,2)*(a(2,0)*a(3,1)-a(2,1)*a(3,0));
      adj(1,0) = -a(0,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2));
      adj(1,0) += a(2,1)*(a(0,2)*a(3,3)-a(3,2)*a(0,3));
      adj(1,0) -= a(3,1)*(a(0,2)*a(2,3)-a(2,2)*a(0,3));
      adj(1,1) = a(0,0)*(a(2,2)*a(3,3)-a(2,3)*a(3,2));
      adj(1,1) -= a(0,2)*(a(2,0)*a(3,3)-a(2,3)*a(3,0));
      adj(1,1) += a(0,3)*(a(2,0)*a(3,2)-a(2,2)*a(3,0));
      adj(1,2) = -a(0,0)*(a(2,1)*a(3,3)-a(2,3)*a(3,1));
      adj(1,2) += a(0,1)*(a(2,0)*a(3,3)-a(3,0)*a(2,3));
      adj(1,2) -= a(0,3)*(a(2,0)*a(3,1)-a(3,0)*a(2,1));
      adj(1,3) = a(0,0)*(a(2,1)*a(3,2)-a(3,1)*a(2,2));
      adj(1,3) -= a(0,1)*(a(2,0)*a(3,2)-a(3,0)*a(2,2));
      adj(1,3) += a(0,2)*(a(2,0)*a(3,1)-a(2,1)*a(3,0));
      adj(2,0) = a(0,1)*(a(1,2)*a(3,3)-a(3,2)*a(1,3));
      adj(2,0) -= a(1,1)*(a(0,2)*a(3,3)-a(3,2)*a(0,3));
      adj(2,0) += a(3,1)*(a(0,2)*a(1,3)-a(1,2)*a(0,3));
      adj(2,1) = -a(0,0)*(a(1,2)*a(3,3)-a(3,2)*a(1,3));
      adj(2,1) += a(1,0)*(a(0,2)*a(3,3)-a(0,3)*a(3,2));
      adj(2,1) -= a(3,0)*(a(0,2)*a(1,3)-a(0,3)*a(1,2));
      adj(2,2) = a(0,0)*(a(1,1)*a(3,3)-a(1,3)*a(3,1));
      adj(2,2) -= a(0,1)*(a(1,0)*a(3,3)-a(1,3)*a(3,0));
      adj(2,2) += a(0,3)*(a(1,0)*a(3,1)-a(1,1)*a(3,0));
      adj(2,3) = -a(0,0)*(a(1,1)*a(3,2)-a(1,2)*a(3,1));
      adj(2,3) += a(0,1)*(a(1,0)*a(3,2)-a(3,0)*a(1,2));
      adj(2,3) -= a(0,2)*(a(1,0)*a(3,1)-a(3,0)*a(1,1));
      adj(3,0) = -a(0,1)*(a(1,2)*a(2,3)-a(2,2)*a(1,3));
      adj(3,0) += a(1,1)*(a(0,2)*a(2,3)-a(2,2)*a(0,3));
      adj(3,0) -= a(2,1)*(a(0,2)*a(1,3)-a(1,2)*a(0,3));
      adj(3,1) = a(0,0)*(a(1,2)*a(2,3)-a(1,3)*a(2,2));
      adj(3,1) -= a(1,0)*(a(0,2)*a(2,3)-a(0,3)*a(2,2));
      adj(3,1) += a(2,0)*(a(0,2)*a(1,3)-a(1,2)*a(0,3));
      adj(3,2) = -a(0,0)*(a(1,1)*a(2,3)-a(2,1)*a(1,3));
      adj(3,2) += a(1,0)*(a(0,1)*a(2,3)-a(0,3)*a(2,1));
      adj(3,2) -= a(2,0)*(a(0,1)*a(1,3)-a(0,3)*a(1,1));
      adj(3,3) = a(0,0)*(a(1,1)*a(2,2)-a(1,2)*a(2,1));
      adj(3,3) -= a(0,1)*(a(1,0)*a(2,2)-a(1,2)*a(2,0));
      adj(3,3) += a(0,2)*(a(1,0)*a(2,1)-a(1,1)*a(2,0));
    }
  }
#ifdef DEBUG
  SquareMatrixInt<R,n> a_copy = SquareMatrixInt<R,n>::identity();
  for (size_t row = 0; row < dim; row++)
    for (size_t col = 0; col < dim; col++)
      a_copy(row,col) = a(row,col);
  R det = a_copy.determinant();
  SquareMatrixInt<R,n>  prod = adj*a;
  for (size_t row = 0; row < dim; row++)
    for (size_t col = 0; col < dim; col++)
      assert(prod(row,col) == ((row==col) ? det : 0));
#endif
  return adj;
}

// static functions

// Here we compute a highly specialized hermite form
// Assuming the matrix is non-singular
// and that we are looking for the hermite form of the matrix
// concatenated to the scalar matrix d*I
template<typename R, size_t n>
inline SquareMatrixInt<R,n> SquareMatrixInt<R,n>::hermiteForm(const R & d) const
{
  R a = 0;
  R b = 0;
  R x = 0;
  R y = 0;
  R g = 0;
  R q_a = 0;
  R q_b = 0;

  SquareMatrixInt<R,n> H = d*SquareMatrixInt<R,n>::identity();
  for (size_t row = 0; row < n; row++) {
    VectorInt<R,n> b_primes = (*this)[row];
    // Here we compute the HNF of H and this row and store it in H
    for (size_t pivot = 0; pivot < n; pivot++) {
      a = H(pivot,pivot);
      b = b_primes[pivot];
      // std::tuple<R,R,R> g_x_y = a.xgcd(b);
      //      g = std::get<0>(g_x_y);
      // x = std::get<1>(g_x_y);
      // y = std::get<2>(g_x_y);
      // for now,  we use gmp gcdext. can replace with a more efficient one for fixed precision.
      Z g_Z,x_Z,y_Z,a_Z,b_Z;
      a_Z = birch_util::convertInteger<R,Z>(a);
      b_Z = birch_util::convertInteger<R,Z>(b);
      mpz_gcdext(g_Z.get_mpz_t(),x_Z.get_mpz_t(),y_Z.get_mpz_t(),a_Z.get_mpz_t(), b_Z.get_mpz_t());
      a = birch_util::convertInteger<Z,R>(a_Z);
      b = birch_util::convertInteger<Z,R>(b_Z);
      g = birch_util::convertInteger<Z,R>(g_Z);
      x = birch_util::convertInteger<Z,R>(x_Z);
      y = birch_util::convertInteger<Z,R>(y_Z);
      q_a = a / g;
      q_b = b / g;
      VectorInt<R,n> g_h_prime = x*H[pivot] + y*b_primes;
      b_primes = q_a*b_primes-q_b*H[pivot];
      for (size_t j = pivot; j < n; j++) {
	R scalar = b_primes[j] / H(j,j);
	b_primes -= scalar*H[j];
      }
      for (size_t col = 0; col < n; col++)
	H(pivot, col) = g_h_prime[col];
    }
    for (size_t pivot = n-1; pivot > 0; pivot--) {
      for (size_t col = pivot; col < n; col++) { 
	R q = H(pivot-1,col) / H(col, col);
	for (size_t j = col; j < n; j++)
	  H(pivot-1, j) -= q*H(col, j);
      }
    }
  }
  return H;
}

// global constants
template<typename R, size_t n>
inline SquareMatrixInt<R,n>
SquareMatrixInt<R,n>::identity()
{
  SquareMatrixInt<R,n> id;
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      id(i,j) = (i == j) ? 1 : 0;
  return id; 
}

// printing
template<typename R, size_t n>
inline std::ostream& operator<<(std::ostream& os, const SquareMatrixInt<R,n>& a)
{
  os << "Matrix(Integers(), " << n << " , (";
  for (size_t row = 0; row < n-1; row++) {
    for (size_t col = 0; col < n; col++)
      os << a(row, col) << ",";
  }
  for (size_t col = 0; col < n-1; col++)
    os << a(n-1,col) << ",";
  os << a(n-1,n-1) <<  "))";
    
  return os;
}

template<typename R, size_t n>
inline std::ostream & SquareMatrixInt<R,n>::prettyPrint(std::ostream & os,
							size_t upTo) const
{
  for (size_t row = 0; row < upTo-1; row++) {
    for (size_t col = 0; col < upTo; col++)
      os << (*this)(row, col) << " ";
    os << std::endl;
  }
  for (size_t col = 0; col < upTo-1; col++)
    os << (*this)(upTo-1,col) << " ";
  os << (*this)(upTo-1,upTo-1) <<  std::endl;
    
  return os;
}

template<typename R, size_t n>
inline void SquareMatrixInt<R,n>::setIdentity(void)
{
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      (*this)(i,j) = (i == j) ? 1 : 0;
  
  return;
}

