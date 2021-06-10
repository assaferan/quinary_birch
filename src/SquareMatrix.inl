// implementation file for SquareMatrix.h

#include <cassert>

// SquareMatrix

template<class R, class Parent, size_t n>
void SquareMatrix<R,Parent,n>::deepCopy(const R mat[n][n])
{
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      this->mat[i][j] = mat[i][j];
}

// c-tors
template<class R, class Parent, size_t n>
SquareMatrix<R,Parent,n>::SquareMatrix(const R mat[n][n])
{
  deep_copy(mat);
  if (n > 0)
    _base = mat[0][0].parent();
}

template<class R, class Parent, size_t n>
SquareMatrix<R,Parent,n>::SquareMatrix(const SquareMatrix<R,Parent,n> & other)
{
  deepCopy(other.mat);
  _base = other._base;
}

// access
template<class R, class Parent, size_t n>
Vector<R,Parent,n> SquareMatrix<R,Parent,n>::operator[](size_t i) const
{
  assert(i < n);
  Vector<R,Parent,n> v(_base);
  for (size_t j = 0; j < n; j++)
    v[j] = (*this)(i,j);
  return v;
}

// assignment
template<class R, class Parent, size_t n>
SquareMatrix<R,Parent,n> &
SquareMatrix<R,Parent,n>::operator=(const SquareMatrix<R,Parent,n> & other)
{
  if (this !=  &other) {
    deepCopy(other.mat);
    _base = other._base;
  }
  return (*this);
}

// arithmetic

template<class R, class Parent, size_t n>
SquareMatrix<R,Parent,n>
SquareMatrix<R,Parent,n>::operator*(const SquareMatrix<R,Parent,n>& other) const
{
  SquareMatrix<R,Parent,n> prod(_base);
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++) {
      prod(i,j) = _base->zero();
      for (size_t k = 0; k < n; k++)
	prod(i,j) += this->mat[i][k]*other(k,j);
    }
  return prod;
}

template<class R, class Parent, size_t n>
Vector<R,Parent,n> SquareMatrix<R,Parent,n>::operator*(const Vector<R,Parent,n>& vec) const
{
  Vector<R,Parent,n> prod(_base);
  for (size_t i = 0; i < n; i++) {
    prod[i] = _base->zero();
    for (size_t j = 0; j < n; j++)
      prod[i] += this->mat[i][j] * vec[j];
  }
  return prod;
}

template<class R, class Parent, size_t n>
SquareMatrix<R,Parent,n> SquareMatrix<R,Parent,n>::operator*(const R & scalar) const {
  SquareMatrix<R,Parent,n> prod(_base);
  for (size_t row = 0; row < n; row++)
    for (size_t col = 0; col < n; col++)
      prod(row,col) = scalar*this->mat[row][col];
  return prod;
}

template<class R, class Parent, size_t n>
SquareMatrix<R,Parent,n>  SquareMatrix<R,Parent,n>::operator/(const R & scalar) const
{
  assert(scalar != 0);

  SquareMatrix<R,Parent,n> quo;
  for (size_t row = 0; row < n; row++)
    for (size_t col = 0; col < n; col++)
      quo(row,col) = this->mat[row][col] / scalar;
  return quo;
}

// booleans
template<class R, class Parent, size_t n>
bool SquareMatrix<R,Parent,n>::operator==(const SquareMatrix<R,Parent,n>& other) const
{
  for (size_t row = 0; row < n; row++)
    for (size_t col = 0; col < n; col++)
      if (mat[row][col] != other(row, col)) return false;
  return true;
}

// !! - TODO - replace == and < by compare to save time
// Also, this is only relevant for integral matrices
template<class R, class Parent, size_t n>
bool SquareMatrix<R,Parent,n>::operator<(const SquareMatrix<R,Parent,n>& other) const
{
  for (size_t i = 0; i < n; i++) {
    if (mat[i][i] < other(i,i)) return true;
    if (mat[i][i] > other(i,i)) return false;
  }
  for (size_t col = 0; col < n-1; col++)
    for (size_t row = 0; row < n-1-col; row++) {
      if (mat[row][row+col+1].abs() > other(row, row+col+1).abs()) return true;
      if (mat[row][row+col+1].abs() < other(row, row+col+1).abs()) return false;
    }
  for (size_t col = 0; col < n-1; col++)
    for (size_t row = 0; row < n-1-col; row++) {
      if (mat[row][row+col+1] > other(row, row+col+1)) return true;
      if (mat[row][row+col+1] < other(row, row+col+1)) return false;
    }
  if (!isSymmetric()) {
    for (size_t col = 0; col < n-1; col++)
      for (size_t row = 0; row < n-1-col; row++) {
	if (mat[row+col+1][col] > other(row+col+1, col)) return true;
	if (mat[row+col+1][col] < other(row+col+1, col)) return false;
      }
  }
  return false;
}

template<class R, class Parent, size_t n>
bool SquareMatrix<R,Parent,n>::isUpperTriangular() const
{
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < i; j++)
      if (!mat[i][j].isZero()) return false;
  return true;
}

template<class R, class Parent, size_t n>
bool SquareMatrix<R,Parent,n>::isLowerTriangular() const
{
  for (size_t i = 0; i < n; i++)
    for (size_t j = i+1; j < n; j++)
      if (!mat[i][j].isZero()) return false;
  return true;
}

template<class R, class Parent, size_t n>
bool SquareMatrix<R,Parent,n>::isSymmetric() const
{
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < i; j++)
      if (mat[i][j] != mat[j][i]) return false;
  return true;
}

template<class R, class Parent, size_t n>
bool SquareMatrix<R,Parent,n>::isPositiveDefinite() const
{
  SquareMatrix<R,Parent,n> L;
  Vector<R,Parent,n> D;
  return ldl(L,D);
}
  
// basic operations
template<class R, class Parent, size_t n>
SquareMatrix<R,Parent,n> SquareMatrix<R,Parent,n>::transpose(void) const
{
  SquareMatrix<R,Parent,n> trans(_base);
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      trans(i,j) = this->mat[j][i];
  return trans;
}

template<class R, class Parent, size_t n>
template<size_t m>
SquareMatrix<R,Parent,m> SquareMatrix<R,Parent,n>::submatrix(size_t idxs[m]) const
{
  SquareMatrix<R,Parent,m> sub(_base);
  for (size_t i = 0; i < m; i++)
    for (size_t j = 0; j < m; j++)
      sub(i,j) = (*this)(idxs[i], idxs[j]);
  
  return sub;
}

template<class R, class Parent, size_t n>
R SquareMatrix<R,Parent,n>::determinant(void) const
{
  // Instead of the previous ad-hoc method, we use Bareiss algorithm
  // to compute the determinant.
  // TODO - can do Cholesky, will be faster
  // !! TODO - Leibniz should be better when n <= 5 !?
  R sign = _base->one();
  SquareMatrix<R,Parent,n+1> M(_base);
  M(0,0) = _base->one();
  // init
  for (size_t row = 0; row < n; row++)
    for (size_t col = 0; col < n; col++)
      M(row+1, col+1) = this->mat[row][col];
  for (size_t k = 1; k < n; k++) {
    if (M(k,k).isZero()) {
      bool found = false;
      for (size_t i = k+1; i <= n; i++) {
	if (!M(i,k).isZero()) {
	  M.swapRows(k,i);
	  sign = -sign;
	  found = true;
	  break;
	}
      }
      if (!found)
	return _base->zero();
    }
    for (size_t i = k+1; i <= n; i++)
      for (size_t j = k+1; j <= n; j++)
	M(i,j) = (M(i,j)*M(k,k) - M(i,k)*M(k,j))/M(k-1,k-1);
  }
  return sign*M(n,n);
}

// this solves using forward substitution in O(n^2)
// for lower triangular matrices
template<class R, class Parent, size_t n>
Vector<R,Parent,n>
SquareMatrix<R,Parent,n>::forwardSubstitution(const Vector<R,Parent,n> & vec) const
{
  Vector <R,Parent,n> sol;
  R sum = _base->zero();
  assert(isLowerTriangular());
 
  for (size_t i = 0; i < n; i++) {
    sum = _base->zero();
    for (size_t j = 0; j < i; j++)
      sum += mat[i][j] * sol[j];
    sol[i] = (vec[i] - sum) / mat[i][i];
  } 
 
  return sol;
}

// this solves using forward substitution in O(n^2)
// for lower triangular matrices
template<class R, class Parent, size_t n>
SquareMatrix<R,Parent,n>
SquareMatrix<R,Parent,n>::inverseLowerTriangular(void) const
{
  SquareMatrix<R,Parent,n> inv(_base);
  inv.setIdentity();
  R sum = _base->zero();
  assert(isLowerTriangular());

  for (size_t col = 0; col < n; col++) {
    for (size_t i = col; i < n; i++) {
      sum.makeZero();
      for (size_t j = 0; j < i; j++)
	sum += mat[i][j] * inv(j, col);
      R delta = (i == col) ? _base->one() : _base->zero();
      inv(i,col) = (delta - sum) / mat[i][i];
    } 
  }
  return inv;
}

template<class R, class Parent, size_t n>
SquareMatrix<R,Parent,n>
SquareMatrix<R,Parent,n>::inverseUpperTriangular(void) const
{
  SquareMatrix<R,Parent,n> inv(_base);
  inv.setIdentity();
  R sum = _base->zero();
  assert(isUpperTriangular());

  for (size_t col = 0; col < n; col++) {
    for (size_t i = col+1; i > 0; i--) {
      sum = _base->zero();
      for (size_t j = i; j < n; j++)
	sum += mat[i-1][j] * inv(j,col);
      R delta = (i-1 == col) ? _base->one() : _base->zero();
      inv(i-1,col) = (delta - sum) / mat[i-1][i-1];
    } 
  }
  return inv;
}

// this solves using backward substitution in O(n^2)
// for upper triangular matrices
template<class R, class Parent, size_t n>
Vector<R,Parent,n>
SquareMatrix<R,Parent,n>::backwardSubstitution(const Vector<R,Parent,n> & vec) const
{
  Vector <R,Parent,n> sol(_base);
  R sum = _base->zero();
  assert(isUpperTriangular());
  
  for (size_t i = n; i > 0; i--) {
    sum = _base->zero();
    for (size_t j = i; j < n; j++)
      sum += mat[i-1][j] * sol[j];
    sol[i-1] = (vec[i-1] - sum) / mat[i-1][i-1];
  } 
  
  return sol;
}

// returns false if the matrix is not positive definite
template<class R, class Parent, size_t n>
bool
SquareMatrix<R,Parent,n>::cholesky(SquareMatrix<R,Parent,n>& L,
				   Vector<R,Parent,n> & D) const
{
  assert(isSymmetric());
  L.setIdentity();
  R sum = _base->zero();
  for (size_t j = 0; j < n; j++) {
    sum.makeZero();
    for (size_t k = 0; k < j; k++)
      sum += L(j,k)*L(j,k)*D[k];
    D[j] = mat[j][j] - sum;
    if (D[j].isZero()) return false;
    for (size_t i = j+1; i < n; i++) {
      sum.makeZero();
      for (size_t k = 0; k < j; k++)
	sum += L(i,k)*L(j,k)*D[k];
      L(i,j) = (mat[i][j] - sum) / D[j];
    } 
  }
  /*
    #ifdef DEBUG
    // verify that L*Q*L^t = D
    SquareMatrix<R,n> diag;
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
template<class R, class Parent, size_t n>
bool
SquareMatrix<R,Parent,n>::ldl(SquareMatrix<R,Parent,n>& L,  Vector<R,Parent,n> & D) const
{
  assert(isSymmetric());
  R prod_diag = _base->one();
  R d = _base->zero();
  R inner_sum = _base->zero();
  // This works but inefficiently - for some reason we get O(n^4) operations.
  // !! TODO - check it out later
  // Oh I see - we should do the L update in two passes...
  for (size_t i = 0; i < n; i++)
    {
      L(i, i) = prod_diag;
      d = prod_diag;
      for (size_t j = 0; j < i; j++)
	{
	  L(i, j).makeZero();
	  for (size_t k = j; k < i; k++)
	    {
	      inner_sum.makeZero();
	      for (size_t r = 0; r <= k; r++)
		inner_sum += L(k, r)*((*this)(i,r))*L(k,j);
	      inner_sum *= -L(i, i) / D[k];
	      L(i,j) += inner_sum;
	    }
	  d = d.gcd(L(i, j));
	}
      for (size_t j = 0; j <= i; j++)
	L(i,j) /= d;
      D[i].makeZero();
      for (size_t j = 0; j <= i; j++)
	for (size_t k = 0; k <= i; k++)
	  D[i] += L(i, j)*((*this)(j,k))*L(i, k);
      if (D[i].isZero()) return false;
      prod_diag = prod_diag.lcm(D[i]);
      for (size_t j = i+1; j < n; j++)
	L(i, j).makeZero();
    }
  
#ifdef DEBUG
  // verify that L*Q*L^t = D
  SquareMatrix<R,Parent,n> diag(_base);
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      diag(i,j) = (i == j) ? D[i] : _base->zero();
  assert(L*(*this)*L.transpose() == diag);
#endif
  return true;
}

// This solves only for symmetric positive definite matrices
// using LDL
template<class R, class Parent, size_t n>
Vector<R,Parent,n>
SquareMatrix<R,Parent,n>::solve(const Vector<R,Parent,n> & vec) const
{
  assert(isSymmetric());
  Vector<R,Parent,n> sol;
  SquareMatrix<R,Parent,n> L;
  Vector<R,Parent,n> D;
  bool is_positive_definite = cholesky(L, D);
  assert(is_positive_definite);
  sol = L.forwardSubstitution(vec);
  for (size_t i = 0; i < n; i++)
    sol[i] /= D[i];
  sol = L.transpose().backwardSubstitution(sol);
  return sol;
}

template<class R, class Parent, size_t n>
void SquareMatrix<R,Parent,n>::swapRows(size_t row1, size_t row2)
{
  assert((row1 < n) && (row2 < n));
  
  Vector<R,Parent,n> temp_row(_base);

  for (size_t col = 0; col < n; col++)
    temp_row[col] = mat[row1][col];

  for (size_t col = 0; col < n; col++)
    mat[row1][col] = mat[row2][col];

  for (size_t col = 0; col < n; col++)
    mat[row2][col] = temp_row[col];

  return;
}

template<class R, class Parent, size_t n>
void SquareMatrix<R,Parent,n>::swapCols(size_t col1, size_t col2)
{
  assert((col1 < n) && (col2 < n));

  Vector<R,Parent,n> temp_col(_base);

  for (size_t row = 0; row < n; row++)
    temp_col[row] = mat[row][col1];

  for (size_t row = 0; row < n; row++)
    mat[row][col1] = mat[row][col2];

  for (size_t row = 0; row < n; row++)
    mat[row][col2] = temp_col[row];

  return;
}

template<class R, class Parent, size_t n>
void SquareMatrix<R,Parent,n>::multiplyRow(size_t row, const R & val)
{
  assert(row < n);

  for (size_t col = 0; col < n; col++)
    mat[row][col] *= val;
  return;
}

template<class R, class Parent, size_t n>
void SquareMatrix<R,Parent,n>::multiplyCol(size_t col, const R & val)
{
  assert(col < n);

  for (size_t row = 0; row < n; row++)
    mat[row][col] *= val;
  return;
}

template<class R, class Parent, size_t n>
void SquareMatrix<R,Parent,n>::addRow(size_t row_to, size_t row_from, const R & val)
{
  assert((row_to < n) && (row_from < n));

  for (size_t col = 0; col < n; col++) {
    mat[row_to][col] += val * mat[row_from][col];
  }
  return;
}

template<class R, class Parent, size_t n>
void SquareMatrix<R,Parent,n>::addCol(size_t col_to, size_t col_from, const R & val)
{
  assert((col_to < n) && (col_from < n));

  for (size_t row = 0; row < n; row++) {
    mat[row][col_to] += val * mat[row][col_from];
  }
  return;
}
  
// a general one, just in case
template<class R, class Parent, size_t n>
SquareMatrix<R,Parent,n> SquareMatrix<R,Parent,n>::inverse(void) const
{  
  assert(this->determinant() != _base->zero());

  if (isLowerTriangular()) return inverseLowerTriangular();
  if (isUpperTriangular()) return inverseUpperTriangular();
  if (isSymmetric()) {
    SquareMatrix<R,Parent,n> L(_base);
    Vector<R,Parent,n> D(_base);
    bool is_positive_definite = cholesky(L,D);
    if (is_positive_definite) {
      SquareMatrix<R,Parent,n> L_inv = L.inverse();
      SquareMatrix<R,Parent,n> L_inv_t = L_inv.transpose();
      for (size_t i = 0; i < n; i++)
	for (size_t j = 0; j < n; j++)
	  L_inv(i,j) /= D[i];
      
      return L_inv_t * L_inv;
    }
  }
  SquareMatrix<R,Parent,n> inv(_base);
  inv.setIdentity();
  SquareMatrix<R,Parent,n> echelon(mat);
  size_t pivot_row = 0;
  size_t pivot_col = 0;
  size_t row_max;
  R max_val = _base->zero();
  R factor = _base->one();
  
  while ((pivot_row < n) && (pivot_col < n)) {
    row_max = pivot_row;
    max_val = echelon(row_max, pivot_col);
    if (max_val.isZero()) {
      pivot_col++;
    }
    else {
      echelon.swapRows(pivot_row, row_max);
      inv.swapRows(pivot_row, row_max);

      assert(inv*(*this) == echelon);

      R scalar = _base->one() / echelon(pivot_row, pivot_col);
      echelon.multiplyRow(pivot_row, scalar);
      inv.multiplyRow(pivot_row, scalar);

      assert(inv*(*this) == echelon);

      // for reduced row echelon form we need also the rows before
      for (size_t row = 0; row < pivot_row; row++) {
	factor = echelon(row, pivot_col);
	echelon(row, pivot_col).makeZero();
	for (size_t col = pivot_col + 1; col < n; col++) {
	  echelon(row,col) -= factor * echelon(pivot_row, col);
	}
	for (size_t col = 0; col < n; col++) {
	  inv(row, col) -= factor * inv(pivot_row, col);
	}

	assert(inv*(*this) == echelon);

      }
      for (size_t row = pivot_row+1; row < n; row++) {
	// factor = echelon(row,pivot_col) / echelon(pivot_row, pivot_col);
	factor = echelon(row, pivot_col);
	echelon(row, pivot_col).makeZero();
	for (size_t col = pivot_col + 1; col < n; col++) {
	  echelon(row,col) -= factor * echelon(pivot_row, col);
	}
	for (size_t col = 0; col < n; col++) {
	  inv(row, col) -= factor * inv(pivot_row, col);
	}
	assert(inv*(*this) == echelon);
      }
      
      pivot_row++;
      pivot_col++;
    }
  }
  assert(inv*(*this) == echelon);
  assert((*this)*inv == identity(_base));

  return inv;
}

template<class R, class Parent, size_t n>
SquareMatrix<R,Parent,n> SquareMatrix<R,Parent,n>::adjugate(size_t dim) const
{
  SquareMatrix<R,Parent,n> adj(_base);
#ifdef DEBUG
  adj = SquareMatrix<R,Parent,n>::identity(_base);
#endif
  // We will use this only for dim <= 4
  // and we write it down explicitly for each case
  // !! TODO - This is not the most effective way to do this
  const SquareMatrix<R,Parent,n> & a = (*this);
  if (dim <= n) {
    if (dim == 1) {
      adj(0,0).makeOne();
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
  SquareMatrix<R,Parent,n> a_copy = SquareMatrix<R,Parent,n>::identity(_base);
  for (size_t row = 0; row < dim; row++)
    for (size_t col = 0; col < dim; col++)
      a_copy(row,col) = a(row,col);
  R det = a_copy.determinant();
  SquareMatrix<R,Parent,n>  prod = adj*a;
  for (size_t row = 0; row < dim; row++)
    for (size_t col = 0; col < dim; col++)
      assert(prod(row,col) == ((row==col) ? det : _base->zero()));
#endif
  return adj;
}

// static functions
template<class R, class Parent, size_t n>
template<class F, class FParent>
F SquareMatrix<R,Parent,n>::innerProduct(const SquareMatrix<R,Parent,n> & G,
					 const SquareMatrix<F,FParent,n> & S,
					 size_t idx1, size_t idx2)
{
  assert((idx1 < n) && (idx2 < n));
  
  F ans = S._base->zero();
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      ans += S(idx1, i) * G(i,j) * S(idx2, j);
  return ans;
}

// Here we compute a highly specialized hermite form
// Assuming the matrix is non-singular
// and that we are looking for the hermite form of the matrix
// concatenated to the scalar matrix d*I
template<class R, class Parent, size_t n>
SquareMatrix<R,Parent,n> SquareMatrix<R,Parent,n>::hermiteForm(const R & d) const
{
  R a = _base->zero();
  R b = _base->zero();
  R x = _base->zero();
  R y = _base->zero();
  R g = _base->zero();
  R q_a = _base->zero();
  R q_b = _base->zero();

  SquareMatrix<R,Parent,n> H = d*SquareMatrix<R,Parent,n>::identity(_base);
  for (size_t row = 0; row < n; row++) {
    Vector<R,Parent,n> b_primes = (*this)[row];
    // Here we compute the HNF of H and this row and store it in H
    for (size_t pivot = 0; pivot < n; pivot++) {
      a = H(pivot,pivot);
      b = b_primes[pivot];
      std::tuple<R,R,R> g_x_y = a.xgcd(b);
      g = std::get<0>(g_x_y);
      x = std::get<1>(g_x_y);
      y = std::get<2>(g_x_y);
      q_a = a / g;
      q_b = b / g;
      Vector<R,Parent,n> g_h_prime = x*H[pivot] + y*b_primes;
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
template<class R, class Parent, size_t n>
SquareMatrix<R,Parent,n>
SquareMatrix<R,Parent,n>::identity(std::shared_ptr<const Parent> ring)
{
  SquareMatrix<R,Parent,n> id(ring);
  R one = ring->one();
  R zero = ring->zero();
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      id(i,j) = (i == j) ? one : zero;
  return id; 
}

// printing
template<class R, class Parent, size_t n>
std::ostream& operator<<(std::ostream& os, const SquareMatrix<R,Parent,n>& a)
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

template<class R, class Parent, size_t n>
std::ostream & SquareMatrix<R,Parent,n>::prettyPrint(std::ostream & os,
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

template<class R, class Parent, size_t n>
void SquareMatrix<R,Parent,n>::setIdentity(void)
{
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      (*this)(i,j) = (i == j) ? _base->one() : _base->zero();
  
  return;
}
