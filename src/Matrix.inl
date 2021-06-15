// For implementations

template<class R, class Parent>
template <size_t n>
Matrix<R,Parent>::Matrix(const R data[n][n])
  : _nrows(n), _ncols(n), _data(n*n)
{
  size_t idx = 0;
  for (size_t row = 0; row < _nrows; row++)
    for (size_t col = 0; col < _ncols; col++)
      _data[idx++] = data[row][col];

  if (n > 0)
    _base = data[0][0].parent();
}

template<class R, class Parent>
template <size_t n>
Matrix<R,Parent>::Matrix(const SquareMatrix<R,Parent,n> & mat)
  : _nrows(n), _ncols(n), _data(n*n), _base(mat.baseRing())
{
  size_t idx = 0;
  for (size_t row = 0; row < _nrows; row++)
    for (size_t col = 0; col < _ncols; col++)
      _data[idx++] = mat(row, col);
}

// return the i-th row
template<class R, class Parent>
inline std::vector<R> Matrix<R,Parent>::operator[](size_t i) const
{
  assert(i < _nrows);

  std::vector<R> v(_ncols, _base->zero());
  for (size_t j = 0; j < _ncols; j++)
    v[j] = (*this)(i,j);
  return v;
}

template<class R, class Parent>
inline R Matrix<R,Parent>::determinant(void) const
{		
  assert(_nrows == _ncols);
  size_t n = _nrows;
  Matrix<R,Parent> M(_base, n+1, n+1);
  M(0,0).makeOne();
  // init
  for (size_t row = 0; row < n; row++)
    for (size_t col = 0; col < n; col++)
      M(row+1, col+1) = (*this)(row, col);
  for (size_t k = 1; k < n; k++)
    for (size_t i = k+1; i <= n; i++)
      for (size_t j = k+1; j <= n; j++) {
	if (M(k-1,k-1).isZero()) return M(k-1,k-1);
	M(i,j) = (M(i,j)*M(k,k) - M(i,k)*M(k,j))/M(k-1,k-1);
      }
  return M(n,n);
}

template<class R, class Parent>
inline void Matrix<R,Parent>::swapRows(size_t row1, size_t row2)
{
  R tmp = _base->zero();
  for (size_t col = 0; col < this->_ncols; col++) {
    tmp = (*this)(row1, col);
    (*this)(row1, col) = (*this)(row2, col);
    (*this)(row2, col) = tmp;
  }
  return;
}

// in place echelon form, returns the rank and trans is the transformation
template<class R, class Parent>
inline size_t Matrix<R,Parent>::rowEchelon(Matrix<R,Parent> & echelon, Matrix<R,Parent>& trans)
{
  // This one assumes R is a field
  // !! TODO - think how to make this method appear only for fields
  static_assert(std::is_base_of<FieldElement<R,Parent>, R>::value);
  
  // trans = identity(echelon.nrows());
  for (size_t row = 0; row < trans.nrows(); row++)
    for (size_t col = 0; col < trans.ncols(); col++)
      trans(row, col) = (row == col) ? 1 : 0;
  
  size_t pivot_row = 0;
  size_t pivot_col = 0;
  // we don't use size pivoting because the field might not have a valuation
  // !! TODO - add it in this case
  size_t row_max;
  R max_val = echelon._base->zero();
  
  while ((pivot_row < echelon.nrows()) && (pivot_col < echelon.ncols())) {
    for (row_max = pivot_row; row_max < echelon.nrows();) {
      max_val = echelon(row_max, pivot_col);
      if (max_val.isZero())
	row_max++;
      else
	break;
    }
    if (max_val.isZero()) {
      pivot_col++;
    }
    else {
      echelon.swapRows(pivot_row, row_max);
      trans.swapRows(pivot_row, row_max);
      for (size_t row = pivot_row+1; row < echelon.nrows(); row++) {
	R factor = echelon(row,pivot_col) / echelon(pivot_row, pivot_col);
	echelon(row, pivot_col).makeZero();
	for (size_t col = pivot_col + 1; col < echelon.ncols(); col++) {
	  echelon(row,col) -= factor * echelon(pivot_row, col);
	}
	for (size_t col = 0; col < trans.ncols(); col++) {
	  trans(row,col) -= factor * trans(pivot_row, col);
	}
      }
      
      pivot_row++;
      pivot_col++;
    }
  }
  return pivot_row;
}

template<class R, class Parent>
inline size_t Matrix<R,Parent>::rank(void) const
{  
  Matrix<R,Parent> echelon((*this));
  Matrix<R,Parent> trans(_base, echelon.nrows(), echelon.nrows());
  return Matrix<R,Parent>::rowEchelon(echelon, trans);
}

template<class R, class Parent>
inline Matrix<R,Parent> Matrix<R,Parent>::kernel(void) const {
  return this->transpose().leftKernel();
}

template<class R, class Parent>
inline Matrix<R,Parent> Matrix<R,Parent>::leftKernel(void) const {
  Matrix<R,Parent> echelon((*this));
  Matrix<R,Parent> trans(_base, _nrows, _nrows);
  size_t rank = Matrix<R,Parent>::rowEchelon(echelon, trans);
  // getting the zero rows
  Matrix<R,Parent> kernel(_base, _nrows - rank, _nrows);
   for (size_t row = rank; row < _nrows; row++)
    for (size_t col = 0; col < _nrows; col++)
      kernel(row-rank,col) = trans(row, col);
  return kernel;
}

template<class R, class Parent>
inline R Matrix<R,Parent>::trace(void) const
{
  // can only take trace of a square matrix
  assert(this->nrows() == this->ncols());

  R ret = _base->zero();
  for (size_t i = 0; i < this->nrows(); i++)
    ret += (*this)(i,i);
  
  return ret;
}

// We implement Faddeev-LeVerrier here, as this is
// not expected to be a bottleneck

template<class R, class Parent>
inline UnivariatePolyInt<Z> Matrix<R,Parent>::charPoly(void) const
{
  // can only compute characteristic polynomial for a square matrix
  assert(this->nrows() == this->ncols());

  size_t n = this->nrows();
  std::vector< Integer<Z> > c(n+1);
  MatrixInt<Z> z(n,n);
  MatrixInt<Z> I = MatrixInt<Z>::identity(this->_base,n);
  std::vector< MatrixInt<Z> > M(n+1, z);
  MatrixInt<Z> A(n,n);
  for (size_t row = 0; row < n; row++)
    for (size_t col = 0; col < n; col++)
      A(row,col) = birch_util::convertInteger<R,Z>((*this)(row,col));

  c[n] = this->_base->one();
  for (size_t k = 1; k <= n; k++) {
    M[k] = A*M[k-1]+c[n-k+1]*I;
    Z k_Z = k;
    c[n-k] = - (A*M[k]).trace() / k_Z;
  }

  UnivariatePolyInt<Z> p(c);
  return p;
}

template<class R, class Parent>
inline Matrix<R,Parent> Matrix<R,Parent>::restrict(const Matrix<R,Parent> & basis) const
{
  assert(basis.nrows() == basis.rank());
  
  Matrix<R,Parent> echelon = basis;
  Matrix<R,Parent> trans(_base, echelon.nrows(), echelon.nrows());
  rowEchelon(echelon, trans);

  return basis * (*this) * trans.transpose();
}

template<class R, class Parent>
inline Matrix<R,Parent> Matrix<R,Parent>::diagonalJoin(const std::vector< Matrix<R,Parent> > & mats)
{
  size_t nrows = 0;
  size_t ncols = 0;
  for (Matrix<R,Parent> mat : mats) {
    nrows += mat.nrows();
    ncols += mat.ncols();
  }
  assert(mats.size() != 0);
  Matrix<R,Parent> diag(mats[0]._base, nrows, ncols);
  size_t big_row = 0;
  size_t big_col = 0;
  for (Matrix<R,Parent> mat : mats) {
    for (size_t row = 0; row < mat.nrows(); row++)
      for (size_t col = 0; col < mat.ncols(); col++)
	diag(big_row + row, big_col + col) = mat(row, col);
    big_row += mat.nrows();
    big_col += mat.ncols();
  }
  return diag;
}

template<class R, class Parent>
inline Matrix<R,Parent> Matrix<R,Parent>::identity(std::shared_ptr< const Parent > ring, size_t n) {
  Matrix<R,Parent> id(ring, n,n);
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      id(i,j) = (i == j) ? ring->one() : ring->zero();
  return id;
}

// TODO - just change access resolution to the same vector instead
template<class R, class Parent>
inline Matrix<R,Parent> Matrix<R,Parent>::transpose(void) const
{
  std::vector<R> _datat(_nrows*_ncols, _base->zero());
  size_t idx = 0;
  for (size_t col = 0; col < _ncols; col++)
    for (size_t row = 0; row < _nrows; row++)
      _datat[idx++] = (*this)(row, col);
  Matrix<R,Parent> tr(_datat, _ncols, _nrows);
  return tr;
}

template<class R, class Parent>
inline Matrix<R,Parent> Matrix<R,Parent>::operator*(const Matrix<R,Parent> & other) const
{
  size_t nrows = this->_nrows;
  size_t ncols = other._ncols;
  assert( this->_ncols == other._nrows );
  Matrix<R,Parent> prod(_base, nrows, ncols);
    
  for (size_t row = 0; row < nrows; row++)
    for (size_t col = 0; col < ncols; col++) {
      prod(row,col) = _base->zero();
      for (size_t j = 0; j < this->_ncols; j++)
	prod(row,col) += (*this)(row,j)*other(j,col);
    }
  return prod;
}

template<class R, class Parent>
inline Matrix<R,Parent>& Matrix<R,Parent>::operator+=(const Matrix<R,Parent> & other)
{
  for (size_t row = 0; row < this->nrows(); row++)
    for (size_t col = 0; col < this->ncols(); col++)
      (*this)(row, col) += other(row,col);
  
  return (*this);
}

template<class R, class Parent>
inline Matrix<R,Parent>& Matrix<R,Parent>::operator-=(const Matrix<R,Parent> & other)
{
  for (size_t row = 0; row < this->nrows(); row++)
    for (size_t col = 0; col < this->ncols(); col++)
      (*this)(row, col) -= other(row,col);
  
  return (*this);
}

template<class R, class Parent>
inline Matrix<R,Parent>& Matrix<R,Parent>::operator*=(const Matrix<R,Parent> & other)
{
  (*this) = (*this)*other;
  
  return (*this);
}

template<class R, class Parent>
inline Matrix<R,Parent> Matrix<R,Parent>::operator+(const Matrix<R,Parent> & other) const
{
  Matrix<R,Parent> sum(this->_base, this->nrows(), this->ncols());
  for (size_t row = 0; row < this->nrows(); row++)
    for (size_t col = 0; col < this->ncols(); col++)
      sum(row, col) = (*this)(row,col)+other(row,col);
  
  return sum;
}

template<class R, class Parent>
inline Matrix<R,Parent> Matrix<R,Parent>::operator*(const R & a) const
{
  Matrix<R,Parent> prod(this->_base, this->nrows(), this->ncols());
  for (size_t row = 0; row < this->nrows(); row++)
    for (size_t col = 0; col < this->ncols(); col++)
      prod(row,col) = a*((*this)(row,col));

  return prod;
}

template<class R, class Parent>
inline bool Matrix<R,Parent>::isZero(void) const
{
  for (size_t i = 0; i < _data.size(); i++)
    if (!_data[i].isZero())
      return false;
  return true;
}

template<class R, class Parent>
inline bool Matrix<R,Parent>::isOne(void) const
{
  if (_nrows != _ncols) return false;
  for (size_t row = 0; row < _nrows; row++)
    for (size_t col = 0; col < _ncols; col++) {
      if (col == row) {
	if (!((*this)(row,col).isOne())) {
	  return false;
	}
      }
      else {
	if (!((*this)(row,col).isZero())) {
	  return false;
	}
      }
	
    }
  
  return true;
}

template<class R, class Parent>
inline Matrix<R,Parent>& Matrix<R,Parent>::makeZero(void)
{
  for (size_t i = 0; i < _data.size(); i++)
    _data[i].makeZero();

  return (*this);
}
