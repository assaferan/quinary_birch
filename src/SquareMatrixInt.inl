template<typename R, size_t n>
SquareMatrixInt<R,n>::SquareMatrixInt(const R mat[n][n])
  : SquareMatrix< Integer<R>, IntegerRing<R>, n>(mat)
{}

template <typename R, size_t n>
SquareMatrixInt<R,n>::SquareMatrixInt(const SquareMatrixInt<R,n> & mat)
  : SquareMatrix< Integer<R>, IntegerRing<R>, n>(mat)
{}

template<typename R, size_t n>
SquareMatrixInt<R,n> SquareMatrixInt<R,n>::operator*(const SquareMatrixInt & other) const
{
  // Here we avoid construction of multiple Integer<R> objects
  R sum;
  SquareMatrixInt<R,n> prod(this->baseRing());
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++) {
      sum = 0;
      for (size_t k = 0; k < n; k++)
	sum += (*this)(i,k).num()*other(k,j).num();
      prod(i,j) = sum;
    }
  return prod;
}

// For now we just use the transpose from the base class and our conversion
template<typename R, size_t n>
SquareMatrixInt<R,n> SquareMatrixInt<R,n>::transpose(void) const
{
  return SquareMatrix< Integer<R>, IntegerRing<R>, n>::transpose();
}

template<typename R, size_t n>
SquareMatrixInt<R,n>  SquareMatrixInt<R,n>::operator/(const Integer<R> & scalar) const
{
  assert (!scalar.isZero());
  
  SquareMatrixInt<R,n> quo(this->baseRing());
  for (size_t row = 0; row < n; row++)
    for (size_t col = 0; col < n; col++)
      quo(row,col) = R((*this)(row,col).num() / scalar.num());
  return quo;
}

