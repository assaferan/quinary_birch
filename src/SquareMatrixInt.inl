template<typename R, size_t n>
SquareMatrixInt<R,n>::SquareMatrixInt(const R mat[n][n])
  : SquareMatrix< Integer<R>, IntegerRing<R>, n>(mat)
{}

template <typename R, size_t n>
SquareMatrixInt<R,n>::SquareMatrixInt(const SquareMatrixInt<R,n> & mat)
  : SquareMatrix< Integer<R>, IntegerRing<R>, n>(mat)
{}

template<typename R, size_t n>
SquareMatrixInt SquareMatrixInt<R,n>::operator*(const SquareMatrixInt & other) const
{
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
