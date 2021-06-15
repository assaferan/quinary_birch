// Implementation of templated functions from QuadForm.h

// c-tors
template<class R, class Parent, size_t n>
QuadForm<R,Parent,n>::QuadForm(const typename QuadForm<R,Parent,n>::SymVec& coeffs)
  : _B(coeffs[0].parent())
{
  size_t idx = 0;
  for (size_t row = 0; row < n; row++)  {
    for (size_t col = 0; col < row; col++)	{	
      this->_B(row,col) = coeffs[idx];
      this->_B(col,row) = coeffs[idx++];
    }
    this->_B(row,row) = coeffs[idx++];
  }
}

// assignment
template<class R, class Parent, size_t n>
inline QuadForm<R,Parent,n>&
QuadForm<R,Parent,n>::operator=(const QuadForm<R,Parent,n> & other)
{
  if (this != &other) {
    this->_B = other._B;
  }
  return *this;
}

// When n is odd we return the half-discriminant
template<class R, class Parent, size_t n>
inline R QuadForm<R,Parent,n>::discriminant(void) const
{
  R det = this->_B.determinant();
  R one = _B.baseRing()->one();
  R two = one+one;
  // !! TODO - add discriminant in the characteristic 2 case
  assert(!two.isZero());
  return (n % 2 == 0) ? det : det/two;
}

template<class R, class Parent, size_t n>
inline std::ostream& operator<<(std::ostream& os, const QuadForm<R,Parent,n> & q)
{
  os << q.bilinearForm();
  return os;
}
