// implementation for Isometry.h

template<typename R, size_t n>
inline SquareMatrixInt<R,n>
Isometry<R,n>::transform(const SquareMatrixInt<R,n>& from) const
{
  return (this->_a).transpose()*from*(this->_a) / (this->_scale * this->_scale);
}

template<typename R, size_t n>
inline bool Isometry<R,n>::isIsometry(const QuadFormInt<R,n>& from,
				      const QuadFormInt<R,n>& to) const
{
  R scalar = this->_scale * this->_scale;
  Integer<R> val;
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++) {
      val = 0;
      for (size_t k = 0; k < n; k++)
	for (size_t l = 0; l < n; l++)
	  val += this->_a(k,i)*from.bilinearForm()(k,l)*
	    this->_a(l,j);
      if (val != to.bilinearForm()(i,j) * scalar)
	return false;
    }
  return true;
}

template<typename R, size_t n>
inline void Isometry<R,n>::updatePerm(const VectorInt<size_t, n> & perm) {
  SquareMatrixInt<R,n> temp;
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      temp(i,j) = this->_a(i,j);
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      this->_a(perm[i],j) = temp(i,j);
  return;
}

template<typename R, size_t n>
inline Isometry<R,n> Isometry<R,n>::inverse(void) const
{
  std::shared_ptr< const RationalField<R> >
    QQ = std::make_shared< const RationalField<R> >();
  // !! TODO - should be able to invert without using rationals
  // for example, can always track back (save the inverse for the ride)
  SquareMatrixRat<R,n> a_rat(QQ);
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++) {
      a_rat(i,j) = this->_a(i,j);
      a_rat(i,j) /= this->_scale;
    }
  a_rat = a_rat.inverse();
  
  SquareMatrixInt<R,n> a_inv;
  // Since this is an isometry, the inverse should be integral
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      a_inv(i,j) = (this->_scale * a_rat(i,j)).floor().num();
#ifdef DEBUG
  R scale2 = (this->_scale)*(this->_scale);
  assert((a_inv * (this->_a) == scale2*SquareMatrixInt<R,n>::identity()));
#endif
  return Isometry(a_inv, this->_scale);
}

template<typename R, size_t n>
inline void Isometry<R,n>::rescale(void)
{
  Integer<R> d = this->_scale;
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      d = d.gcd(this->_a(i,j));
  
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      this->_a(i,j) /= d.num();

  this->_scale /= d.num();
  return;
}
