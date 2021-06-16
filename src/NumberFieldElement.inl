template<typename R>
inline NumberFieldElement<R> & NumberFieldElement<R>::operator=(const R & a)
{
  this->_elt = a;
  return (*this);
}

template<typename R>
inline NumberFieldElement<R> NumberFieldElement<R>::operator-(void) const
{
  NumberFieldElement<R> neg(this->_K);
  neg._elt = -(this->_elt);
  
  return neg;
}

template<typename R>
inline NumberFieldElement<R>
NumberFieldElement<R>::operator+(const NumberFieldElement<R> & other) const
{
  NumberFieldElement<R> sum(this->_K);
  sum._elt = (this->_elt) + other._elt;
  
  return sum;
}

template<typename R>
inline NumberFieldElement<R>
NumberFieldElement<R>::operator-(const NumberFieldElement<R> & other) const
{
  return (*this)+(-other);
}

template<typename R>
inline NumberFieldElement<R>
NumberFieldElement<R>::operator*(const NumberFieldElement<R> & other) const
{
  NumberFieldElement<R> prod(this->_K);
  prod._elt = (this->_elt)  * other._elt % (this->_K->modulus());

  return prod;
}

template<typename R>
inline NumberFieldElement<R>
NumberFieldElement<R>::operator/(const NumberFieldElement<R> & other) const
{
  return (*this) * other.inverse();
}

template<typename R>
inline NumberFieldElement<R> NumberFieldElement<R>::operator*(const R & a) const
{
  NumberFieldElement<R> prod(this->_K);
  prod._elt = (this->_elt) * a;

  return prod;
}

template<typename R>
inline NumberFieldElement<R> NumberFieldElement<R>::operator/(const R & a) const
{
  NumberFieldElement<R> quo(this->_K);
  quo._elt = (this->_elt) / a;

  return quo;
}

template<typename R>
inline NumberFieldElement<R> &
NumberFieldElement<R>::operator+=(const NumberFieldElement<R> & other)
{
  this->_elt += other._elt;
  return (*this);
}

template<typename R>
inline NumberFieldElement<R> &
NumberFieldElement<R>::operator-=(const NumberFieldElement<R> & other)
{
  this->_elt -= other._elt;
  return (*this);
}

template<typename R>
inline NumberFieldElement<R> &
NumberFieldElement<R>::operator*=(const NumberFieldElement<R> & other)
{
  *this = (*this)*other;
  return (*this);
}

template<typename R>
inline NumberFieldElement<R> &
NumberFieldElement<R>::operator/=(const NumberFieldElement<R> & other)
{
  *this = (*this)/other;
  return (*this);
}

template<typename R>
inline NumberFieldElement<R>& NumberFieldElement<R>::operator*=(const R & a)
{
  this->_elt *= a;
  return (*this);
}

template<typename R>
inline NumberFieldElement<R>& NumberFieldElement<R>::operator/=(const R & a)
{
  this->_elt /= a;
  return (*this);
}

template<typename R>
inline NumberFieldElement<R>
NumberFieldElement<R>::inverse(void) const
{
  std::shared_ptr<const RationalField<R> > QQ = std::make_shared< const RationalField<R> >();
  UnivariatePolyRat<R> s(QQ);
  UnivariatePolyRat<R> t(QQ);
  
  UnivariatePolyRat<R>::xgcd(this->_elt, this->_K->modulus(), s, t);
  
  NumberFieldElement<R> inv(this->_K, s);
  return inv;
}
