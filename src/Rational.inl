template <typename R>
Rational<R> Rational<R>::operator+(const Rational<R> &other) const
{
  Rational<R> sum;
  sum._num = (this->_num) * other._denom + (this->_denom) * other._num;
  sum._denom = (this->_denom) * other._denom;
  sum.reduce();
  return sum;
}

template <typename R>
Rational<R> Rational<R>::operator*(const Rational<R> &other) const
{
  Rational<R> prod;
  prod._num = (this->_num) * other._num;
  prod._denom = (this->_denom) * other._denom;
  prod.reduce();
  return prod;
}

template <typename R>
Rational<R> Rational<R>::operator/(const Rational<R> &other) const
{
  Rational<R> prod;
  prod._num = (this->_num) * other._denom;
  prod._denom = (this->_denom) * other._num;
  prod.reduce();
  return prod;
}

template <typename R>
bool Rational<R>::operator==(const Rational<R> &other) const
{
  return (this->_num) * other._denom == other._num * (this->_denom); 
}

template <typename R>
bool Rational<R>::operator<(const Rational<R> &other) const
{
  Rational<R> diff = (*this)-other;
  return diff.num() * diff.denom() < 0; 
}

template<typename R>
void Rational<R>::reduce(void)
{
  Integer<R> d = _num.gcd(_denom);
  _num /= d;
  _denom /= d;
  return;
}
