template <typename R>
inline Rational<R> Rational<R>::operator+(const Rational<R> &other) const
{
  Rational<R> sum;
  sum._num = (this->_num) * other._denom + (this->_denom) * other._num;
  sum._denom = (this->_denom) * other._denom;
  sum._reduce();
  return sum;
}

template <typename R>
inline Rational<R> Rational<R>::operator*(const Rational<R> &other) const
{
  Rational<R> prod;
  prod._num = (this->_num) * other._num;
  prod._denom = (this->_denom) * other._denom;
  prod._reduce();
  return prod;
}

template <typename R>
inline Rational<R> Rational<R>::operator/(const Rational<R> &other) const
{
  Rational<R> prod;
  prod._num = (this->_num) * other._denom;
  prod._denom = (this->_denom) * other._num;
  prod._reduce();
  return prod;
}

template <typename R>
inline bool Rational<R>::operator==(const Rational<R> &other) const
{
  return (this->_num) * other._denom == other._num * (this->_denom); 
}

template <typename R>
inline bool Rational<R>::operator<(const Rational<R> &other) const
{
  Rational<R> diff = (*this)-other;
  return diff.num() * diff.denom() < 0; 
}

template<typename R>
inline void Rational<R>::_reduce(void)
{
  Integer<R> d = _num.gcd(_denom);
  _num /= d;
  _denom /= d;
  return;
}

// collect all of the along the way.
// We are using Akiyama and Tanigawa's algorithm
// It's not the fastest, but it is one of the simplest.
template <typename R>
inline std::vector< Rational<R> > Rational<R>::_bernoulliUpTo(size_t n)
{
  std::shared_ptr<const RationalField<R> > QQ = std::make_shared<const RationalField<R> >();
  std::vector< Rational<R> > a(n+1, QQ);
  std::vector< Rational<R> > b(n+1, QQ);
  for (size_t i = 0; i <= n; i++) {
    Rational<R> r(1, i+1);
    a[i] = r;
  }
  b[0] = a[0];
  for (size_t i = 0; i < n; i++)
    {
      for (size_t j = 0; j < n - i; j++) {
	Integer<R> mult = R(j + 1);
	a[j] = mult*(a[j] - a[j+1]);
      }
      b[i+1] = a[0];
    }
  
  return b;
  
}

template <typename R>
inline Rational<R> Rational<R>::bernoulliNumber(size_t n)
{
  std::vector< Rational<R> > b = _bernoulliUpTo(n);
  return b[n];
}


template <typename R>
inline std::vector< Rational<R> > Rational<R>::_bernoulliPoly(size_t n)
{
  std::vector< Rational<R> > b = _bernoulliUpTo(n);
  std::reverse(b.begin(), b.end());
  for (size_t k = 0; k <= n; k++) {
    Integer<R> n_int = n;
    b[k] *= n_int.binomialCoefficient(R(k));
  }
  return b;
}

// B_{n. chi} where chi is the quadratic character corresponding to
// quadratic field with discrminant d (Is it? verify we are working
// with the correct discriminant (d or 4d maybe?))
template <typename R>
inline Rational<R> Rational<R>::bernoulliNumber(size_t n, const Integer<R> & d)
{
  std::vector< Rational<R> > b = Rational<R>::_bernoulliPoly(n);
  Integer<R> d_pow = d^n;

  Rational<R> b_chi = Rational<R>::zero();
  for (Integer<R> a = Integer<R>::zero(); a < d; a++)
    {
      int chi_a = a.kroneckerSymbol(n);      
      Integer<R> a_pow = Integer<R>::one();
      Rational<R> s = Rational<R>::zero();
      for (size_t k = 0; k <= n; k++)
	{
	  d_pow /= d;
	  s += b[k]*a_pow*d_pow;
	  a_pow *= a;
	}
      s *= R(chi_a);
      b_chi += s;
    }
  return b_chi;
}
