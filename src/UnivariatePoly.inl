
#include "birch_util.h"
#include "FpElement.h"
#include "Integer.h"

// create the constant polynomial
template<class R, class Parent>
UnivariatePoly<R,Parent>::UnivariatePoly(const R & a)
: _base(a.parent())			     
{
  if (!a.isZero()) {
    this->_coeffs.clear();
    this->_coeffs.push_back(a);
  }
}

// create polynomial from coefficients
template<class R, class Parent>
UnivariatePoly<R,Parent>::UnivariatePoly(std::shared_ptr<const Parent> base_ring, const std::vector<R> & vec)
  : _base(base_ring), _coeffs(vec)
{}

// create the polynomial x^i
template<class R, class Parent>
inline UnivariatePoly<R,Parent> UnivariatePoly<R,Parent>::x(std::shared_ptr<const Parent> base_ring, size_t i)
{
  UnivariatePoly<R,Parent> p(base_ring);
  p._coeffs.resize(i+1);
  for (size_t j = 0; j < i; j++)
    p._coeffs[j] = base_ring->zero();
  p._coeffs[i] = base_ring->one();
  return p;
}

// coefficient of x^i
template<class R, class Parent>
inline R UnivariatePoly<R,Parent>::coefficient(size_t i) const
{
  if (i < this->_coeffs.size())
    return this->_coeffs[i];
  return this->_base->zero();
}

template<class R, class Parent>
inline R UnivariatePoly<R,Parent>::content(void) const
{
  R c = this->_base->zero();
  for (size_t i = 0; i < this->_coeffs.size(); i++)
    c = c.gcd(this->_coeffs[i]);

  return c;
}

template<class R, class Parent>
inline UnivariatePoly<R,Parent> & UnivariatePoly<R,Parent>::operator=(const R & a)
{
  this->_coeffs.clear();
  if (!a.isZero())
    this->_coeffs.push_back(a);

  return (*this);
}
  
// arithmetic
template<class R, class Parent>
inline UnivariatePoly<R,Parent> UnivariatePoly<R,Parent>::operator-() const
{
  UnivariatePoly<R,Parent> neg(this->_base);
  neg._coeffs.resize(this->_coeffs.size());
  for (size_t i = 0; i < this->_coeffs.size(); i++)
    neg._coeffs[i] = -this->_coeffs[i];
  return neg;
}

template<class R, class Parent>
inline void UnivariatePoly<R,Parent>::_eliminateDeg()
{
  // eliminate redundant zeros
  
  size_t i = this->_coeffs.size();
  while((i > 0) && (this->_coeffs[i-1].isZero())) i--;
  this->_coeffs.resize(i);

  return;
}

template<class R, class Parent>
inline UnivariatePoly<R,Parent>
UnivariatePoly<R,Parent>::operator+(const UnivariatePoly<R,Parent> & other) const
{
  UnivariatePoly<R,Parent> sum(this->_base);
  if (this->_coeffs.size() < other._coeffs.size())
    return other + (*this);
  // here we may assume this is the polynomial with the larger degree
  sum._coeffs.resize(this->_coeffs.size());
  size_t i;
  for (i = 0; i < other._coeffs.size(); i++)
    sum._coeffs[i] = this->_coeffs[i] + other._coeffs[i];
  for (; i < this->_coeffs.size(); i++)
    sum._coeffs[i] = this->_coeffs[i];

  if (this->_coeffs.size() == other._coeffs.size())
    sum._eliminateDeg();
  
  return sum;
}

template<class R, class Parent>
inline UnivariatePoly<R,Parent>
UnivariatePoly<R,Parent>::operator-(const UnivariatePoly<R,Parent> & other) const
{
  return (*this) + (-other);
}

template<class R, class Parent>
inline UnivariatePoly<R,Parent>
UnivariatePoly<R,Parent>::operator*(const UnivariatePoly<R,Parent> & other) const
{
  if (this->isZero())
    return (*this);
  if (other.isZero())
    return other;
  UnivariatePoly<R,Parent> prod(this->_base);
  prod._coeffs.resize(this->degree()+other.degree()+1);
  std::fill(prod._coeffs.begin(), prod._coeffs.end(), this->_base->zero());
  size_t i, j;
  for (i = 0; i < this->_coeffs.size(); i++)
    for (j = 0; j < other._coeffs.size(); j++)
      prod._coeffs[i+j] += this->_coeffs[i] * other._coeffs[j];

  // only needed if R is not an integral domain - we assume R is Euclidean
  // prod._eliminateDeg();
  
  return prod;
}

template<class R, class Parent>
inline UnivariatePoly<R,Parent>
UnivariatePoly<R,Parent>::operator/(const UnivariatePoly<R,Parent> & other) const
{
  UnivariatePoly<R,Parent> q(this->_base);
  UnivariatePoly<R,Parent> r(this->_base);
  
  divRem((*this),other,q,r);

  return q;
}

template<class R, class Parent>
inline UnivariatePoly<R,Parent>
UnivariatePoly<R,Parent>::operator%(const UnivariatePoly<R,Parent> & other) const
{
  UnivariatePoly<R,Parent> q(this->_base);
  UnivariatePoly<R,Parent> r(this->_base);
  
  divRem((*this),other,q,r);

  return r;
}

template<class R, class Parent>
inline UnivariatePoly<R,Parent> UnivariatePoly<R,Parent>::operator*(const R & a) const
{
  if (a.isZero())
    return UnivariatePoly<R,Parent>::zero(this->_base);
  
  UnivariatePoly<R,Parent> prod(this->_base);
  prod._coeffs.resize(this->_coeffs.size());
  for (size_t i = 0; i < this->_coeffs.size(); i++)
    prod._coeffs[i] = a * this->_coeffs[i];

  // only relevant if a == 0
  //   prod._eliminateDeg();
  
  return prod;
}

template<class R, class Parent>
inline UnivariatePoly<R,Parent> UnivariatePoly<R,Parent>::operator/(const R & a) const
{
  UnivariatePoly<R,Parent> prod(this->_base);
  prod._coeffs.resize(this->_coeffs.size());
  for (size_t i = 0; i < this->_coeffs.size(); i++)
    prod._coeffs[i] = this->_coeffs[i] / a;
  
  return prod;
}

template<class R, class Parent>
inline UnivariatePoly<R,Parent> UnivariatePoly<R,Parent>::operator%(const R & a) const
{
  UnivariatePoly<R,Parent> res(this->_base);
  res._coeffs.resize(this->_coeffs.size());
  for (size_t i = 0; i < this->_coeffs.size(); i++)
    res._coeffs[i] = this->_coeffs[i] % a;

  res._eliminateDeg();
  
  return res;
}

template<class R, class Parent>
inline UnivariatePoly<R,Parent> &
UnivariatePoly<R,Parent>::operator+=(const UnivariatePoly<R,Parent> & other)
{
  size_t i, deg;
  deg = this->_coeffs.size();
  if (deg < other._coeffs.size()) {
    this->_coeffs.resize(other._coeffs.size());
    for (i = 0; i < deg; i++)
      this->_coeffs[i] += other._coeffs[i];
    for (; i < other._coeffs.size(); i++)
      this->_coeffs[i] = other._coeffs[i];
  }
  else {
    for (i = 0; i < other._coeffs.size(); i++)
      this->_coeffs[i] += other._coeffs[i];
     // eliminate redundant zeros
    if (this->_coeffs.size() == other._coeffs.size()) {
      this->_eliminateDeg();
    }
  }
  return (*this);
}

template<class R, class Parent>
inline UnivariatePoly<R,Parent> &
UnivariatePoly<R,Parent>::operator-=(const UnivariatePoly<R,Parent> & other)
{
  return ((*this) += (-other));
}

template<class R, class Parent>
inline UnivariatePoly<R,Parent> &
UnivariatePoly<R,Parent>::operator*=(const UnivariatePoly<R,Parent> & other)
{
  (*this) = (*this)*other;
  return (*this);
}

template<class R, class Parent>
inline UnivariatePoly<R,Parent> &
UnivariatePoly<R,Parent>::operator/=(const UnivariatePoly<R,Parent> & other)
{
  (*this) = (*this)/other;
  return (*this);
}

template<class R, class Parent>
inline UnivariatePoly<R,Parent> &
UnivariatePoly<R,Parent>::operator%=(const UnivariatePoly<R,Parent> & other)
{
  (*this) = (*this)%other;
  return (*this);
}

template<class R, class Parent>
inline UnivariatePoly<R,Parent>& UnivariatePoly<R,Parent>::operator*=(const R & a)
{
  if (a.isZero())
    return makeZero();
  
  for (size_t i = 0; i < this->_coeffs.size(); i++)
    this->_coeffs[i] *= a;

  // This is only relevant if a == 0
  // this->_eliminateDeg();
  
  return (*this);
}

template<class R, class Parent>
inline UnivariatePoly<R,Parent>& UnivariatePoly<R,Parent>::operator/=(const R & a)
{
  for (size_t i = 0; i < this->_coeffs.size(); i++)
    this->_coeffs[i] /= a;
  
  return (*this);
}

template<class R, class Parent>
inline UnivariatePoly<R,Parent>& UnivariatePoly<R,Parent>::operator%=(const R & a)
{
  for (size_t i = 0; i < this->_coeffs.size(); i++)
    this->_coeffs[i] %= a;

  this->_eliminateDeg();
  
  return (*this);
}

template<class R, class Parent>
inline UnivariatePoly<R,Parent>
UnivariatePoly<R,Parent>::evaluate(const UnivariatePoly<R,Parent> & f) const
{
  UnivariatePoly<R,Parent> comp(this->_base);
  UnivariatePoly<R,Parent> f_i(this->_base->one());
  for (size_t i = 0; i < this->_coeffs.size(); i++) {
    comp += this->_coeffs[i]*f_i;
    f_i *= f;
  }
  return comp;
}

template<class R, class Parent>
inline R UnivariatePoly<R,Parent>::evaluate(const R & a) const
{
  R res;
  R a_i = this->_base->one();
  for (size_t i = 0; i < this->_coeffs.size(); i++) {
    res += this->_coeffs[i]*a_i;
    a_i *= a;
  }
  return res;
}

template<class R, class Parent>
template<class S, class SParent>
inline Matrix<S,SParent> UnivariatePoly<R,Parent>::evaluate(const Matrix<S,SParent> & a) const
{
  assert(a.nrows() == a.ncols());

  Matrix<S,SParent> res(a.baseRing(), a.nrows(), a.nrows());
  Matrix<S,SParent> a_i = Matrix<S,SParent>::identity(a.baseRing(), a.nrows());
  for (size_t i = 0; i < this->_coeffs.size(); i++) {
    res += birch_util::convert<R,S>(this->_coeffs[i])*a_i;
    a_i *= a;
  }
  return res;
}

// booleans
template<class R, class Parent>
inline bool UnivariatePoly<R,Parent>::operator==(const UnivariatePoly<R,Parent> & other) const
{
  if (this->_coeffs.size() != other._coeffs.size())
    return false;
  for (size_t i = 0; i < this->_coeffs.size(); i++)
    if(this->_coeffs[i] != other._coeffs[i])
      return false;
  
  return true;
}

template<class R, class Parent>
inline bool UnivariatePoly<R,Parent>::operator!=(const UnivariatePoly<R,Parent> & other) const
{
  return !((*this)==other);
}

template<class R, class Parent>
inline bool UnivariatePoly<R,Parent>::operator==(const R & a) const
{
  if (a.isZero()) return this->isZero();
  if (this->_coeffs.size() != 1)
    return false;

  return (this->_coeffs[0] == a);
}

template<class R, class Parent>
inline bool UnivariatePoly<R,Parent>::operator!=(const R & a) const
{
  return !((*this)==a);
}

template<class R, class Parent>
inline std::ostream& operator<<(std::ostream& os, const UnivariatePoly<R,Parent> & p)
{
  p.print(os);
  return os;
}

template<class R, class Parent>
inline UnivariatePoly<R,Parent> UnivariatePoly<R,Parent>::derivative() const
{
  UnivariatePoly<R,Parent> f_prime(this->_base);
  f_prime._coeffs.resize(this->degree());
  R a = this->_base->one();
  for (size_t i = 1; i < this->_coeffs.size(); i++,a++)
    f_prime._coeffs[i-1] = a * this->_coeffs[i];
  
  return f_prime;
}

template<class R, class Parent>
inline void UnivariatePoly<R,Parent>::divRem(const UnivariatePoly<R,Parent> & f,
					     const UnivariatePoly<R,Parent> & g,
					     UnivariatePoly<R,Parent> & q,
					     UnivariatePoly<R,Parent> & r)
{
  assert(!g.isZero());

  UnivariatePoly<R,Parent> t(f.baseRing());
  
  q = f.baseRing()->zero();
  r = f;

  // we will use pseudo-remainder
  if (f.degree() >= g.degree()) {
    size_t d = f.degree() + 1 - g.degree();
    R g_d = g.lead()^d;
    r *= g_d;
  }

  while ((!r.isZero()) && (r.degree() >= g.degree())) {
    R lc = r.lead() / g.lead();
    t = lc * UnivariatePoly<R,Parent>::x(f.baseRing(),r.degree()-g.degree());
    q += t;
    r -= t*g;
  }

  return;
}

template<class R, class Parent>
inline UnivariatePoly<R,Parent> UnivariatePoly<R,Parent>::gcd(const UnivariatePoly<R,Parent> & f,
							      const UnivariatePoly<R,Parent> & g)
{
  UnivariatePoly<R,Parent> q(f.baseRing());
  UnivariatePoly<R,Parent> r_minus(f.baseRing());
  UnivariatePoly<R,Parent> r(f.baseRing());
  UnivariatePoly<R,Parent> r_plus(f.baseRing());

  r_minus = f / f.content();
  r = g / g.content();
  
  while (!r.isZero()) {
     UnivariatePoly<R,Parent>::divRem(r_minus, r, q, r_plus);
    r_minus = r;
    r = r_plus / r_plus.content();
  }
  
  return r_minus;
}

template<class R, class Parent>
inline UnivariatePoly<R,Parent> UnivariatePoly<R,Parent>::xgcd(const UnivariatePoly<R,Parent> & f,
							       const UnivariatePoly<R,Parent> & g,
							       UnivariatePoly<R,Parent> & s,
							       UnivariatePoly<R,Parent> & t)
{
  UnivariatePoly<R,Parent> q(f.baseRing());
  UnivariatePoly<R,Parent> r_minus(f.baseRing());
  UnivariatePoly<R,Parent> s_minus(f.baseRing());
  UnivariatePoly<R,Parent> t_minus(f.baseRing());
  UnivariatePoly<R,Parent> r(f.baseRing());
  UnivariatePoly<R,Parent> r_plus(f.baseRing());
  UnivariatePoly<R,Parent> s_plus(f.baseRing());
  UnivariatePoly<R,Parent> t_plus(f.baseRing());
  
  s = f.baseRing()->one();
  s_minus = f.baseRing()->zero();
  t = f.baseRing()->zero();
  t_minus = f.baseRing()->one();
  
  r_minus = f / f.content();
  r = g / g.content();
  
  while (r != f.baseRing()->zero()) {
    divRem(r_minus, r, q, r_plus);
    
    R c = r_plus.content();
    W64 e = r_minus.degree()+1-r.degree();
    R a = r.lead()^e;
    
    r_minus = r;
    r = r_plus / c;
    s_plus = (a*s_minus - q*s) / c;
    t_plus = (a*t_minus - q*t) / c;
    s_minus = s;
    t_minus = t;
    s = s_plus;
    t = t_plus;
  }

  // finalize
  s = s_minus;
  t = t_minus;
  
  return r_minus;
}

template<class R, class Parent>
inline UnivariatePoly<R,Parent> &  UnivariatePoly<R,Parent>::makeZero(void)
{
  this->_coeffs.clear();
  return *this;
}

template<class R, class Parent>
inline UnivariatePoly<R,Parent> & UnivariatePoly<R,Parent>::makeOne(void)
{
  this->_coeffs.resize(1);
  this->_coeffs[0] = this->_base->one();
  return *this;
}

template<class R, class Parent>
void UnivariatePoly<R,Parent>::print(std::ostream & os) const
{
  const UnivariatePoly<R,Parent> & p = *this;
  size_t deg = p.degree();
  for (size_t i = deg+1; i > 0; i--) {
    R coeff = p.coefficient(i-1);
    if (!coeff.isZero()) {
      if ((i <= deg) && (coeff > p.baseRing()->zero()))
	os << '+';
      else if ((coeff == -p.baseRing()->one()) && (i != 1))
	os << '-';
      else if (coeff != p.baseRing()->one())
	os << coeff;
      if (i > 1)
	os << "x";
      if (i > 2)
	os << "^" << (i-1);
    }
  }
  return;
}

