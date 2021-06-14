// implementation of methods for Integer<R>

#include <cassert>

// assignment
template<typename R>
inline Integer<R> & Integer<R>::operator=(const Integer<R> & b)
{
  if (this != &b) {
    _num = b._num;
  }
  return (*this);
}

// arithmetic

template <typename R>
inline Integer<R> & Integer<R>::operator+=(const Integer<R> &other)
{
  this->_num += other._num;
  return (*this);
}

template <typename R>
inline Integer<R> & Integer<R>::operator-=(const Integer<R> &other)
{
  this->_num -= other._num;
  return (*this);
}

template <typename R>
inline Integer<R> & Integer<R>::operator*=(const Integer<R> &other)
{
  this->_num *= other._num;
 
  return (*this);
}

template <typename R>
inline bool Integer<R>::operator==(const Integer<R> &other) const
{
  return (this->_num == other._num); 
}

template <typename R>
inline bool Integer<R>::operator<(const Integer<R> &other) const
{
  return (this->_num < other._num); 
}

template <typename R>
inline void Integer<R>::print(std::ostream& os) const
{
  os << this->_num;
}

// Euclidean division
template <typename R>
inline typename EuclideanDomainElement<Integer<R>, IntegerRing<R> >::DivRes
Integer<R>::euclideanDivision(const Integer<R>& other) const
{
  assert(!other.isZero());
  
  R a = this->_num;
  R b = other._num;
  if (a < 0) {
    a = -a-1;
  }
  if (b < 0) {
    b = -b;
  }
  
  R q = a / b;
  R r = a % b;
  
  if (this->_num < 0) {
    q = -q-1;
    r = b-1-r;
  }

  if (other._num < 0) {
    q = -q;
  }

  Integer<R> q_int(q);
  Integer<R> r_int(r);
  
  assert((r_int._num >= 0) && (r_int._num < b));
  assert(*this == q_int*other+r_int);
	 
  return std::make_pair(q_int,r_int);
}

// We do only trial divison, since our numbers are always small enough
// (in all our use cases, num will be in the order of 1000 at most)

template <typename R>
inline typename Integer<R>::FactorData Integer<R>::factorization(void) const
{
  Integer<R> num = *this;
  assert(num < R(10000));
  typename Integer<R>::FactorData factors;
  Integer<R> temp_num = num.abs();
  size_t exp;
  Integer<R> a = R(2);
  while (!temp_num.isOne())
    {
      if ((temp_num % a).isZero())
	{
	  exp = 1; 
	  temp_num /= a;
	  while ((temp_num % a).isZero())
	    {
	      exp++;
	      temp_num /= a;
	    }
	  factors.push_back(std::make_pair(a, exp));
	}
      a++;
    }
  return factors;
}

template <typename R>
inline size_t Integer<R>::valuation(const Integer<R>& p) const
{
  assert(!(this->isZero()));
  Integer<R> t = *this;
  
  size_t exp = 0;

  while ((t % p).isZero())
    {
      exp++;
      t /= p;
    }

  return exp;
}

template <typename R>
inline int Integer<R>::kroneckerSymbol(const Integer<R> & n) const
{
  const Integer<R> & a = *this;
  // extremal cases
  if (n.isZero()) return (this->abs().isOne()) ? 1 : 0;
  if ((-n).isOne()) return (*this < Integer<R>::zero()) ? -1 : 1;
  if (n.isOne()) return 1;
  if (n == R(2)) {
    if (((*this) % R(2)).isZero()) return 0;
    R val = (this->num() % 8)/2;
    if ((val == 0) || (val == 3))
      return 1;
    return -1;
  }
  if ((-a).isOne()) {
    R n_prime = n.num();
    while (n_prime % 2 == 0) n_prime /= 2;
    return ((n_prime / 2) % 2 == 0) ? 1 : -1;
  }
  if ((this->num() == 2) && ((n % R(2)).isOne())) {
    return (((n^2) / R(8)) % R(2)).isZero() ? 1 : -1;
  }
  // multiplicativity
  if (n < Integer<R>::zero()) return this->kroneckerSymbol(-Integer<R>::one())*this->kroneckerSymbol(-n);
  if (a < Integer<R>::zero())
    return (-Integer<R>::one()).kroneckerSymbol(n)*(-a).kroneckerSymbol(n);

  // now may assume n >= 3, a >= 0
 
  // quadratic reciprocity
  if (a < n) {
    Integer<R> n_star;
    R n_prime = n.num();
    while (n_prime % 2 == 0) n_prime /= 2;
    n_star = ((n_prime / 2) % 2 == 0) ? n : -n;
    return n_star.kroneckerSymbol(*this);
  }

  // now we may also assume a ge n

  // if n = 2 mod 4, we can't reduce, use multiplicativity again
  if (n.num() % 4 == 2) return kroneckerSymbol(n/R(2))*kroneckerSymbol(R(2));
  // now we can reduce
  return (a % n).kroneckerSymbol(n);
}


template <typename R>
inline bool Integer<R>::isLocalSquare(const Integer<R>& p) const
{
  if (this->isZero()) return true;
  size_t v = this->valuation(p);
  if (v % 2 == 1) return false;
  Integer<R> a0 = *this;
  for (size_t i = 0; i < v; i++) a0 /= p;
  bool ok = (a0.kroneckerSymbol(p) != -1);
  if (p != R(2)) return ok;
  size_t w = (a0-Integer<R>::one()).valuation(p);
  assert(w >= 1);
  size_t ee = 2;

  while ((w < ee) && (w % 2 == 0)) {
    R ww = (1+ (1<< (w/2)))*(1+ (1<< (w/2)));
    a0 /= ww;
    w = (a0-Integer<R>::one()).valuation(p);
  }
  return ((w > ee) || ((w == ee) && (a0 % R(8)).isOne()));
}

template<typename R>
inline Integer<R> Integer<R>::binomialCoefficient(const Integer<R> & k) const
{
  const Integer<R> & n = *this;
  Integer<R> res = Integer<R>::one();
  if (k > n - k)
    return n.binomialCoefficient(n-k);
  for (Integer<R> i = 0; i < k; i++) {
    res *= (n-i);
    res /= (i+Integer<R>::one());
  }
  return res;
}

template <typename R>
inline Integer<R> Integer<R>::nextPrime(void) const
{
  Z p = birch_util::convert_Integer<R,Z>(this->num());
  mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
  return birch_util::convert_Integer<Z,R>(p); 
}
