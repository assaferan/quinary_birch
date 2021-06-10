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
inline typename Integer<R>::FactorData Integer<R>::factorization() const
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
