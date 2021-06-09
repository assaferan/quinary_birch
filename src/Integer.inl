// implementation of methods for Integer<R>

#include <cassert>

// assignment
template<typename R>
Integer<R> & Integer<R>::operator=(const Integer<R> & b)
{
  if (this != &b) {
    _num = b._num;
  }
  return (*this);
}

// arithmetic

template <typename R>
Integer<R> & Integer<R>::operator+=(const Integer<R> &other)
{
  this->_num += other._num;
  return (*this);
}

template <typename R>
Integer<R> & Integer<R>::operator-=(const Integer<R> &other)
{
  this->_num -= other._num;
  return (*this);
}

template <typename R>
Integer<R> & Integer<R>::operator*=(const Integer<R> &other)
{
  this->_num *= other._num;
 
  return (*this);
}

template <typename R>
bool Integer<R>::operator==(const Integer<R> &other) const
{
  return (this->_num == other._num); 
}

template <typename R>
bool Integer<R>::operator<(const Integer<R> &other) const
{
  return (this->_num < other._num); 
}

template <typename R>
void Integer<R>::print(std::ostream& os) const
{
  os << this->_num;
}

// Euclidean division
template <typename R>
typename EuclideanDomain<Integer<R> >::DivRes
Integer<R>::euclideanDivision(const Integer<R>& other) const
{
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
  
  assert((r_int._num >= 0) && (r_int._num < abs(other._num)));
  assert(*this == q_int*other+r_int);
	 
  return std::make_pair(q_int,r_int);
}
