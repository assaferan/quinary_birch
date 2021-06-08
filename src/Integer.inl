// implementation of methods for Integer<R>

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
Integer<R>::euclideanDivision(const Integer<R>& b) const
{
  Integer<R> q(this->_num / b._num);
  Integer<R> r(this->_num % b._num);
  
  //return std::pair< Integer<R>, Integer<R> >(q,r);
  return std::make_pair(q,r);
}
