// VectorInt

template<typename R, size_t n>
inline bool VectorInt<R,n>::operator==(const VectorInt<R,n> & other) const
{
  for (size_t i = 0; i < n; i++)
    if (this->_v[i] != other[i]) return false;
  return true;
}

// for comparison in standard containers such as std::set and std::map
template<typename R, size_t n>
inline bool VectorInt<R,n>::operator<(const VectorInt<R,n> & other) const
{
  for (size_t i = 0; i < n; i++) {
    if (this->_v[i] > other[i]) return false;
    if (this->_v[i] < other[i]) return true;
  }
  return false;
}

template<typename R, size_t n>
inline VectorInt<R,n> VectorInt<R,n>::operator-() const {
  VectorInt<R,n> neg;
  for (size_t i = 0; i < n; i++)
    neg[i] = -this->_v[i];
  return neg;
}

template<typename R, size_t n>
inline VectorInt<R,n> VectorInt<R,n>::operator+(const VectorInt<R,n> & other) const {
  VectorInt<R,n> sum;
  for (size_t i = 0; i < n; i++)
    sum[i] = this->_v[i] + other[i];
  return sum;
}

template<typename R, size_t n>
inline VectorInt<R,n> VectorInt<R,n>::operator-(const VectorInt<R,n> & other) const {
  VectorInt<R,n> diff;
  for (size_t i = 0; i < n; i++)
    diff[i] = this->_v[i] - other[i];
  return diff;
}

template<typename R, size_t n>
inline VectorInt<R,n> VectorInt<R,n>::operator*(const R & a) const {
  VectorInt<R,n> prod;
  for (size_t i = 0; i < n; i++)
    prod[i] = a * this->_v[i];
  return prod;
}

template<typename R, size_t n>
inline VectorInt<R,n> & VectorInt<R,n>::operator+=(const VectorInt<R,n> & other)
{
  for (size_t i = 0; i < n; i++)
    this->_v[i] += other[i];
  return (*this);
}

template<typename R, size_t n>
inline VectorInt<R,n> & VectorInt<R,n>::operator-=(const VectorInt<R,n> & other)
{
  for (size_t i = 0; i < n; i++)
    this->_v[i] -= other[i];
  return (*this);
}

template<typename R, size_t n>
inline VectorInt<R,n> VectorInt<R,n>::operator*(const SquareMatrixInt<R,n>& mat) const
{
  VectorInt<R,n> prod;
  for (size_t i = 0; i < n; i++) {
    prod[i] = 0;
    for (size_t j = 0; j < n; j++)
      prod[i] += mat(j, i) * this->_v[j];
  }
  return prod;
}

template<typename R, size_t n>
inline R VectorInt<R,n>::innerProduct(const VectorInt<R,n> & vec1,
				      const VectorInt<R,n> & vec2)
{
  R prod = 0;

  for (size_t i = 0; i < n; i++)
    prod += vec1[i]*vec2[i];
  
  return prod;
}

// printing
template<typename R, size_t n>
inline std::ostream& operator<<(std::ostream& os, const VectorInt<R,n>& vec)
{
  os << "Vector(";
  for (size_t i = 0; i < n-1; i++)
    os << vec[i] << ",";
  os << vec[n-1] <<  ")";
  return os;
}

template<typename R, size_t n>
inline std::ostream & VectorInt<R,n>::prettyPrint(std::ostream & os,
						  size_t upTo) const
{
  for (size_t i = 0; i < upTo; i++) {
    os << (*this)[i] << " ";
  }
  os << std::endl;
 
  return os;
}
