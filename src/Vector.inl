// Vector

template<class R, class Parent, size_t n>
bool Vector<R,Parent,n>::operator==(const Vector<R,Parent,n> & other) const
{
  for (size_t i = 0; i < n; i++)
    if (this->v[i] != other[i]) return false;
  return true;
}

// for comparison in standard containers such as std::set and std::map
template<class R, class Parent, size_t n>
bool Vector<R,Parent,n>::operator<(const Vector<R,Parent,n> & other) const
{
  for (size_t i = 0; i < n; i++) {
    if (this->v[i] > other[i]) return false;
    if (this->v[i] < other[i]) return true;
  }
  return false;
}

template<class R, class Parent, size_t n>
Vector<R,Parent,n> Vector<R,Parent,n>::operator-() const {
  Vector<R,Parent,n> neg(_base);
  for (size_t i = 0; i < n; i++)
    neg[i] = -this->v[i];
  return neg;
}

template<class R, class Parent, size_t n>
Vector<R,Parent,n> Vector<R,Parent,n>::operator+(const Vector<R,Parent,n> & other) const {
  Vector<R,Parent,n> sum(_base);
  for (size_t i = 0; i < n; i++)
    sum[i] = this->v[i] + other[i];
  return sum;
}

template<class R, class Parent, size_t n>
Vector<R,Parent,n> Vector<R,Parent,n>::operator-(const Vector<R,Parent,n> & other) const {
  Vector<R,Parent,n> diff(_base);
  for (size_t i = 0; i < n; i++)
    diff[i] = this->v[i] - other[i];
  return diff;
}

template<class R, class Parent, size_t n>
Vector<R,Parent,n> Vector<R,Parent,n>::operator*(const R & a) const {
  Vector<R,Parent,n> prod(_base);
  for (size_t i = 0; i < n; i++)
    prod[i] = a * this->v[i];
  return prod;
}

template<class R, class Parent, size_t n>
Vector<R,Parent,n> & Vector<R,Parent,n>::operator+=(const Vector<R,Parent,n> & other)
{
  for (size_t i = 0; i < n; i++)
    this->v[i] += other[i];
  return (*this);
}

template<class R, class Parent, size_t n>
Vector<R,Parent,n> & Vector<R,Parent,n>::operator-=(const Vector<R,Parent,n> & other)
{
  for (size_t i = 0; i < n; i++)
    this->v[i] -= other[i];
  return (*this);
}

template<class R, class Parent, size_t n>
Vector<R,Parent,n> Vector<R,Parent,n>::operator*(const SquareMatrix<R, n>& mat) const
{
  Vector<R,Parent,n> prod;
  for (size_t i = 0; i < n; i++) {
    prod[i] = _base->zero();
    for (size_t j = 0; j < n; j++)
      prod[i] += mat(j, i) * this->v[j];
  }
  return prod;
}

template<class R, class Parent, size_t n>
R Vector<R,Parent,n>::innerProduct(const Vector<R,Parent,n> & vec1,
				    const Vector<R,Parent,n> & vec2)
{
  R prod = _base->zero();

  for (size_t i = 0; i < n; i++)
    prod += vec1[i]*vec2[i];
  
  return prod;
}

// printing
template<class R, class Parent, size_t n>
std::ostream& operator<<(std::ostream& os, const Vector<R,Parent,n>& vec)
{
  os << "Vector(";
  for (size_t i = 0; i < n-1; i++)
    os << vec[i] << ",";
  os << vec[n-1] <<  ")";
  return os;
}

template<class R, class Parent, size_t n>
std::ostream & Vector<R,Parent,n>::prettyPrint(std::ostream & os,
					 size_t upTo) const
{
  for (size_t i = 0; i < upTo; i++) {
    os << (*this)[i] << " ";
  }
  os << std::endl;
 
  return os;
}
