// Integers - the ring of integers

template<typename R>
Integers<R>& Integers<R>::getInstance()
{
  static Integers<R> instance;

  return instance;
}

template<typename R>
Integer<R> Integers<R>::zero() const
{
  R z = 0;
  return z;
}

template<typename R>
Integer<R> Integers<R>::one() const
{
  R z = 1;
  return z;
}

// implementation of methods for Integer<R>

// c-tors
template<typename R>
Integer<R>::Integer(const R& num, const R& denom)
  : RingElement(Integers<R>::getInstance()), num_(num)
{ this->reduce(); }

template<typename R>
Integer<R>::Integer(const R & num)
  : RingElement(Integers<R>::getInstance()), num_(num)
{}
  
// default c-tor
template<typename R>
Integer<R>::Integer() : RingElement(Integers<R>::getInstance()), num_(0)
{}

// copy c-tor
template<typename R>
Integer<R>::Integer(const Integer<R> & other)
  : RingElement(other.parent()), num_(other.num_), denom_(other.denom_)
{}

// assignment
template<typename R>
Integer<R> & Integer<R>::operator=(const Integer<R> & b)
{
  if (this != &b) {
    num_ = b.num_;
  }
  return (*this);
}

// arithmetic

template <typename R>
Integer<R> Integer<R>::operator+(const Integer<R> &other) const
{
  Integer<R> sum;
  sum.num_ = this->num_ + other.num_;
  
  return sum;
}

template <typename R>
Integer<R> Integer<R>::operator*(const Integer<R> &other) const
{
  Integer<R> prod;
  prod.num_ = (this->num_) * other.num_;
 
  return prod;
}

template <typename R>
Integer<R> Integer<R>::operator*(const R & b) const {
  Integer<R> b_int(b);
  return (*this)*b_int;
}

template <typename R>
Integer<R> Integer<R>::operator/(const Integer<R> &other) const
{
  Integer<R> quo;
  quo.num_ = (this->num_) / other.denom_;
 
  return prod;
}

template <typename R>
Integer<R> Integer<R>::operator/(const R & b) const {
  Integer<R> b_int(b);
  return (*this)/b_int;
}

template <typename R>
bool Integer<R>::operator==(const Integer<R> &other) const
{
  return (this->num_ == other.num_); 
}

template <typename R>
bool Integer<R>::operator<(const Integer<R> &other) const
{
  return (this->num_ < other.num_); 
}

// other
template <typename R>
Integer<R> operator*(R b, const Integer<R> & r) {
  return r*b;
}

template <typename R>
Integer<R> operator-(R b, const Integer<R> & r) {
  Integer<R> b_int(b);
  return b_int-r;
}

template <typename R>
Integer<R> operator+(R b, const Integer<R> & r) {
  Integer<R> b_int(b);
  return b_int+r;
}

template <typename R>
Integer<R> operator/(R b, const Integer<R> & r) {
  Integer<R> b_int(b);
  return b_int/r;
}

template <typename R>
std::ostream& operator<<(std::ostream & os, const Integer<R> & r)
{
  os << r.num();
  return os;
}
