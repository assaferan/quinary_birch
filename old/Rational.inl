// Rationals - the field of rational numbers

template<typename R>
Rationals<R>& Rationals<R>::getInstance()
{
  static Rationals<R> instance;

  return instance;
}

template<typename R>
Rational<R> Rationals<R>::zero() const
{
  R z = 0;
  return z;
}

template<typename R>
Rational<R> Rationals<R>::one() const
{
  R z = 1;
  return z;
}

// implementation of methods for Rational<R>

// c-tors
template<typename R>
Rational<R>::Rational(const R& num, const R& denom)
  : FieldElement(Rationals<R>::getInstance()), num_(num), denom_(denom)
{ this->reduce(); }

template<typename R>
Rational<R>::Rational(const R & num)
  : FieldElement(Rationals<R>::getInstance()), num_(num), denom_(1)
{}
  
// default c-tor
template<typename R>
Rational<R>::Rational() : FieldElement(Rationals<R>::getInstance()), num_(0), denom_(1)
{}

// copy c-tor
template<typename R>
Rational<R>::Rational(const Rational<R> & other)
  : FieldElement(other.parent()), num_(other.num_), denom_(other.denom_)
{}

// assignment
template<typename R>
Rational<R> & Rational<R>::operator=(const Rational<R> & b)
{
  if (this != &b) {
    num_ = b.num_;
    denom_ = b.denom_;
  }
  return (*this);
}

// arithmetic

template <typename R>
Rational<R> Rational<R>::operator+(const Rational<R> &other) const
{
  Rational<R> sum;
  sum.num_ = (this->num_) * other.denom_ + (this->denom_) * other.num_;
  sum.denom_ = (this->denom_) * other.denom_;
  sum.reduce();
  return sum;
}

template <typename R>
Rational<R> Rational<R>::operator*(const Rational<R> &other) const
{
  Rational<R> prod;
  prod.num_ = (this->num_) * other.num_;
  prod.denom_ = (this->denom_) * other.denom_;
  prod.reduce();
  return prod;
}

template <typename R>
Rational<R> Rational<R>::operator*(const R & b) const {
  Rational<R> b_rat(b);
  return (*this)*b_rat;
}

template <typename R>
Rational<R> Rational<R>::operator/(const Rational<R> &other) const
  {
  Rational<R> prod;
  prod.num_ = (this->num_) * other.denom_;
  prod.denom_ = (this->denom_) * other.num_;
  prod.reduce();
  return prod;
}

template <typename R>
Rational<R> Rational<R>::operator/(const R & b) const {
  Rational<R> b_rat(b);
  return (*this)/b_rat;
}

template <typename R>
bool Rational<R>::operator==(const Rational<R> &other) const
{
  return (this->num_) * other.denom_ == other.num_ * (this->denom_); 
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
  R d = Math<R>::gcd(num_, denom_);
  num_ /= d;
  denom_ /= d;
  return;
}

// other
template<typename R>
R Rational<R>::floor() const    
{
  R num = num_;
  R denom = denom_;
    
  if (denom < 0) {
    denom = -denom;
    num = - num;
  }
      
  return ((num >= parent_.zero()) ? num : (num - denom + 1)) / denom;
}

template<typename R>
R Rational<R>::ceiling() const
{
  R num = num_;
  R denom = denom_;
    
  if (denom < 0) {
    denom = -denom;
    num = - num;
  }
      
  return ((num >= parent_.zero()) ? (num + denom - 1) : num) / denom;

}

// other
template <typename R>
Rational<R> operator*(R b, const Rational<R> & r) {
  return r*b;
}

template <typename R>
Rational<R> operator-(R b, const Rational<R> & r) {
  Rational<R> b_rat(b);
  return b_rat-r;
}

template <typename R>
Rational<R> operator+(R b, const Rational<R> & r) {
  Rational<R> b_rat(b);
  return b_rat+r;
}

template <typename R>
Rational<R> operator/(R b, const Rational<R> & r) {
  Rational<R> b_rat(b);
  return b_rat/r;
}

template <typename R>
std::ostream& operator<<(std::ostream & os, const Rational<R> & r)
{
  R one = 1;
  if (r.denom() == one) return os << r.num();
  if (r.denom() == -one) return os << -r.num();
  os << r.num() << "/" << r.denom();
  return os;
}
