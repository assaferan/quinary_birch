// PolynomialFp

// create the zero polynomial
template<typename R, typename S>
PolynomialFp<R,S>::PolynomialFp(std::shared_ptr<const Fp<R,S>> GF)
{
  this->_GF = GF;
}

// create the constant polynomial 
template<typename R, typename S>
PolynomialFp<R,S>::PolynomialFp(const FpElement<R,S> & a)
{
  this->_GF = a.parent();
  std::multiset<size_t> empty_set;
  this->_mons[empty_set] = a;
}

// create the constant polynomial 
template<typename R, typename S>
PolynomialFp<R,S>::PolynomialFp(std::shared_ptr<const Fp<R,S>> GF, 
				  const R & a)
{
  this->_GF = GF;
  std::multiset<size_t> empty_set;
  FpElement<R,S> elt(GF, a);
  this->_mons[empty_set] = elt;
}

// create the polynomial x_i
template<typename R, typename S>
inline PolynomialFp<R,S>
PolynomialFp<R,S>::x(std::shared_ptr<const Fp<R,S>> GF, size_t i)
{
  PolynomialFp<R,S> x_i(GF);

  std::multiset<size_t> singleton;
  singleton.insert(i);
  FpElement<R,S> one(GF,1);
  x_i._mons[singleton] = one;

  return x_i;
}

template<typename R, typename S>
template<size_t n>
PolynomialFp<R,S>::PolynomialFp(const SquareMatrixFp<R,S,n> & q)
{
  this->_GF = q.baseRing();
  
  for (size_t i = 0; i < n; i++)
    for (size_t j = i; j < n; j++) {
      std::multiset<size_t> mon;
      mon.insert(i);
      mon.insert(j);
      this->_mons[mon] = q(i,j);
    }
  if (this->_GF->prime() != 2) {
    FpElement<R,S> two(this->_GF,2);
    for (size_t i = 0; i < n; i++) {
      std::multiset<size_t> mon;
      mon.insert(i);
      mon.insert(i);
      this->_mons[mon] /= two;
    }
  }
}

// returns the coefficient of a monomial
template<typename R, typename S>
inline FpElement<R,S>
PolynomialFp<R,S>::coefficient(const std::multiset<size_t> & mon) const
{
  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator it;
  
  it = this->_mons.find(mon);
  if (it == this->_mons.end()) {
    FpElement<R,S> zero(this->_GF,0);
    return zero;
  }
  return it->second;
}

template<typename R, typename S>
inline FpElement<R, S> PolynomialFp<R,S>::constCoefficient(void) const {
  std::multiset<size_t> empty_set;
  return this->coefficient(empty_set);
}

// coefficient of x_i
template<typename R, typename S>
inline FpElement<R,S> PolynomialFp<R,S>::coefficient(size_t i) const {
  std::multiset<size_t> singleton;
  singleton.insert(i);
  return this->coefficient(singleton);
}

// coefficient of x_i x_j
template<typename R, typename S>
inline FpElement<R,S> PolynomialFp<R,S>::coefficient(size_t i, size_t j) const
{
  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator it;
  std::multiset<size_t> mon;
  mon.insert(i);
  mon.insert(j);
  return this->coefficient(mon);
}

template<typename R, typename S>
inline PolynomialFp<R,S> PolynomialFp<R,S>::quadraticPart(void) const
{
  PolynomialFp<R,S> quad(this->_GF);

  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator i;
  
  for (i = this->_mons.begin(); i != this->_mons.end(); i++) {
    if ((i->first).size() == 2)
      quad._mons[i->first] = i->second;
  }
  
  return quad;
}

template<typename R, typename S>
inline std::vector< FpElement<R,S> > PolynomialFp<R,S>::linearPart(size_t rank) const
{
  std::vector< FpElement<R,S> > linear;
  for (size_t i = 0; i < rank; i++)
    linear.push_back(this->coefficient(i));
  return linear;
}

template<typename R, typename S>
inline int PolynomialFp<R,S>::degree(size_t i) const
{
  int deg = -1;
  int mon_deg;
  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator it;
  
  for (it = this->_mons.begin(); it != this->_mons.end(); it++) {
    if (it->second != 0) {
      mon_deg = it->first.count(i);
      if (deg < mon_deg)
	deg = mon_deg;
    }
  }

  return deg;
}

template<typename R, typename S>
inline PolynomialFp<R,S> & PolynomialFp<R,S>::operator=(const PolynomialFp<R,S> & other)
{
  if (this != (&other)) {
    this->_GF = other._GF;
    this->_mons = other._mons;
  }
  return (*this);
}

template<typename R, typename S>
inline PolynomialFp<R,S> & PolynomialFp<R,S>::operator=(const FpElement<R,S> & a)
{
  this->_GF = a.parent();
  this->_mons.clear();
  std::multiset<size_t> empty_set;
  this->mons[empty_set] = a;
  return (*this);
}

template<typename R, typename S>
inline PolynomialFp<R,S> PolynomialFp<R,S>::operator-() const
{
  PolynomialFp<R,S> neg(this->_GF);
  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator it;
  for (it = this->_mons.begin(); it != this->_mons.end(); it++) {
    neg.mons[it->first] = -it->second;
  }
  
  return neg;
}

template<typename R, typename S>
inline PolynomialFp<R,S> PolynomialFp<R,S>::operator+(const PolynomialFp<R,S> & other) const
{
  PolynomialFp<R,S> sum(this->_GF);
  FpElement<R,S> zero(this->_GF, 0);
  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator it, it2;
  for (it = this->_mons.begin(); it != this->_mons.end(); it++) {
    sum._mons[it->first] = it->second;
  }
  for (it = other._mons.begin(); it != other._mons.end(); it++) {
    it2 = sum._mons.find(it->first);
    if (it2 == sum._mons.end())
      sum._mons[it->first] = zero;
    sum._mons[it->first] += it->second;
  }
  
  return sum;
}

template<typename R, typename S>
inline PolynomialFp<R,S> PolynomialFp<R,S>::operator-(const PolynomialFp<R,S> & other) const
{
  PolynomialFp<R,S> diff(this->_GF);
  FpElement<R,S> zero(this->_GF, 0);
  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator it, it2;
  for (it = this->_mons.begin(); it != this->_mons.end(); it++) {
    diff._mons[it->first] = it->second;
  }
  for (it = other._mons.begin(); it != other._mons.end(); it++) {
    it2 = diff._mons.find(it->first);
    if (it2 == diff._mons.end())
      diff._mons[it->first] = zero;
    diff._mons[it->first] -= it->second;
  }
  
  return diff;
}

template<typename R, typename S>
inline PolynomialFp<R,S> PolynomialFp<R,S>::operator*(const FpElement<R,S> & a) const
{
  PolynomialFp<R,S> prod(this->_GF);

  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator it;
  for (it = this->_mons.begin(); it != this->_mons.end(); it++) {
    prod._mons[it->first] = a*it->second;
  }
  
  return prod;
}

template<typename R, typename S>
inline PolynomialFp<R,S> PolynomialFp<R,S>::operator*(const PolynomialFp<R,S> & other) const
{
  PolynomialFp<R,S> prod(this->_GF);
  FpElement<R,S> zero(this->_GF, 0);

  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator i,j,loc;
  for (i = this->_mons.begin(); i != this->_mons.end(); i++)
    for (j = other._mons.begin(); j != other._mons.end(); j++) {
      std::multiset<size_t> mon;
      mon.insert(i->first.begin(), i->first.end());
      mon.insert(j->first.begin(), j->first.end());
      loc = prod._mons.find(mon);
      if (loc == prod._mons.end())
	prod._mons[mon] = zero;
      prod._mons[mon] += (i->second)*(j->second);
  }
  
  return prod;
}

template<typename R, typename S>
inline PolynomialFp<R,S> & PolynomialFp<R,S>::operator+=(const PolynomialFp<R,S> & other)
{
  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator it, it2;
  
  for (it = other._mons.begin(); it != other._mons.end(); it++) {
    it2 = this->_mons.find(it->first);
    if (it2 == this->_mons.end())
      this->_mons[it->first] = it->second;
    else
      this->_mons[it->first] += it->second;
  }
  return (*this);
}

template<typename R, typename S>
inline PolynomialFp<R,S> & PolynomialFp<R,S>::operator-=(const PolynomialFp<R,S> & other)
{
  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator it, it2;
  
  for (it = other._mons.begin(); it != other._mons.end(); it++) {
    it2 = this->_mons.find(it->first);
    if (it2 == this->_mons.end())
      this->_mons[it->first] = -it->second;
    else
      this->_mons[it->first] -= it->second;
  }
  return (*this);
}

template<typename R, typename S>
inline PolynomialFp<R,S> & PolynomialFp<R,S>::operator*=(const PolynomialFp<R,S> & other)
{
  // Here we have no advantage doing it in place)
  (*this) = (*this)*other;
  return (*this);
}

template<typename R, typename S>
inline PolynomialFp<R,S> & PolynomialFp<R,S>::operator*=(const FpElement<R,S> & a)
{
  typename std::map<std::multiset<size_t>, FpElement<R,S> >::iterator it;
  for (it = this->_mons.begin(); it != this->_mons.end(); it++) {
    it->second *= a;
  }
  return (*this);
}

template<typename R, typename S>
inline FpElement<R,S>
PolynomialFp<R,S>::evaluate(const std::vector<FpElement<R,S> > & vec) const
{
  FpElement<R,S> res(this->_GF, 0);

  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator i;
  
  for (i = this->_mons.begin(); i != this->_mons.end(); i++) {
    FpElement<R,S> prod = i->second;
    std::multiset<size_t>::const_iterator j;
    for (j = i->first.begin(); j != i->first.end(); j++) {
      assert((*j) < vec.size());
      prod *= vec[*j];
    }
    res += prod;
  }

  return res;
}

template<typename R, typename S>
inline PolynomialFp<R,S>
PolynomialFp<R,S>::evaluate(const std::vector<PolynomialFp<R,S> > & vec) const
{
  PolynomialFp<R,S> res(this->_GF);

  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator i;
  
  for (i = this->_mons.begin(); i != this->_mons.end(); i++) {
    PolynomialFp<R,S> prod = i->second;
    std::multiset<size_t>::const_iterator j;
    for (j = i->first.begin(); j != i->first.end(); j++) {
      assert((*j) < vec.size());
      prod *= vec[*j];
    }
    res += prod;
  }

  return res;
}

// booleans

template<typename R, typename S>
inline bool PolynomialFp<R,S>::isZero() const
{
  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator i;
  for (i = this->_mons.begin(); i != this->_mons.end(); i++)
    if (i->second != 0)
      return false;

  return true;
}

template<typename R, typename S>
inline bool PolynomialFp<R,S>::operator==(const PolynomialFp<R,S> & other) const
{
  return ((*this)-other).isZero();
}

template<typename R, typename S>
inline bool PolynomialFp<R,S>::operator!=(const PolynomialFp<R,S> & other) const
{
  return !((*this)==other);
}

template<typename R, typename S>
inline bool PolynomialFp<R,S>::operator==(const FpElement<R,S> & a) const
{
  PolynomialFp<R,S> poly(a);
  return ((*this)==poly);
}

template<typename R, typename S>
inline bool PolynomialFp<R,S>::operator!=(const FpElement<R,S> & a) const
{
  return !((*this)==a);
}

template<typename R, typename S>
inline std::ostream& operator<<(std::ostream& os, const PolynomialFp<R,S>& poly)
{
  bool first = true;
  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator i;
  
  for (i = poly.monomials().begin(); i != poly.monomials().end(); i++) {
    if (i->second == 0)
      continue;
    if (!first)
      os << "+";
    if ((i->second != 1) || (i->first).empty())
      os << i->second;
    std::multiset<size_t>::const_iterator j;
    bool inner_first = true;
    for (j = i->first.begin(); j != i->first.end(); j++) {
      if (!inner_first)
	os << "*";
      os << "x_" << (*j);
      inner_first = false;
    }
    first = false;
  }

  if (first)
    os << "0";

  return os;
}
