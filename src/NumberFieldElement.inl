#include <unordered_map>

template<typename R>
inline NumberFieldElement<R> & NumberFieldElement<R>::operator=(const R & a)
{
  this->_elt = a;
  return (*this);
}

template<typename R>
inline NumberFieldElement<R> NumberFieldElement<R>::operator-(void) const
{
  NumberFieldElement<R> neg(this->_K);
  neg._elt = -(this->_elt);
  
  return neg;
}

template<typename R>
inline NumberFieldElement<R>
NumberFieldElement<R>::operator+(const NumberFieldElement<R> & other) const
{
  NumberFieldElement<R> sum(this->_K);
  sum._elt = (this->_elt) + other._elt;
  
  return sum;
}

template<typename R>
inline NumberFieldElement<R>
NumberFieldElement<R>::operator-(const NumberFieldElement<R> & other) const
{
  return (*this)+(-other);
}

template<typename R>
inline NumberFieldElement<R>
NumberFieldElement<R>::operator*(const NumberFieldElement<R> & other) const
{
  NumberFieldElement<R> prod(this->_K);
  prod._elt = (this->_elt)  * other._elt % (this->_K->modulus());

  return prod;
}

template<typename R>
inline NumberFieldElement<R>
NumberFieldElement<R>::operator/(const NumberFieldElement<R> & other) const
{
  return (*this) * other.inverse();
}

template<typename R>
inline NumberFieldElement<R> NumberFieldElement<R>::operator*(const R & a) const
{
  NumberFieldElement<R> prod(this->_K);
  prod._elt = (this->_elt) * a;

  return prod;
}

template<typename R>
inline NumberFieldElement<R> NumberFieldElement<R>::operator/(const R & a) const
{
  NumberFieldElement<R> quo(this->_K);
  quo._elt = (this->_elt) / a;

  return quo;
}

template<typename R>
inline NumberFieldElement<R> &
NumberFieldElement<R>::operator+=(const NumberFieldElement<R> & other)
{
  this->_elt += other._elt;
  return (*this);
}

template<typename R>
inline NumberFieldElement<R> &
NumberFieldElement<R>::operator-=(const NumberFieldElement<R> & other)
{
  this->_elt -= other._elt;
  return (*this);
}

template<typename R>
inline NumberFieldElement<R> &
NumberFieldElement<R>::operator*=(const NumberFieldElement<R> & other)
{
  *this = (*this)*other;
  return (*this);
}

template<typename R>
inline NumberFieldElement<R> &
NumberFieldElement<R>::operator/=(const NumberFieldElement<R> & other)
{
  *this = (*this)/other;
  return (*this);
}

template<typename R>
inline NumberFieldElement<R>& NumberFieldElement<R>::operator*=(const R & a)
{
  this->_elt *= a;
  return (*this);
}

template<typename R>
inline NumberFieldElement<R>& NumberFieldElement<R>::operator/=(const R & a)
{
  this->_elt /= a;
  return (*this);
}

template<typename R>
inline NumberFieldElement<R>
NumberFieldElement<R>::inverse(void) const
{
  std::shared_ptr<const RationalField<R> > QQ = this->_elt.baseRing();
  UnivariatePolyRat<R> s(QQ);
  UnivariatePolyRat<R> t(QQ);
  
  UnivariatePolyRat<R>::xgcd(this->_elt, this->_K->modulus(), s, t);
  
  NumberFieldElement<R> inv(this->_K, s);

  assert(((*this)*inv).isOne());
  return inv;
}

template<typename R>
inline NumberFieldElement<R>
NumberFieldElement<R>::zero(std::shared_ptr<const NumberField<R> > fld)
{
  NumberFieldElement z(fld);
  z.makeZero();

  return z;
}

template<typename R>
inline NumberFieldElement<R>
NumberFieldElement<R>::one(std::shared_ptr<const NumberField<R> > fld)
{
  NumberFieldElement z(fld);
  z.makeOne();

  return z;
}

template<typename R>
inline bool NumberFieldElement<R>::isZero(void) const
{
  return _elt.isZero();
}

template<typename R>
inline bool NumberFieldElement<R>::isOne(void) const
{
  return _elt.isOne();
}

template<typename R>
inline NumberFieldElement<R> & NumberFieldElement<R>::operator=(const NumberFieldElement<R> & other)
{
  if (this != &other) {
    this->_K = other._K;
    this->_elt = other._elt;
  }
  return *this;
}

template<typename R>
inline MatrixRat<R> NumberFieldElement<R>::_multByMatrix(void) const
{
  int d = _K->modulus().degree();
  // creating the basis {1,x,....,x^{d-1}}
  std::vector< NumberFieldElement<R> > basis;
  for (int i = 0; i < d; i++) {
    NumberFieldElement<R> a(_K, UnivariatePolyRat<R>::x(_elt.baseRing(), i));
    basis.push_back(a);
  }

  MatrixRat<R> mat(_elt.baseRing(), d, d);

  for (int i = 0; i < d; i++) {
    UnivariatePolyRat<R> mul = ((*this)*basis[i])._elt;
    for (int j = 0; j < d; i++) {
      mat(i,j) = mul.coefficient(j);
    }
  }

  return mat;

}

template<typename R>
inline UnivariatePolyInt<R> NumberFieldElement<R>::minimalPolynomial(void) const
{
  MatrixRat<R> mult_mat = this->_multByMatrix();
  UnivariatePolyRat<R> char_poly = mult_mat.charPoly();
  
  Integer<R> denom = Integer<R>::one();
  std::vector< Integer<R> > coeffs_int;
  for (int i = 0; i <= char_poly.degree(); i++)
    denom = denom.lcm(char_poly.coefficient(i).denom());
  char_poly *= denom;
  for (int i = 0; i <= char_poly.degree(); i++) 
    coeffs_int.push_back(char_poly.coefficient(i).floor());
  
  UnivariatePolyInt<R> char_poly_int(coeffs_int);
  
  std::unordered_map< UnivariatePolyInt<R>, size_t > fac = char_poly_int.factor();

  UnivariatePolyInt<R> min_poly = Integer<R>::one();

  // Here we use the fact that all algebraic numbers are separable
  for (std::pair< UnivariatePolyInt<R>, size_t > fa : fac) {
    min_poly *= fa.first;
  }

  return min_poly;
}
