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
    for (int j = 0; j < d; j++) {
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

template<typename R>
void NumberFieldElement<R>::_initAntic(void) {
  nf_elem_init(_nf_elt_antic, _K->antic());
  
  if (!(_elt.isZero())) {
    fmpq_poly_t poly;
    mpq_t* c_mpq;

    fmpq_poly_init(poly);
    c_mpq = (mpq_t *)flint_malloc((_elt.degree()+1)*sizeof(mpq_t));
    for (int i = 0; i <= _elt.degree(); i++) {
      Rational<R> c = _elt.coefficient(i);
      mpq_init(c_mpq[i]);
      mpq_set_si(c_mpq[i], birch_util::convertInteger<R,Z64>(c.num().num()), birch_util::convertInteger<R,W64>(c.denom().num()));
    }
    
    fmpq_poly_set_array_mpq(poly, (const mpq_t *)c_mpq, _elt.degree()+1);
    nf_elem_set_fmpq_poly(_nf_elt_antic, poly, _K->antic());

    fmpq_poly_clear(poly);
    
    for (int i = 0; i <= _elt.degree(); i++)
      mpq_clear(c_mpq[i]);
    flint_free(c_mpq);
  }
  else {
    nf_elem_zero(_nf_elt_antic, _K->antic());
  }

  return;
}
