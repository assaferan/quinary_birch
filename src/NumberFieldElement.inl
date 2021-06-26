#include <unordered_map>

#include "antic/nf_elem.h"

#include "birch_util.h"

template<typename R>
NumberFieldElement<R>::NumberFieldElement(const NumberFieldElement<R> & other)
  : _K(other._K)
{
  nf_elem_init(_nf_elt_antic, _K->antic());
  nf_elem_set(this->_nf_elt_antic, other._nf_elt_antic, this->_K->antic());
}

template<typename R>
inline NumberFieldElement<R> & NumberFieldElement<R>::operator=(const R & a)
{
  nf_elem_set_si(this->_nf_elt_antic, birch_util::convertInteger<R,slong>(a),this->_K->antic());
  
  return (*this);
}

template<typename R>
inline NumberFieldElement<R> NumberFieldElement<R>::operator-(void) const
{
  NumberFieldElement<R> neg(this->_K);

  nf_elem_neg(neg._nf_elt_antic, this->_nf_elt_antic, this->_K->antic());
  
  return neg;
}

template<typename R>
inline NumberFieldElement<R>
NumberFieldElement<R>::operator+(const NumberFieldElement<R> & other) const
{
  NumberFieldElement<R> sum(this->_K);
  
  nf_elem_add(sum._nf_elt_antic, this->_nf_elt_antic, other._nf_elt_antic, this->_K->antic());
  
  return sum;
}

template<typename R>
inline NumberFieldElement<R>
NumberFieldElement<R>::operator-(const NumberFieldElement<R> & other) const
{
  // return (*this)+(-other);
  NumberFieldElement<R> diff(this->_K);

  nf_elem_sub(diff._nf_elt_antic, this->_nf_elt_antic, other._nf_elt_antic, this->_K->antic());
  
  return diff;
}

template<typename R>
inline NumberFieldElement<R>
NumberFieldElement<R>::operator*(const NumberFieldElement<R> & other) const
{
  NumberFieldElement<R> prod(this->_K);

  nf_elem_mul(prod._nf_elt_antic, this->_nf_elt_antic, other._nf_elt_antic, this->_K->antic());
  
  return prod;
}

template<typename R>
inline NumberFieldElement<R>
NumberFieldElement<R>::operator/(const NumberFieldElement<R> & other) const
{
  NumberFieldElement<R> quo(this->_K);
  
  nf_elem_div(quo._nf_elt_antic, this->_nf_elt_antic, other._nf_elt_antic, this->_K->antic());
  
  return quo;
}

template<typename R>
inline NumberFieldElement<R> NumberFieldElement<R>::operator*(const R & a) const
{
  NumberFieldElement<R> prod(this->_K);

  nf_elem_scalar_mul_si(prod._nf_elt_antic, this->_nf_elt_antic, birch_util::convertInteger<R,slong>(a), this->_K->antic());

  return prod;
}

template<typename R>
inline NumberFieldElement<R> NumberFieldElement<R>::operator/(const R & a) const
{
  NumberFieldElement<R> quo(this->_K);

  nf_elem_scalar_div_si(quo._nf_elt_antic, this->_nf_elt_antic, birch_util::convertInteger<R,slong>(a), this->_K->antic());
 
  return quo;
}

template<typename R>
inline NumberFieldElement<R> &
NumberFieldElement<R>::operator+=(const NumberFieldElement<R> & other)
{
  nf_elem_add(this->_nf_elt_antic, this->_nf_elt_antic, other._nf_elt_antic, this->_K->antic());
  
  return (*this);
}

template<typename R>
inline NumberFieldElement<R> &
NumberFieldElement<R>::operator-=(const NumberFieldElement<R> & other)
{
  nf_elem_sub(this->_nf_elt_antic, this->_nf_elt_antic, other._nf_elt_antic, this->_K->antic());
  
  return (*this);
}

template<typename R>
inline NumberFieldElement<R> &
NumberFieldElement<R>::operator*=(const NumberFieldElement<R> & other)
{
  nf_elem_mul(this->_nf_elt_antic, this->_nf_elt_antic, other._nf_elt_antic, this->_K->antic());
  
  return (*this);
}

template<typename R>
inline NumberFieldElement<R> &
NumberFieldElement<R>::operator/=(const NumberFieldElement<R> & other)
{
  nf_elem_div(this->_nf_elt_antic, this->_nf_elt_antic, other._nf_elt_antic, this->_K->antic());
  
  return (*this);
}

template<typename R>
inline NumberFieldElement<R>& NumberFieldElement<R>::operator*=(const R & a)
{
  nf_elem_scalar_mul_si(this->_nf_elt_antic, this->_nf_elt_antic, birch_util::convertInteger<R,slong>(a), this->_K->antic());
   
  return (*this);
}

template<typename R>
inline NumberFieldElement<R>& NumberFieldElement<R>::operator/=(const R & a)
{
  nf_elem_scalar_div_si(this->_nf_elt_antic, this->_nf_elt_antic, birch_util::convertInteger<R,slong>(a), this->_K->antic());
  
  return (*this);
}

template<typename R>
inline NumberFieldElement<R>
NumberFieldElement<R>::inverse(void) const
{
  NumberFieldElement<R> inv(this->_K);

  nf_elem_inv(inv._nf_elt_antic, this->_nf_elt_antic, this->_K->antic());

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
  return nf_elem_is_zero(_nf_elt_antic, _K->antic());
}

template<typename R>
inline bool NumberFieldElement<R>::isOne(void) const
{
  return nf_elem_is_one(_nf_elt_antic, _K->antic());
}

template<typename R>
inline NumberFieldElement<R> & NumberFieldElement<R>::operator=(const NumberFieldElement<R> & other)
{
  if (this != &other) {
    this->_K = other._K;
    nf_elem_set(this->_nf_elt_antic, other._nf_elt_antic, this->_K->antic());
  }
  return (*this);
}

// !! TODO - replace this by the appropriate method from antic nf_elem
template<typename R>
inline MatrixRat<R> NumberFieldElement<R>::_multByMatrix(void) const
{
  int d = _K->modulus().degree();
  // creating the basis {1,x,....,x^{d-1}}
  std::vector< NumberFieldElement<R> > basis;
  for (int i = 0; i < d; i++) {
    NumberFieldElement<R> a(_K, UnivariatePolyRat<R>::x(_K->modulus().baseRing(), i));
    basis.push_back(a);
  }

  MatrixRat<R> mat(_K->modulus().baseRing(), d, d);

  fmpq_t fcoeff;
  Z num, denom;
  fmpq_init(fcoeff);
  for (int i = 0; i < d; i++) {
    nf_elem_t mul = ((*this)*basis[i])._nf_elt_antic;
    for (int j = 0; j < d; j++) {
      nf_elem_get_coeff_fmpq(fcoeff, mul, j, _K->antic());
      fmpq_get_mpz_frac(num.get_mpz_t(), denom.get_mpz_t(), fcoeff);
      Rational<R> c(birch_util::convertInteger<Z,R>(num), birch_util::convertInteger<Z,R>(denom));
      mat(i,j) = c;
    }
  }
  fmpq_clear(fcoeff);

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
}

template<typename R>
void NumberFieldElement<R>::_initAntic(const Rational<R> & a) {
  UnivariatePolyRat<R> f(a);
  _initAntic(f);
  return;
}

template<typename R>
void NumberFieldElement<R>::_initAntic(const UnivariatePolyInt<R> & f) {
  UnivariatePolyRat<R> f_int(f);
  _initAntic(f_int);
  return;
}

template<typename R>
void NumberFieldElement<R>::_initAntic(const UnivariatePolyRat<R> & f) {
  nf_elem_init(_nf_elt_antic, _K->antic());

  if (!(f.isZero())) {
    fmpq_poly_t poly;
    mpq_t* c_mpq;

    fmpq_poly_init(poly);
    c_mpq = (mpq_t *)flint_malloc((f.degree()+1)*sizeof(mpq_t));
    for (int i = 0; i <= f.degree(); i++) {
      Rational<R> c = f.coefficient(i);
      mpq_init(c_mpq[i]);
      mpq_set_si(c_mpq[i], birch_util::convertInteger<R,slong>(c.num().num()), birch_util::convertInteger<R,ulong>(c.denom().num()));
    }
    
    fmpq_poly_set_array_mpq(poly, (const mpq_t *)c_mpq, f.degree()+1);
    nf_elem_set_fmpq_poly(_nf_elt_antic, poly, _K->antic());

    fmpq_poly_clear(poly);
    
    for (int i = 0; i <= f.degree(); i++)
      mpq_clear(c_mpq[i]);
    flint_free(c_mpq);
  }
  else {
    nf_elem_zero(_nf_elt_antic, _K->antic());
  }

  return;
}
