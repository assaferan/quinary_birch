#include <cassert>

template<typename R>
NumberField<R>::NumberField(const UnivariatePolyInt<R> & mod)
  : _f(std::make_shared< const RationalField<R> >())
{
  std::shared_ptr<const RationalField<R> > QQ = this->_f.baseRing();
  for (int i = 0; i <= mod.degree(); i++) {
    Rational<R> coeff(mod.coefficient(i));
    this->_f += coeff*UnivariatePolyRat<R>::x(QQ, i);
  }
  _initAntic();
}

template<typename R>
void NumberField<R>::_initAntic(void)
{
  fmpq_poly_t poly;
  mpq_t* c_mpq;
  
  fmpq_poly_init(poly);
  assert(!_f.isZero());
      
  c_mpq = (mpq_t *)flint_malloc((_f.degree()+1)*sizeof(mpq_t));
    
  for (int i = 0; i <= _f.degree(); i++) {
    Rational<R> c = _f.coefficient(i);
    mpq_init(c_mpq[i]);
    mpq_set_si(c_mpq[i], birch_util::convertInteger<R,Z64>(c.num().num()), c.denom().num());
  }
    
  fmpq_poly_set_array_mpq(poly, (const mpq_t *)c_mpq, _f.degree()+1);
  nf_init(_nf_antic, poly);

  fmpq_poly_clear(poly);
  for (int i = 0; i <= _f.degree(); i++)
    mpq_clear(c_mpq[i]);
  flint_free(c_mpq);
}





