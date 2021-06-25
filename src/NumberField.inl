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
  fmpq_poly_init(poly);
  assert(!_f.isZero());
      
  mpq_t* c_mpq = new mpq_t[_f.degree()+1];
    
  for (int i = 0; i <= _f.degree(); i++) {
    Rational<R> c = _f.coefficient(i);
    fmpq_set_si(c_mpq[i], c.num().num(), c.denom().num());
  }
    
  fmpq_poly_set_array_mpq(poly, c_mpq, _f.degree()+1);
  nf_init(_nf_antic, poly);

  fmpq_poly_clear(poly);
  delete[] c_mpq;
}





