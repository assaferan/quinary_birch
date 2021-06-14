template<typename R>
NumberField<R>::NumberField(const UnivariatePolyInt<R> & mod)
{
  for (int i = 0; i <= mod.degree(); i++) {
    Rational<R> coeff(mod.coefficient(i));
    this->_f += coeff*UnivariatePolyRat<R>::x(i);
  }
}







