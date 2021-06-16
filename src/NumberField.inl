template<typename R>
NumberField<R>::NumberField(const UnivariatePolyInt<R> & mod)
{
  std::shared_ptr<const RationalField<R> > QQ = std::make_shared< const RationalField<R> >();
  for (int i = 0; i <= mod.degree(); i++) {
    Rational<R> coeff(mod.coefficient(i));
    this->_f += coeff*UnivariatePolyRat<R>::x(QQ, i);
  }
}







