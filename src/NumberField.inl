template<typename R>
NumberField<R>::NumberField(const UnivariatePolyInt<R> & mod)
  : _f(std::make_shared< const RationalField<R> >())
{
  std::shared_ptr<const RationalField<R> > QQ = this->_f.baseRing();
  for (int i = 0; i <= mod.degree(); i++) {
    Rational<R> coeff(mod.coefficient(i));
    this->_f += coeff*UnivariatePolyRat<R>::x(QQ, i);
  }
}







