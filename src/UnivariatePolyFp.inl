
template<typename R, typename S>
inline UnivariatePolyInt<R> UnivariatePolyFp<R,S>::lift(void) const
{
  UnivariatePolyInt<R> ret;
  for (size_t i = 0; i < this->_coeffs.size(); i++)
    ret._coeffs[i] = this->coefficient(i).lift();

  return ret;
}

template<typename R, typename S>
inline UnivariatePolyFp<R,S>
UnivariatePolyFp<R,S>::powMod(size_t m, const UnivariatePolyFp<R,S> & f) const
{
  FpElement<R,S> one(this->baseRing(), 1); 
  UnivariatePolyFp<R,S> res(one);
  UnivariatePolyFp<R,S> q(this->baseRing());
  UnivariatePolyFp<R,S> r(this->baseRing());

  for (size_t i = 0; i < m; i++) {
    res = (*this)*res;
    UnivariatePolyFp<R,S>::divRem(res, f, q, r);
    res = r;
  }
  
  return res;
}

template<typename R, typename S>
inline std::vector< UnivariatePolyFp<R,S> >
UnivariatePolyFp<R,S>::_czEqDegPartialFactor(size_t r) const
{
  
  if (this->degree() == static_cast<int>(r)) {
    std::vector< UnivariatePolyFp<R,S> > ret(1, *this);
    return ret;
  }
  
  while (true) {
    UnivariatePolyFp<R,S> b(this->baseRing());
    
    for (int i = 0; i < this->degree(); i++)
      b._coeffs.push_back(this->baseRing()->random());

    b._eliminateDeg();
    
    Integer<R> p_r = this->baseRing()->prime()^r;
    size_t m = (birch_util::convert_Integer<R,size_t>(p_r.num()) - 1) / 2;

    UnivariatePolyFp<R,S> b_m = b.powMod(m, *this);
    UnivariatePolyFp<R,S> factor(this->baseRing());
    UnivariatePolyFp<R,S> b_m_shifted(this->baseRing());
    FpElement<R,S> shift(this->baseRing(), this->baseRing()->prime()-1);
    for (size_t i = 0; i < 3; i++) {
      b_m_shifted = b_m + shift;
      shift++;
      factor = gcd(b_m_shifted, *this);
      if ((factor.degree() != 0) && (factor.degree() != this->degree()))
	return factor._czEqDegFactor(r);
    }
  }
}

template<typename R, typename S>
inline std::vector< UnivariatePolyFp<R,S> >
UnivariatePolyFp<R,S>::_czEqDegFactor(size_t r) const
{
  std::vector< UnivariatePolyFp<R,S> > facs, partial;
  UnivariatePolyFp<R,S> f = (*this);

  for (size_t i = 0; (r*i < this->_coeffs.size()) && (f.degree() > 0); i++) {
    if (f.degree() == static_cast<int>(r)) {
      facs.push_back(f);
      return facs;
    }
    partial = f._czEqDegPartialFactor(r);
    for ( UnivariatePolyFp<R,S> g : partial) {
      f /= g;
      facs.push_back(g);
    }
  }
  return facs;
}

template<typename R, typename S>
inline std::vector< UnivariatePolyFp<R,S> >
UnivariatePolyFp<R,S>::_czDistinctDegFactor(void) const
{
  FpElement<R,S> one = this->_base->one();
  size_t n = this->degree();
  std::vector< UnivariatePolyFp<R,S> > facs(n, one);
  
  
  size_t beta = n / 2;
  size_t l = floor(sqrt(beta));
  if (l == 0) l = 1;
  size_t m = (beta + l - 1) / l;
  Integer<R> p = this->_base->prime();
  Integer<R> p_l = p^l;

  std::vector< UnivariatePolyFp<R,S> > h, H, I;

  UnivariatePolyFp<R,S> x_p_i = UnivariatePolyFp<R,S>::x(this->_base);
  for (size_t i = 0; i < l; i++) {
    h.push_back(x_p_i);
    x_p_i = x_p_i.powMod(birch_util::convert_Integer<R,size_t>(p.num()), *this);
  }

  x_p_i = x(this->baseRing());
  for (size_t i = 0; i <= m; i++) {
    x_p_i = x_p_i.powMod(birch_util::convert_Integer<R, size_t>(p_l.num()), *this);
    H.push_back(x_p_i);
  }

  UnivariatePolyFp<R,S> prod(this->baseRing());
  UnivariatePolyFp<R,S> q(this->baseRing());
  UnivariatePolyFp<R,S> r(this->baseRing());
  UnivariatePolyFp<R,S> g(this->baseRing());
  UnivariatePolyFp<R,S> mul(this->baseRing());
  UnivariatePolyFp<R,S> diff(this->baseRing());
  
  for (size_t i = 0; i <= m; i++) {
    prod = one;
    for (size_t j = 0; j < l; j++) {
      diff = H[i]-h[j];
      mul = prod*diff;
      divRem(mul, *this, q, r);
      prod = r;
    }
    I.push_back(prod);
  }

  UnivariatePolyFp<R, S> f = *this;
  
  for (size_t i = 0; i <= m; i++) {
    g = gcd(f, I[i]);
    f /= g;
    for (size_t j = l; j > 0; j--) {
      diff = H[i] - h[j-1];
      assert( (j <= l*(i+1)) && (l*(i+1)-j < facs.size()));
      facs[l*(i+1)-j] = UnivariatePolyFp<R,S>::gcd(g, diff);
      g /= facs[l*(i+1)-j];
    }
  }

  if (f.degree() >= 1)
    facs[f.degree()-1] = f;

  return facs;
}

template<typename R, typename S>
inline std::vector< UnivariatePolyFp<R,S> >
UnivariatePolyFp<R,S>::_sqfFactor(void) const
{
  std::vector< UnivariatePolyFp<R,S> > fac;

  std::vector< UnivariatePolyFp<R,S> > dist = this->_czDistinctDegFactor();

  for (size_t r = 0; r < dist.size(); r++) {
    std::vector< UnivariatePolyFp<R,S> > eq_deg =
      dist[r]._czEqDegFactor(r+1);

    fac.insert(fac.end(), eq_deg.begin(), eq_deg.end());
    
  }
  
  return fac;
}

template<typename R, typename S>
inline UnivariatePolyFp<R,S>
UnivariatePolyFp<R,S>::gcd(const UnivariatePolyFp<R,S> & f,
			   const UnivariatePolyFp<R,S> & g)
{
  UnivariatePolyFp<R,S> q(f.baseRing());
  UnivariatePolyFp<R,S> r(f.baseRing());
  UnivariatePolyFp<R,S> r_minus(f.baseRing());
  UnivariatePolyFp<R,S> r_plus(f.baseRing());
  
  r_minus = f;
  r = g;
  
  while (!r.isZero()) {
    divRem(r_minus, r, q, r_plus);
    r_minus = r;
    r = r_plus;
  }

  // alwasy return a monic factor
  return r_minus / r_minus.lead();
}

template<typename R, typename S>
inline UnivariatePolyFp<R,S>
UnivariatePolyFp<R,S>::xgcd(const UnivariatePolyFp<R,S> & f,
			    const UnivariatePolyFp<R,S> & g,
			    UnivariatePolyFp<R,S> & s,
			    UnivariatePolyFp<R,S> & t)
{
  UnivariatePolyFp<R,S> q(f.baseRing());
  UnivariatePolyFp<R,S> r(f.baseRing());
  UnivariatePolyFp<R,S> r_minus(f.baseRing());
  UnivariatePolyFp<R,S> r_plus(f.baseRing());
  UnivariatePolyFp<R,S> s_minus(f.baseRing());
  UnivariatePolyFp<R,S> s_plus(f.baseRing());
  UnivariatePolyFp<R,S> t_minus(f.baseRing());
  UnivariatePolyFp<R,S> t_plus(f.baseRing());
  
  FpElement<R,S> zero = f.baseRing()->zero();
  FpElement<R,S> one = f.baseRing()->one();
  s = zero;
  s_minus = one;
  t = one;
  t_minus = zero;
  
  r_minus = f;
  r = g;
  
  while (r != zero) {
    div_rem(r_minus, r, q, r_plus);

    assert(r_minus == q*r+r_plus);
    assert(s*f + t*g == r);
    assert(s_minus*f + t_minus*g == r_minus);

    r_minus = r;
    r = r_plus;
    s_plus = (s_minus - q*s);
    t_plus = (t_minus - q*t);
    s_minus = s;
    t_minus = t;
    s = s_plus;
    t = t_plus;
  }

  // finalize
  s = s_minus / r_minus.lead();
  t = t_minus / r_minus.lead();

  assert(r_minus / r_minus.lead() == s*f+t*g);
  
  return r_minus / r_minus.lead();

}

template<typename R, typename S>
inline void UnivariatePolyFp<R,S>::divRem(const UnivariatePolyFp<R,S> & f,
					  const UnivariatePolyFp<R,S> & g,
					  UnivariatePolyFp<R,S> & q,
					  UnivariatePolyFp<R,S> & r)
{

  assert(!g.isZero());

  UnivariatePolyFp<R,S> t(f.baseRing());
  q.makeZero();
  r = f;

  while ((!r.isZero()) && (r.degree() >= g.degree())) {
    FpElement<R,S> lc = r.lead() / g.lead();
    t = lc * UnivariatePolyFp<R,S>::x(f.baseRing(), r.degree()-g.degree());
    q += t;
    r -= t*g;
  }
  
  assert(f == q*g+r);

  return;
}
