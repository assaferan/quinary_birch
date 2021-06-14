#include <memory>
#include <set>
#include <vector>

#include "Fp.h"
#include "Integer.h"
#include "UnivariatePolyFp.h"

template<typename R>
template<typename S, typename T>
inline UnivariatePolyFp<S,T>
UnivariatePolyInt<R>::mod(std::shared_ptr< const Fp<S,T> > GF) const
{
  std::vector< FpElement<S,T> > vec;
  for (size_t i = 0; i < this->_coeffs.size(); i++)
    vec.push_back(GF->mod(this->_coeffs[i].num()));

  UnivariatePolyFp<S,T> ret(vec);
  
  return ret;
}

// We employ Yun's algorithm for this part
template<typename R>
inline std::vector< UnivariatePolyInt<R> >
UnivariatePolyInt<R>::_squarefreeFactor(void) const
{
  std::vector< UnivariatePolyInt<R> > fac;

  const UnivariatePolyInt<R> & f = *this;
  UnivariatePolyInt<R> f_prime = f.derivative();
  
  UnivariatePolyInt<R> a,b,c,d;
  b = f;
  c = f_prime;
  d = f_prime;
   
  while (b != Integer<R>::one()) {
    a = UnivariatePolyInt<R>::gcd(b,d);
    if (a == -Integer<R>::one())
      a = Integer<R>::one();
    fac.push_back(a);
    b /= a;
    c = d / a;
    d = c - b.derivative();
  }
  
  return fac;
}

// only makes sense for integers
template<typename R>
inline Integer<R> UnivariatePolyInt<R>::_landauMignotte(void) const
{
  Integer<R> d = R(this->degree() / 2);
  Integer<R> B = (d-Integer<R>::one()).binomialCoefficient(d/R(2)-Integer<R>::one());

  Integer<R> norm = this->_base->zero();
  for (size_t i = 0; i < this->_coeffs.size(); i++)
    norm += this->coefficient(i)*this->coefficient(i);

  // we might need ceiling here
  R sqrt_norm = sqrt(norm.num());
  norm = (d-Integer<R>::one()).binomialCoefficient(d/R(2))*sqrt_norm;

  return B + norm;
}

// Here we asume f = prod(u) mod p^i
// and sum((prod(u)/u_i) v_i = 1 mod p^i
// Want to lift to the same but mod p^(i+1)

template<typename R>
template<typename S, typename T>
inline void UnivariatePolyInt<R>::_henselStep(std::vector<UnivariatePolyInt<R> > & u,
					      std::vector<UnivariatePolyInt<R> > & v,
					      std::shared_ptr< const Fp<S,T> > GF,
					      size_t i) const
{
  Integer<R> p = birch_util::convert_Integer<S,R>(GF->prime());
  Integer<R> p_i = p^i;
  UnivariatePolyInt<R> prod = this->_base->one();
  for (size_t j = 0; j < u.size(); j++) {
    prod *= u[j];
  }
  UnivariatePolyInt<R> sum = this->_base->zero();
  
#ifdef DEBUG
  assert(this->lead() % p != 0);
  assert((this->lead() - u[0].lead()) % p_i == 0);
  assert(u.size() == v.size());
  UnivariatePolyInt<R> sum2 = this->_base->zero();
  for (size_t j = 0; j < u.size(); j++) {
    sum2 += (prod / u[j])*v[j];
    if (j > 0)
      assert(u[j].lead() == this->_base->one());
  }
  assert( ((*this)-prod) % p_i == 0);
#endif  

  // step 1 - lift the u_j
  
  u[0].lead() = this->lead();
  
  UnivariatePolyInt<R> t = ((*this) - prod) / p_i;
  UnivariatePolyInt<R> r;

  UnivariatePolyFp<S,T> t_p = t.mod(GF);
  UnivariatePolyFp<S,T> tv_bar(GF);
  UnivariatePolyFp<S,T> u_bar(GF);
  UnivariatePolyFp<S,T> q_bar(GF);
  UnivariatePolyFp<S,T> r_bar(GF);
  prod = this->_base->one();
  for (size_t j = 0; j < u.size(); j++) {
    u_bar = u[j].mod(GF);
    tv_bar = t_p*v[j].mod(GF);
    UnivariatePolyFp<S,T>::divRem(tv_bar, u_bar, q_bar, r_bar);
    r = r_bar.lift();
    u[j] += p_i * r;
    prod *= u[j];
  }

  return;
}



// here we lift all the way: f = prod(g) mod p
// and we lift to f = prod(g_lift) mod p^a

template<typename R>
template<typename S, typename T>
inline std::vector< UnivariatePolyInt<R> >
UnivariatePolyInt<R>::_henselLift(const std::vector<UnivariatePolyFp<S,T> > & g,
				  size_t a) const
{
  R p = g[0].baseRing()->prime();
  std::vector< UnivariatePolyInt<R> > u, v;
  std::vector< UnivariatePolyFp<S,T> > v_bar;
  UnivariatePolyFp<S,T> t(g[0].baseRing());
  UnivariatePolyFp<S,T> prod = g[0];
  FpElement<S,T> one = g[0].baseRing()->one();
  
  if (g.size() == 1) {
    u.push_back(*this);
    return u;
  }
 
  v_bar.push_back(one);
  for (size_t i = 1; i < g.size(); i++) {
    // we just push something so we will be able to use it
    v_bar.push_back(one);
    UnivariatePolyFp<S,T>::xgcd(prod,g[i],v_bar[i],t);
    for (size_t j = 0; j < i; j++)
      v_bar[j] *= t;
    prod *= g[i];
  }
#ifdef DEBUG
  UnivariatePolyFp<S,T> s = -one;
  for (size_t i = 0; i < g.size(); i++)
    s += (prod / g[i]) * v_bar[i];
  assert(s.isZero());
#endif
  v.resize(g.size());
  u.resize(g.size());
  for (size_t i = 0; i < g.size(); i++) {
    v[i] = v_bar[i].lift();
    u[i] = g[i].lift();
  }

  for (size_t i = 1; i < a; i++) {
    this->_henselStep(u, v, g[0].baseRing(), i);
  }

  return u;
}


template<typename R>
inline std::set< std::set<size_t> >
UnivariatePolyInt<R>::_subsets(const std::set<size_t> & S, size_t k)
{
  std::set< std::set<size_t> > subs;
  if (k > S.size())
    return subs;
  
  if (k == 0) {
    std::set<size_t> emptyset;
    subs.insert(emptyset);
    return subs;
  }
  
  std::set<size_t> S_copy = S;
  size_t i = *S_copy.begin();
  S_copy.erase(i);
  
  subs = UnivariatePolyInt<R>::_subsets(S_copy, k);
  std::set< std::set<size_t> > subs_i = UnivariatePolyInt<R>::_subsets(S_copy, k-1);

  for (std::set<size_t> sub : subs_i) {
    std::set<size_t> sub_i = sub;
    sub_i.insert(i);
    subs.insert(sub_i);
  }

  return subs;
}

template<typename R>
inline R _balance(const R & a, const R & n)
{
  R b = a % n;
  if ((b+b) > n)
    b -= n;
  return b;
}

template<typename R>
inline void
UnivariatePolyInt<R>::_findTrialFactor(const std::vector< UnivariatePolyInt<R> > & u,
				       const Integer<R> & N,
				       size_t & j,
				       std::set<size_t> & C,
				       size_t & s,
				       std::vector< UnivariatePolyInt<R> > & gs )
{
  UnivariatePolyInt<R> g, q, r;
  for (size_t m = j; m <= C.size(); m++) {
    std::set< std::set<size_t> > subs = UnivariatePolyInt<R>::_subsets(C,m);
    for (std::set<size_t> sub : subs) {
      g = this->lead();
      for (size_t i : sub)
	g *= u[i];
      for (size_t i = 0; i < g._coeffs.size(); i++)
	g._coeffs[i] = _balance(g._coeffs[i], N);
      UnivariatePolyInt<R>::divRem((*this), g, q, r);
      if (r.isZero()) {
	s++;
	gs.push_back(g);
	*this = q;
	j = m;
	for (size_t i : sub)
	  C.erase(i);
	return;
      }
    }
  }
  return;
}


template<typename R>
inline std::vector< UnivariatePolyInt<R> >
UnivariatePolyInt<R>::_trialFactor(const std::vector< UnivariatePolyInt<R> > & u,
				   const Integer<R> & N) const
{
  UnivariatePolyInt<R> h = *this;
  size_t r = u.size();
  std::set<size_t> C;
  for (size_t i = 1; i < r; i++)
    C.insert(i);
  size_t j,s,t;

  s = 0;
  j = 1;
  
  std::vector< UnivariatePolyInt<R> > g;

  do {
    t = s;
    h._findTrialFactor(u,N,j,C,s,g);
  } while (t != s);

  g.push_back(h);

  return g;
}
