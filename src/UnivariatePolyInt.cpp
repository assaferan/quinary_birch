#include "Integer.h"
#include "UnivariatePolyInt.h"

template class UnivariatePolyInt<Z>;

template<>
W64 UnivariatePolyInt<Z>::hashValue(void) const
{
  W64 fnv = FNV_OFFSET;
  for (size_t i = 0; i < this->coefficients().size(); i++)
    fnv = (fnv ^ mpz_get_si(this->coefficient(i).num().get_mpz_t())) * FNV_PRIME;
  return fnv;
}

// !! TODO - now that we have PolyInt, we could write it template in R
template<>
std::unordered_map< UnivariatePolyInt<Z>, size_t >
UnivariatePolyInt<Z>::factor(void) const
{
  std::unordered_map< UnivariatePolyInt<Z>, size_t > fac;

  std::vector< UnivariatePolyInt<Z> > sqf = this->_squarefreeFactor();
  std::random_device rd;
  W64 seed = rd();
  
  for (size_t i = 0; i < sqf.size(); i++) {
    UnivariatePolyInt<Z> f = sqf[i];
    if (f == Integer<Z>::one()) continue;
    
    Integer<Z> p = Z(3);
    W16 p_16 = birch_util::convertInteger<Z,W16>(p.num());
    std::shared_ptr< const W16_Fp > GF
      = std::make_shared< W16_Fp >(p_16,seed);
    UnivariatePolyFp<W16,W32> f_p = f.mod(GF);
    UnivariatePolyFp<W16,W32> d =
      UnivariatePolyFp<W16,W32>::gcd(f_p, f_p.derivative());
    while (d.degree() > 0) {
      p = p.nextPrime();
      p_16 = birch_util::convertInteger<Z,W16>(p.num());
      GF = std::make_shared<W16_Fp>(p_16,seed);
      f_p = f.mod(GF);
      d = UnivariatePolyFp<W16,W32>::gcd(f_p, f_p.derivative());
    }
    
    std::vector< UnivariatePolyFp<W16,W32> > fac_p = f_p._sqfFactor();
    Integer<Z> L = f._landauMignotte();
    size_t a = 1;
    Integer<Z> p_a = p;
    while (p_a <= (L+L)) {
      a++;
      p_a *= p;
    }
    std::vector< UnivariatePolyInt<Z> > fac_lift = f._henselLift(fac_p,a);

    fac_lift = f._trialFactor(fac_lift, p_a);
    
    for ( UnivariatePolyInt<Z> g : fac_lift) {
      fac[g] = i;
    }
  }
  
  return fac;
}
