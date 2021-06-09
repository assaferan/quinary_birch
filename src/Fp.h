#ifndef __FP_H_
#define __FP_H_

#include <memory>
#include <random>

#include "birch.h"

template<typename R, typename S>
class FpElement;

template<typename R, typename S>
class Fp : public std::enable_shared_from_this< const Fp<R, S> >
{
public:
  Fp(const R& p, W64 seed, bool use_inverse_lut=false);

  inline const R& prime(void) const { return this->p; }

  template<typename T>
  inline FpElement<R,S> mod(const T& a) const;
    
  inline std::shared_ptr< const Fp<R, S> > getptr() const;

  inline virtual R neg(R a) const;
  inline virtual R mul(R a, R b) const;
  inline virtual R add(R a, R b) const;
  inline virtual R sub(R a, R b) const;

  inline int legendre(R a) const;

  inline virtual R sqrt(R a) const;

  inline virtual R inverse(R a) const;

  inline virtual R inverse(const Z& a) const;

  inline virtual R inverse(const Z64& a) const;

  inline FpElement<R, S> random(void);

private:
  R p;
  R kp;
  R kp_inv;
  static constexpr int bits = 8 * sizeof(R);
  bool use_inverse_lut;
  std::vector<R> inverse_lut;

  // Random number generator.
  std::unique_ptr<std::mt19937> rng;
  std::unique_ptr<std::uniform_int_distribution<>> distr;

  inline virtual R inv(R a) const;
  void inverseLutPopulate(Z32 offset, Z32 len);

  void makeInverseLut(void);
};

template<typename R, typename S>
class F2 : public Fp<R,S>
{
public:
  F2(const R& p, W64 seed) : Fp<R,S>(p, seed, false) {}

  inline R mul(R a, R b) const override
  { return ((a & b) & 1); }

  inline R add(R a, R b) const override
  { return ((a ^ b) & 1); }

  inline R sub(R a, R b) const override
  { return ((a ^ b) & 1);}

  inline R pow(R a, Z64 e) const override
  { return e == 0 ? 1 : (a & 1);}

  inline R sqrt(R a) const override
  { return (a & 1); }

  inline R inverse(R a) const override
  { return (a & 1);}

  inline R inverse(const Z& a) const override
  {return (mpz_get_ui(a.get_mpz_t()) & 1); }

  inline R inverse(const Z64& a) const override
  { return (a & 1);}

  inline R neg(R a) const override
  { return (a & 1);}
  
private:
  inline R inv(R a) const override
  { return (a & 1); }
};

template<>
template<>
FpElement<W16, W32> W16_Fp::mod(const Z& a) const;

template<>
template<>
FpElement<W32, W64> W32_Fp::mod(const Z& a) const;

template<>
template<>
FpElement<W64, W128> W64_Fp::mod(const Z& a) const;

#include "Fp.inl"

#endif // __FP_H_
