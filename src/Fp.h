#ifndef __FP_H_
#define __FP_H_

#include <memory>
#include <random>

#include "birch.h"
#include "Ring.h"

template<typename R, typename S>
class FpElement;

template<typename R, typename S>
class Fp : public virtual Ring< Fp<R,S>, FpElement<R,S> >
{
public:
  Fp(const R& p, W64 seed, bool use_inverse_lut=false);

  inline const R& prime(void) const { return this->p; }

  template<typename T>
  FpElement<R,S> mod(const T& a) const;
    
  std::shared_ptr< const Fp<R, S> > getPtr() const override;

  virtual R neg(R a) const;
  virtual R mul(R a, R b) const;
  virtual R add(R a, R b) const;
  virtual R sub(R a, R b) const;

  int legendre(R a) const;

  virtual R inverse(R a) const;

  virtual R inverse(const Z& a) const;

  virtual R inverse(const Z64& a) const;

  FpElement<R,S> random(void);

  inline FpElement<R,S> zero(void) const override
  {return FpElement<R,S>::zero(this->getPtr()); }
  
  inline FpElement<R,S> one(void) const override
  {return FpElement<R,S>::one(this->getPtr()); }

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

  virtual R inv(R a) const;
  void inverseLutPopulate(Z32 offset, Z32 len);

  void makeInverseLut(void);
};

template<typename R, typename S>
class F2 : public virtual Fp<R,S>
{
public:
  F2(const R& p, W64 seed) : Fp<R,S>(p, seed, false) {}

  inline R mul(R a, R b) const override
  { return ((a & b) & 1); }

  inline R add(R a, R b) const override
  { return ((a ^ b) & 1); }

  inline R sub(R a, R b) const override
  { return ((a ^ b) & 1);}

  inline R pow(R a, Z64 e) const
  { return e == 0 ? 1 : (a & 1);}

  inline R sqrt(R a) const
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
