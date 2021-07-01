#include <cassert>
#include <memory>
#include <random>

#include "birch_util.h"

template<typename R, typename S>
Fp<R,S>::Fp(const R& p, W64 seed, bool use_inverse_lut)
{
  std::random_device rd;
  this->_rng = std::unique_ptr<std::mt19937>( new std::mt19937(seed) );
  this->_distr = std::unique_ptr<std::uniform_int_distribution<>>
    (new std::uniform_int_distribution<>(0, p-1));

  this->_p = p;
  if (this->_p != 2)
    {
      this->_kp = (((R)-1)/p)*p;
      this->_kp_inv = ((S)-1)/(this->_kp);
      this->_use_inverse_lut = use_inverse_lut;
      if (use_inverse_lut) this->_makeInverseLut();
    }
}

template<typename R, typename S>
template<typename T>
inline FpElement<R,S> Fp<R,S>::mod(const T& a) const
{
  // we relax this for now to enable Z128 support
  //        static_assert(std::is_integral<T>::value, "Undefined type.");
  // std::cerr << "a = " << a << std::endl;
  T value = (T)a % this->_p;
  // std::cerr << "value = " << value << std::endl;
  R r_val = (value < 0) ? (R)(value+this->_p) : (R)value;
  // std::cerr << "r_val = " << r_val << std::endl;
  return FpElement<R,S>(this->getPtr(), r_val);
}

template<typename R, typename S>
std::shared_ptr< const Fp<R, S> > Fp<R,S>::getPtr() const {
  return std::enable_shared_from_this< const Fp<R, S> >::shared_from_this();
}

template<typename R, typename S>
inline R Fp<R,S>::neg(R a) const
{
  return this->_kp-a;
}

template<typename R, typename S>
inline R Fp<R,S>::mul(R a, R b) const
{
  S rem = ((S)a)*b;
  R hi = rem >> _bits;
  R t = (((S)hi*(S)this->_kp_inv) >> _bits) + hi;
  rem -= (S)t * this->_kp;
  rem = (rem >= this->_kp) ? rem-this->_kp : rem;
  rem = (rem >= this->_kp) ? rem-this->_kp : rem;
  rem = (rem >= this->_kp) ? rem-this->_kp : rem;

  assert( ((S)a*(S)b)%(this->_p) == rem%(this->_p) );

  return rem;
}

template<typename R, typename S>
inline R Fp<R,S>::add(R a, R b) const
{
  R neg = this->_kp-a;
  return (b >= neg) ? b-neg : this->_kp-(neg-b);
}

template<typename R, typename S>
inline R Fp<R,S>::sub(R a, R b) const
{
  return add(a, this->_kp-b);
}

template<typename R, typename S>
inline int Fp<R,S>::legendre(R a) const
{
  Z aa = birch_util::convertInteger<R,Z>(a);
  Z pp = birch_util::convertInteger<R,Z>(this->_p);
  return mpz_legendre(aa.get_mpz_t(), pp.get_mpz_t());
}

template<typename R, typename S>
inline R Fp<R,S>::inverse(R a) const
{
  if (this->_use_inverse_lut) return this->_inverse_lut[a];
  else return this->_inv(a);
}

template<typename R, typename S>
inline R Fp<R,S>::inverse(const Z& a) const
{
  R inv = mpz_get_ui(a.get_mpz_t());
  return this->inverse(inv);
}

template<typename R, typename S>
inline R Fp<R,S>::inverse(const Z64& a) const
{
  R inv = (R)a;
  return this->inverse(inv);
}

template<typename R, typename S>
inline FpElement<R, S> Fp<R,S>::random(void) const
{
  return FpElement<R,S>(this->getPtr(), (R)(*this->_distr)(*this->_rng));
}

template<typename R, typename S>
inline R Fp<R,S>::_inv(R a) const
{
  if (a == 0) return 0;
  Z aa = birch_util::convertInteger<R,Z>(a);
  Z pp = birch_util::convertInteger<R,Z>(this->_p);
  mpz_invert(aa.get_mpz_t(), aa.get_mpz_t(), pp.get_mpz_t());
  R ainv = mpz_get_ui(aa.get_mpz_t());

  assert( ((S)a * (S)ainv) % (this->_p) == 1 );

  return ainv;
}

template<typename R, typename S>
inline void Fp<R,S>::_inverseLutPopulate(Z32 offset, Z32 len)
{
  Z32 datalen = len<<1;

  // Set the initial values.
  std::iota(this->_inverse_lut.begin()+offset,
	    this->_inverse_lut.begin()+offset+len, offset);

  // Multiply up to the root node...
  for (Z32 i=0, j=len; i<j; i+=2, j++)
    {
      R a = this->_inverse_lut[offset+i];
      R b = this->_inverse_lut[offset+i+1];
      this->_inverse_lut[offset+j] = this->mul(a, b);
    }

  // ...invert the root node...
  R ainv = this->_inv(this->_inverse_lut[offset+datalen-2]);
  this->_inverse_lut[offset+datalen-2] = ainv;

  // ...and then backtrack to the inverse.
  for (Z32 i=datalen-4, j=datalen-2; i>=0; i-=2, --j)
    {
      R temp = this->_inverse_lut[offset+i];
      R a = this->_inverse_lut[offset+i+1];
      R b = this->_inverse_lut[offset+j];
      this->_inverse_lut[offset+i] = this->mul(a, b);
      this->_inverse_lut[offset+i+1] = this->mul(temp, b);
    }
}

template<typename R, typename S>
inline void Fp<R,S>::_makeInverseLut(void)
{
  // TODO: The following could probably be replaced with clz.
  Z32 len = 1;
  R q = this->_p;
  while (q != 1)
    {
      len <<= 1;
      q >>= 1;
    }

  // Make enough room in the lut to grow the tree.
  _inverse_lut.resize(1+(len<<1));

  Z32 offset = 1;
  q = this->_p;
  while (q > 1)
    {
      this->_inverseLutPopulate(offset, len);

      offset |= len;
      q ^= len;

      R r = q;
      len = 1;
      while (r > 1)
	{
	  len <<= 1;
	  r >>= 1;
	}
    }

  // Shrink the lut to the appropriate size.
  this->_inverse_lut.resize(this->_p);

  for (R i=1; i<(this->_p); i++)
    {
      assert( this->mul(i, this->_inverse_lut[i]) % (this->_p) == 1 );
    }
}
