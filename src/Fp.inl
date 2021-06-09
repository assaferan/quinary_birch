#include <random>

template<typename R, typename S>
Fp<R,S>::Fp(const R& p, W64 seed, bool use_inverse_lut=false)
{
  std::random_device rd;
  this->rng = std::unique_ptr<std::mt19937>( new std::mt19937(seed) );
  this->distr = std::unique_ptr<std::uniform_int_distribution<>>(0, p-1);

  this->p = p;
  if (this->p != 2)
    {
      this->kp = (((R)-1)/p)*p;
      this->kp_inv = ((S)-1)/kp;
      this->use_inverse_lut = use_inverse_lut;
      if (use_inverse_lut) this->make_inverse_lut();
    }
}

template<typename R, typename S>
template<typename T>
inline FpElement<R,S> Fp<R,S>::mod(const T& a) const
{
  // we relax this for now to enable Z128 support
  //        static_assert(std::is_integral<T>::value, "Undefined type.");
  T value = (T)a % this->p;
  R r_val = (value < 0) ? (R)(value+this->p) : (R)value;
  return FpElement<R,S>(this->getptr(), r_val);
}

template<typename R, typename S>
std::shared_ptr< const Fp<R, S> > Fp<R,S>::getptr() const {
  return std::enable_shared_from_this< const Fp<R, S> >::shared_from_this();
}

template<typename R, typename S>
inline virtual R Fp<R,S>::neg(R a) const
{
  return this->kp-a;
}

template<typename R, typename S>
inline virtual R Fp<R,S>::mul(R a, R b) const
{
  S rem = ((S)a)*b;
  R hi = rem >> bits;
  R t = (((S)hi*(S)this->kp_inv) >> bits) + hi;
  rem -= (S)t * this->kp;
  rem = (rem >= this->kp) ? rem-this->kp : rem;
  rem = (rem >= this->kp) ? rem-this->kp : rem;
  rem = (rem >= this->kp) ? rem-this->kp : rem;

  assert( ((S)a*(S)b)%p == rem%p );

  return rem;
}

template<typename R, typename S>
inline virtual R Fp<R,S>::add(R a, R b) const
{
  R neg = this->kp-a;
  return (b >= neg) ? b-neg : this->kp-(neg-b);
}

template<typename R, typename S>
inline virtual R Fp<R,S>::sub(R a, R b) const
{
  return add(a, this->kp-b);
}

template<typename R, typename S>
inline int Fp<R,S>::legendre(R a) const
{
  Z aa(a);
  Z pp(p);
  return mpz_legendre(aa.get_mpz_t(), pp.get_mpz_t());
}

template<typename R, typename S>
inline virtual R Fp<R,S>::sqrt(R a) const
{
  a = a % p;
  if (a == 1) return 1;
  if (a == 0) return 0;
  if (this->legendre(a) != 1) return 0;

  R q = p-1;
  R s = 0;
  while(q % 2 == 0) { q >>= 1; s++; }

  if (s == 1) return this->pow(a, (p+1)/4);

  R z = this->random().lift();
  while (z == 0 || this->legendre(z) == 1)
    {
      z = this->random().lift();
    }

  int m = s;
  R c = this->pow(z, q);
  R r = this->pow(a, (q+1)/2);
  R t = this->pow(a, q);
  if (t >= p) t = this->mod(t).lift();

  while (1)
    {
      if (t == 1) return r;
      int i = 0;
      R t1 = t;
      while (t1 != 1)
	{
	  t1 = this->mul(t1, t1);
	  if (t1 >= p) t1 = this->mod(t1).lift();
	  i++;
	}

      int e = 1;
      for (int j=0; j<m-i-1; j++) e <<= 1;
      R b = this->pow(c, e);
      r = this->mul(r, b);
      c = this->mul(b, b);
      t = this->mul(t, c);
      if (t >= p) t %= p;
      m = i;
    }

  return 0;
}

template<typename R, typename S>
inline virtual R Fp<R,S>::inverse(R a) const
{
  if (this->use_inverse_lut) return this->inverse_lut[a];
  else return this->inv(a);
}

template<typename R, typename S>
inline virtual R Fp<R,S>::inverse(const Z& a) const
{
  R inv = mpz_get_ui(a.get_mpz_t());
  return this->inverse(inv);
}

template<typename R, typename S>
inline virtual R Fp<R,S>::inverse(const Z64& a) const
{
  R inv = (R)a;
  return this->inverse(inv);
}

template<typename R, typename S>
inline FpElement<R, S> Fp<R,S>::random(void) const
{
  return FpElement<R,S>(this->getptr(), (R)(*this->distr)(*this->rng));
}

template<typename R, typename S>
inline virtual R Fp<R,S>::inv(R a) const
{
  if (a == 0) return 0;
  Z aa(a);
  Z pp(p);
  mpz_invert(aa.get_mpz_t(), aa.get_mpz_t(), pp.get_mpz_t());
  R ainv = mpz_get_ui(aa.get_mpz_t());

  assert( ((S)a * (S)ainv) % p == 1 );

  return ainv;
}

template<typename R, typename S>
void Fp<R,S>::inverseLutPopulate(Z32 offset, Z32 len)
{
  Z32 datalen = len<<1;

  // Set the initial values.
  std::iota(this->inverse_lut.begin()+offset,
	    this->inverse_lut.begin()+offset+len, offset);

  // Multiply up to the root node...
  for (Z32 i=0, j=len; i<j; i+=2, j++)
    {
      R a = this->inverse_lut[offset+i];
      R b = this->inverse_lut[offset+i+1];
      this->inverse_lut[offset+j] = this->mul(a, b);
    }

  // ...invert the root node...
  R ainv = this->inv(this->inverse_lut[offset+datalen-2]);
  this->inverse_lut[offset+datalen-2] = ainv;

  // ...and then backtrack to the inverse.
  for (Z32 i=datalen-4, j=datalen-2; i>=0; i-=2, --j)
    {
      R temp = this->inverse_lut[offset+i];
      R a = this->inverse_lut[offset+i+1];
      R b = this->inverse_lut[offset+j];
      this->inverse_lut[offset+i] = this->mul(a, b);
      this->inverse_lut[offset+i+1] = this->mul(temp, b);
    }
}

template<typename R, typename S>
void Fp<R,S>::makeInverseLut(void)
{
  // TODO: The following could probably be replaced with clz.
  Z32 len = 1;
  R q = p;
  while (q != 1)
    {
      len <<= 1;
      q >>= 1;
    }

  // Make enough room in the lut to grow the tree.
  inverse_lut.resize(1+(len<<1));

  Z32 offset = 1;
  q = p;
  while (q > 1)
    {
      this->inverse_lut_populate(offset, len);

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
  this->inverse_lut.resize(p);

  for (Z32 i=1; i<p; i++)
    {
      assert( this->mul(i, this->inverse_lut[i]) % p == 1 );
    }
}
