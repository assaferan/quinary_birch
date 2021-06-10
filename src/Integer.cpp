#include "birch_util.h"
#include "Integer.h"

template class Integer<Z>;
template class Integer<Z64>;
template class Integer<Z128>;
template class Integer<W16>;

const std::vector<int> hilbert_lut_odd = { 1, 1, 1, 1,
                                           1, 1,-1,-1,
                                           1,-1, 1,-1,
                                           1,-1,-1, 1 };

const std::vector<int> hilbert_lut_p2 = { 1, 1, 1, 1, 1, 1, 1, 1,
                                          1,-1, 1,-1,-1, 1,-1, 1,
                                          1, 1, 1, 1,-1,-1,-1,-1,
                                          1,-1, 1,-1, 1,-1, 1,-1,
                                          1,-1,-1, 1, 1,-1,-1, 1,
                                          1, 1,-1,-1,-1,-1, 1, 1,
                                          1,-1,-1, 1,-1, 1, 1,-1,
                                          1, 1,-1,-1, 1, 1,-1,-1 };

int hilbertSymbol(const Z & x, const Z & y, const Z & p) {
  Z a = x;
  Z b = y;
  int a_val = 0;
  while (a % p == 0)
    {
      ++a_val;
      a /= p;
    }

  int b_val = 0;
  while (b % p == 0)
    {
      ++b_val;
      b /= p;
    }

  if (p == 2)
    {
      int aa = (mpz_class(a%8).get_si() >> 1) & 0x3;
      int bb = (mpz_class(b%8).get_si() >> 1) & 0x3;
      int index = ((a_val&0x1)<<5) | (aa << 3) | ((b_val&0x1)<<2) | bb;
      return hilbert_lut_p2[index];
    }

  int a_notsqr = mpz_legendre(a.get_mpz_t(), p.get_mpz_t()) == -1;
  int b_notsqr = mpz_legendre(b.get_mpz_t(), p.get_mpz_t()) == -1;

  int index = ((a_val&0x1)<<3) | (a_notsqr<<2) | ((b_val&0x1)<<1) | b_notsqr;
  if (((index & 0xa) == 0xa) && ((p%4) == 0x3))
    {
      return -hilbert_lut_odd[index];
    }
  else
    {
      return hilbert_lut_odd[index];
    }
}

template<>
int Integer<Z>::hilbertSymbol(const Integer<Z>& other, const Integer<Z>& p) const
{
  return hilbertSymbol(this->_num, other._num, p._num);
}


template<>
int Integer<Z64>::hilbertSymbol(const Integer<Z64> b, const Integer<Z64>& p)
{
  return hilbertSymbol(birch_util::convert_Integer<Z64,Z>(this->_num),
		       birch_util::convert_Integer<Z64,Z>(b._num),
		       birch_util::convert_Integer<Z64,Z>(p._num));
}

template<>
int Integer<Z128>::hilbert_symbol(const Integer<Z128> b, const Integer<Z128>& p)
{
  return hilbertSymbol(birch_util::convert_Integer<Z128,Z>(this->_num),
		       birch_util::convert_Integer<Z128,Z>(b._num),
		       birch_util::convert_Integer<Z128,Z>(p._num));
}
