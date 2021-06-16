#ifndef __SPINOR_H_
#define __SPINOR_H_

#include <random>
#include <unordered_map>

#include "birch.h"
#include "NeighborManager.h"
#include "SquareMatrix.h"

template<typename R>
class Spinor
{
  template<typename T>
  friend class Spinor;

public:
  template <size_t n>
  Spinor(const std::vector<R>& primes, const QuadFormZZ<R,n>& q)
  {
    this->_primes = primes;
    this->_twist = (1LL << this->_primes.size()) - 1;
    std::random_device rd;
    this->_seed = rd();
    for (R prime : primes) {
      std::shared_ptr<W16_Fp> GF = std::make_shared<W16_Fp>(prime, this->_seed, true);
      NeighborManager<W16,W32,R,n> manager(q, GF);
      std::vector< W16_VectorFp<n> > rad = manager.radical();
      W16_MatrixFp rad_mat(GF, rad.size(), n);
      for (size_t row = 0; row < rad.size(); row++)
	for (size_t col = 0; col < n; col++)
	  rad_mat(row,col) = rad[row][col];
      this->_rads.insert(std::make_pair(prime,rad_mat));
    }
  }

  template<size_t n>
  inline Z64 norm(const Isometry<R,n>& s) const
  {
    std::vector< W16_FpElement> dets;
    for (R prime : this->_primes) {
      std::shared_ptr<W16_Fp> GF = std::make_shared<W16_Fp>(prime, this->_seed, true);
      std::shared_ptr< SquareMatrixFp<W16,W32,n> > s_p = mod(s.integralMatrix(),GF);
      W16_MatrixFp s_mat(*s_p);
      W16_MatrixFp rad = this->_rads.at(prime);
      W16_FpElement det = s_mat.transpose().restrict(rad).determinant();
      dets.push_back(det);
    }
    return this->_computeVals(dets);    
  }

  inline const std::vector<R> & primes(void) const
  {
    return this->_primes;
  }

private:
  std::vector<R> _primes;
  Z64 _twist;
  W64 _seed;
  std::unordered_map<R, W16_MatrixFp> _rads;

  inline Z64 _computeVals(const std::vector< W16_FpElement> & a) const
  {
    Z64 val = 0;
    Z64 mask = 1;
    for (size_t i = 0; i < a.size(); i++) {
      if (!a[i].isOne()) {
	assert ((-a[i]).isOne());
	val ^= mask;
      }
      mask <<= 1;
    }
    return val;
  }
};

#endif // __SPINOR_H_
