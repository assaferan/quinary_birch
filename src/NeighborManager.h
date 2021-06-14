#ifndef __NEIGHBOR_MANAGER_H_
#define __NEIGHBOR_MANAGER_H_

#include "birch.h"
#include "Polynomial.h"
#include "QuadForm.h"
#include "Isometry.h"
#include "Fp.h"

template<typename R, typename S, typename T, size_t n>
class NeighborManager
{
public:
  NeighborManager(const QuadFormZZ<T,n>& q, std::shared_ptr<Fp<R,S>> GF,
		  size_t k = 1);

  void nextIsotropicSubspace(void);

  GenusRep<T,n> getReducedNeighborRep(void);

  // to representative of the line
  VectorInt<R,n> transformVector(const GenusRep<T,n>& dst, VectorInt<R,n> src);
  
  void getNextNeighbor(void);

  QuadFormZZ<T,n> buildNeighbor(Isometry<T,n>& ) const;

  const std::vector< VectorInt<T,n> > & getIsotropicSubspace() const
  {return this->X; }

protected:
  std::shared_ptr< Fp<R,S> > _GF;
  QuadFormZZ<T,n> _q;
  T _disc;
  SquareMatrixFp<R,S,n> _b;
  SquareMatrixInt<T,n> _quot_gram;
  std::shared_ptr< SquareMatrixFp<R,S,n> > _p_std_gram;
  std::shared_ptr< SquareMatrixFp<R,S,n> > _p_basis;
  std::shared_ptr< PolynomialFp<R,S> > _p_q_std;
  // dimension of the radical
  size_t _rad_dim;
  // dimension of the anisotropic subspace
  size_t _aniso_dim;
  // the Witt index (number of hyperbolic planes)
  size_t _witt_index;

  VectorFp<R,S,n> _vec;
  std::vector< std::vector< size_t> > _pivots;
  size_t _pivot_ptr;
  size_t _k; // dimension of the isotropic subspace
  size_t _skew_dim;
  std::shared_ptr< MatrixFp<R,S> > _p_skew;
  std::vector<size_t> _free_vars;
  std::vector<FpElement<R,S> > _params;
  std::shared_ptr<Matrix<PolynomialFp<R,S>, PolynomialRingFp<R,S> > > _p_isotropic_param;
  std::vector< VectorFp<R,S,n> > _iso_subspace;
  std::vector< VectorInt<T,n> > _X, _Z, _U;

  // The 2-isotropic vectors were stored in binary within each of the
  // coordinates of `vec` and so we use this function to unpack them into
  // actual 2-isotropic vectors.
  VectorInt<R,n> _isotropicVector_p2(R t) const;

  // get all possible pivots
  static std::vector< std::vector<size_t> >
  __pivots(size_t dim, size_t aniso, size_t k);

  void __initializePivots(void);

  SquareMatrixInt<T,n> __gram(const SquareMatrixInt<T,n> & B, bool quot = true) const;
  
  void _liftSubspace();
  void _updateSkewSpace();
  void _updateSkewMatrix(size_t &, size_t &);
};

#include "NeighborManager.inl"

#endif // __NEIGHBOR_MANAGER_H
