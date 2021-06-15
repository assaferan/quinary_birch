#ifndef __GENUS_H_
#define __GENUS_H_

#include <unordered_map>

#include "birch.h"
#include "Eigenvector.h"
#include "HashMap.h"
#include "Isometry.h"
#include "NeighborManager.h"
#include "NumberField.h"
#include "Rational.h"
#include "Spinor.h"

template<typename R, size_t n>
class GenusRep
{
public:
  GenusRep() = default;
  GenusRep(const GenusRep<R, n>& genus) = default;
  GenusRep(GenusRep<R, n>&& genus) = default;

  // comparison
  inline bool operator==(const GenusRep<R,n>& other) const
  { return this->q == other.q; }
  
  QuadFormZZ<R,n> q;
  Isometry<R,n> s;
  Isometry<R,n> sinv;
  Z64 parent;
  R p;
  std::map<R,int> es;
};

template<typename R, size_t n>
class Genus
{
  template<typename T, size_t N>
  friend class Genus;

public:
  // c-tors
  Genus() = default;

  Genus(const QuadFormZZ<R,n>& q,
	const std::vector<PrimeSymbol<R>>& symbols, W64 seed=0);

  // copy c-tor
  template<typename T>
  Genus(const Genus<T,n>& src);

  template<typename T>
  inline static Genus<T,n> convert(const Genus<R,n>& src)
  { return Genus<T,n>(src); }

  // access
  inline size_t size(void) const
  { return this->_hash->keys().size(); }

  inline W64 seed(void) const
  { return this->_seed; }

  inline std::map<R,size_t> dimensionMap(void) const
  {
    std::map<R,size_t> temp;
    size_t num_conductors = this->_conductors.size();
    for (size_t k=0; k<num_conductors; k++)
      {
	temp[this->_conductors[k]] = this->_dims[k];
      }
    return temp;
  }

  inline std::map<R,std::vector<int>> heckeMatrixDense(const R& p) const
  {
    if (this->_disc % p == 0)
      {
	throw std::invalid_argument("Prime must not divide the discriminant.");
      }
    return this->_heckeMatrixDenseInternal(p);
  }

  inline std::map<R,std::vector<std::vector<int>>> heckeMatrixSparse(const R& p) const
  {
    if (this->_disc % p == 0)
      {
	throw std::invalid_argument("Prime must not divide the discriminant.");
      }
    return this->_heckeMatrixSparseInternal(p);
  }

  Eigenvector<R> eigenvector(const std::vector<Z32>&, const R& ) const;

  std::vector<Z32> eigenvalues(EigenvectorManager<R,n>&, const R&) const;

  inline const GenusRep<R,n>& representative(size_t idx) const
  { return this->_hash->get(idx); }

  inline size_t indexof(const GenusRep<R,n>& rep) const
  { return this->_hash->indexof(rep); }

  std::map<R, std::vector< std::vector< NumberFieldElement<Z> > > >
  eigenvectors();

protected:
  R _disc;
  std::vector<R> _prime_divisors;
  std::vector<R> _conductors;
  std::vector<size_t> _dims;
  std::vector<std::vector<size_t>> _num_auts;
  std::vector<std::vector<int>> _lut_positions;
  Rational<Z> _mass;
  std::unique_ptr<HashMap<W16>> _spinor_primes;
  std::unique_ptr<HashMap<GenusRep<R, n>>> _hash;
  std::unique_ptr<Spinor<R>> _spinor;
  // hashing the invariants
  std::unique_ptr<HashMap<GenusRep<R, n>>> _inv_hash;
  // mapping from the invariant hash to the genus representatives
  std::unordered_map< size_t, size_t > _inv_map;
  W64 _seed;

  template<typename S, typename T>
  std::vector<Z32> _eigenvectors(EigenvectorManager<R, n>&,
				 std::shared_ptr<Fp<S,T>>, const R& ) const;

  Rational<Z> _getMass(const QuadFormZZ<R,n>&,
		       const std::vector<PrimeSymbol<R>>&);

  static Rational<Z> _localFactor(const MatrixRat<R> & g,
				  const Integer<R> & p);

  static Rational<Z> _combine(const QuadFormZZ<R,n>& q,
			      const Integer<R> & p);

  std::map<R,std::vector<std::vector<int>>>
  _heckeMatrixSparseInternal(const R& ) const;

  std::map<R,std::vector<int>> _heckeMatrixDenseInternal(const R&) const;

  static std::set<Integer<R> > _wittToHasse(const Integer<R> &,
					    const std::set<std::pair<Integer<R>, int> > &);

  std::vector< MatrixInt<int> > _decomposition(size_t k) const;
  
  std::vector< MatrixInt<int> >
  _decompositionRecurse(const MatrixInt<int> & V_basis,
			const Integer<R> & p, size_t k) const;
};

namespace std
{
  template<typename R, size_t n>
  struct hash<GenusRep<R,n>>
  {
    Z64 operator()(const GenusRep<R,n>& rep) const
    {
      return rep.q.hashValue();
    }
  };
}

#include "Genus.inl"

#endif // __GENUS_H_
