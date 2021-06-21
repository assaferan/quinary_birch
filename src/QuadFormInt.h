#ifndef __QUAD_FORM_INT_H_
#define __QUAD_FORM_INT_H_

#include <memory>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "birch.h"
#include "Integer.h"
#include "ParseNipp.h"
#include "QuadForm.h"
#include "SquareMatrixInt.h"

// quadratic forms over the integers

template<typename R, size_t n>
class QuadFormInt
{
public:
  QuadFormInt()
    : _is_reduced(false),
      _num_aut(0),
      _num_aut_init(false)
  {}

  // We adhere to magma convention - giving the rows
  // up to the diagonal
  QuadFormInt(const typename R_QuadForm<R,n>::SymVec& coeffs)
    : _is_reduced(false), _num_aut(0), _num_aut_init(false) {}

  QuadFormInt(const SquareMatrix<Integer<R>, IntegerRing<R>, n> & B)
    : _is_reduced(false), _num_aut(0), _num_aut_init(false)  {}

  QuadFormInt(const SquareMatrixInt<R,n> &);
  
  // assignment
  QuadFormInt<R,n>& operator=(const QuadFormInt<R,n> &);

  // access
  inline const R & operator()(size_t i, size_t j) const {return _B(i,j); }
  inline R & operator()(size_t i, size_t j) {return _B(i,j); }
  
  R discriminant(void) const;
  
  inline R evaluate(const VectorInt<R,n>& vec) const
  { return VectorInt<R,n>::innerProduct(vec, (this->_B) * vec) / 2; }

  struct jordan_data {
    std::vector< MatrixRat<R> > matrices;
    std::vector< MatrixRat<R> > grams;
    std::vector<size_t> exponents;
  };

  typedef std::set<std::pair<Integer<R>, int> > QFInv;
  
  Integer<R> invariants(std::set< Integer<R> > &, size_t& ) const;
  
  Integer<R> invariants(QFInv &, size_t& ) const;
  
  jordan_data jordanDecomposition(const Integer<R> & p) const;

  template<typename S, typename T>
  std::shared_ptr< QuadFormFp<S,T,n> > mod(std::shared_ptr< Fp<S,T> >) const;

  VectorInt<R,n> orthogonalizeGram(void) const;

  static QuadFormZZ<R,n> reduce(const QuadFormZZ<R,n> & q,
				Isometry<R,n> & isom,
				bool calc_aut = false);

  std::unordered_set<Isometry<R,n>> properAutomorphisms() const;

  static std::vector< QuadFormZZ<R,5> > nippToForms(NippEntry entry);
  
  static std::vector<std::vector< QuadFormZZ<R,5> > >
  getQuinaryForms(const R & disc);

  static void greedy(SquareMatrixInt<R,n>& q, Isometry<R,n>& s, size_t dim = n);

  std::unordered_map< QuadFormZZ<R,n>, Isometry<R,n> > generateOrbit() const;

  inline bool isReduced(void) const { return this->_is_reduced; }

  size_t numAutomorphisms(void) const;

  inline void setNumAut(size_t num_aut)
  { this->_num_aut = num_aut; this->_num_aut_init = true; return;}

  inline void setReduced(void)
  { this->_is_reduced = true; return; }
  
protected:

  SquareMatrixInt<R,n> _B;
  // we save these for quick access
  
  bool _is_reduced;
  size_t _num_aut;
  bool _num_aut_init;
  
  // helper functions
  
  // reduce the form to a Minkowski reduced form
  // This is non-constant because we update the members
  // updates also the automorphism group of the lattice
  
  static size_t _iReduce(SquareMatrixInt<R,n> & qf,
			Isometry<R,n> & isom,
			std::unordered_set< Isometry<R,n> > & auts,
			bool calc_aut = true);

  static bool _permutationReduction(SquareMatrixInt<R,n> & qf,
				   Isometry<R,n> & isom,
				   std::unordered_set< Isometry<R,n> > & auts,
				   bool calc_aut = true);
  
  static bool _signNormalization(SquareMatrixInt<R,n> & qf,
				Isometry<R,n> & isom,
				std::unordered_set< Isometry<R,n> > & auts,
				bool calc_aut = true);
  
  static bool _normEchelon(SquareMatrixInt<R,n> & qf, Isometry<R,n> & isom);
  
  static bool _neighborReduction(SquareMatrixInt<R,n> & qf,
				Isometry<R,n> & isom,
				std::unordered_set< Isometry<R,n> > & auts,
				bool calc_aut = true);
  
  static size_t _generateAuts(std::unordered_set< Isometry<R,n> > & auts);

  static VectorInt<R,n-1> _voronoiBounds(size_t dim = n);
  
  // static helper functions

  static std::vector< std::vector<size_t> > _allPerms(size_t m);
  
  static int _hasse(const VectorInt<R,n>& , const Integer<R> & );

  // update in-place q and iso according to the closest vector
  // to the space spanned by the n-1 first ones
  static void _closestLatticeVector(SquareMatrixInt<R,n> &q,
				   Isometry<R,n> & iso,
				   size_t dim = n);

  static bool _signNormalizationSlow(SquareMatrixInt<R,n> & qf,
				    Isometry<R,n> & isom,
				    std::unordered_set< Isometry<R,n> > & auts);

  static bool _signNormalizationFast(SquareMatrixInt<R,n> & qf,
				    Isometry<R,n> & isom);

  static std::vector<uint8_t> _bitTranspose(const std::vector< uint8_t > & mat);
  static uint8_t _bitEchelonForm(std::vector< uint8_t > & mat,
				std::vector< uint8_t > & trans);
  
  static std::vector<uint8_t> _kernel(const std::vector< uint8_t > & mat);

  std::unordered_map< QuadFormZZ<R,n>, Isometry<R,n> >
  _permutationOrbit(void) const;
  
  std::unordered_map<  QuadFormZZ<R,n>, Isometry<R,n> >
  _signOrbit(void) const;
};

// we need this intermediate class for the partial specialization

template<typename R, size_t n>
class QuadFormZZ : public QuadFormInt<R,n>
{
public:
  QuadFormZZ() : QuadFormInt<R,n>() {}

  // a more general constructor
  // We adhere to magma convention - giving the rows
  // up to the diagonal
  QuadFormZZ(const typename QuadFormInt<R,n>::SymVec& coeffs)
    : QuadFormInt<R,n>(coeffs) {}

  QuadFormZZ(const SquareMatrixInt<R,n> & B)
    : QuadFormInt<R,n>(B) {}

  using QuadFormInt<R,n>::operator==;
    
  using QuadFormInt<R,n>::discriminant;
  
  W64 hashValue(void) const;

  using QuadFormInt<R,n>::evaluate;
  using QuadFormInt<R,n>::reduce;
  using QuadFormInt<R,n>::generateOrbit;

};

namespace std
{
  template<class R, size_t n>
  struct hash<QuadFormZZ<R,n>>
  {
    Z64 operator()(const QuadFormZZ<R,n>& q) const
    {
      return q.hashValue();
    }
  };
}

// Here we find that we must instantiate the following classes due to
// partial template specilization of hashValue.
  
template<size_t n>
class QuadFormZZ<Z,n> : public QuadFormInt<Z,n>
{
public:
  QuadFormZZ() : QuadFormInt<Z,n>() {}

  // a more general constructor
  // We adhere to magma convention - giving the rows
  // up to the diagonal
  QuadFormZZ(const typename R_QuadForm<Z,n>::SymVec& coeffs)
    : QuadFormInt<Z,n>(coeffs) {}

  QuadFormZZ(const SquareMatrixInt<Z,n> & B)
    : QuadFormInt<Z,n>(B) {}

  using QuadFormInt<Z,n>::operator==;
    
  using QuadFormInt<Z,n>::discriminant;
  
  W64 hashValue(void) const;

  using QuadFormInt<Z,n>::evaluate;
  using QuadFormInt<Z,n>::generateOrbit;
  using QuadFormInt<Z,n>::reduce;
  
  static Z_QuadForm<3> getQuadForm(const std::vector<Z_PrimeSymbol>& input);

};


template<size_t n>
class QuadFormZZ<Z64,n> : public QuadFormInt<Z64,n>
{
public:
  QuadFormZZ() : QuadFormInt<Z64,n>() {}

  // a more general constructor
  // We adhere to magma convention - giving the rows
  // up to the diagonal
  QuadFormZZ(const typename R_QuadForm<Z64,n>::SymVec& coeffs)
    : QuadFormInt<Z64,n>(coeffs) {}

  QuadFormZZ(const SquareMatrixInt<Z64, n> & B)
    : QuadFormInt<Z64,n>(B) {}

  using QuadFormInt<Z64,n>::operator==;
    
  using QuadFormInt<Z64,n>::discriminant;
  
  W64 hashValue(void) const;

  using QuadFormInt<Z64,n>::evaluate;
  using QuadFormInt<Z64,n>::reduce;
  
  static Z64_QuadForm<3> getQuadForm(const std::vector<Z64_PrimeSymbol>& input);

  using QuadFormInt<Z64,n>::generateOrbit;

};

template<size_t n>
class QuadFormZZ<Z128,n> : public QuadFormInt<Z128,n>
{
public:
  QuadFormZZ() : QuadFormInt<Z128,n>() {}

  // a more general constructor
  // We adhere to magma convention - giving the rows
  // up to the diagonal
  QuadFormZZ(const typename R_QuadForm<Z128,n>::SymVec& coeffs)
    : QuadFormInt<Z128,n>(coeffs) {}

  QuadFormZZ(const SquareMatrixInt<Z128,n> & B)
    : QuadFormInt<Z128,n>(B) {}

  using QuadFormInt<Z128,n>::operator==;
    
  using QuadFormInt<Z128,n>::discriminant;
  
  W64 hashValue(void) const;

  using QuadFormInt<Z128,n>::evaluate;
  using QuadFormInt<Z128,n>::reduce;
  
  static Z128_QuadForm<3> getQuadForm(const std::vector<Z128_PrimeSymbol>& input);

  using QuadFormInt<Z128,n>::generateOrbit;

};

#include "QuadFormInt.inl"

#endif // __QUAD_FORM_INT_H_
