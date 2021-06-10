#ifndef __QUAD_FORM_INT_H_
#define __QUAD_FORM_INT_H_

#include <memory>
#include <set>
#include <unordered_map>
#include <vector>

#include "birch.h"
#include "ParseNipp.h"
#include "QuadForm.h"

// quadratic forms over the integers

template<typename R, size_t n>
class QuadFormInt : public R_QuadForm<R,n>
{
public:
  QuadFormInt() : R_QuadForm<R,n>(), _is_reduced(false), _num_aut(0), _num_aut_init(false)
  {}

  // We adhere to magma convention - giving the rows
  // up to the diagonal
  QuadFormInt(const typename R_QuadForm<R,n>::SymVec& coeffs)
    : R_QuadForm<R,n>(coeffs), _is_reduced(false), _num_aut(0), _num_aut_init(false) {}

  QuadFormInt(const SquareMatrixInt<R,n> & B)
    : R_QuadForm<R,n>(B), _is_reduced(false), _num_aut(0), _num_aut_init(false)  {}

  // assignment
  QuadFormInt<R,n>& operator=(const QuadFormInt<R,n> &);
    
  using R_QuadForm<R,n>::discriminant;
  using R_QuadForm<R,n>::evaluate;
  using R_QuadForm<R,n>::reduce;
  using R_QuadForm<R,n>::generateOrbit;

  struct jordan_data {
    std::vector< SquareMatrixRat<R,n> > matrices;
    std::vector< SquareMatrixRat<R,n> > grams;
    std::vector<size_t> exponents;
  };

  typedef std::set<std::pair<Integer<R>, int> > QFInv;
  
  Integer<R> invariants(std::set< Integer<R> > &, size_t& ) const;
  
  Integer<R> invariants(QFInv &, size_t& ) const;
  
  jordan_data jordanDecomposition(const Integer<R> & p) const;

  template<typename S, typename T>
  std::shared_ptr< QuadFormFp<S,T,n> > mod(std::shared_ptr< Fp<S,T> >) const;

  VectorInt<R,n> orthogonalizeGram() const;

  static QuadFormInt<R,n> reduce(const QuadFormInt<R,n> & q,
				 Isometry<R,n> & isom,
				 bool calc_aut = false);

  std::set<Isometry<R,n>> properAutomorphisms() const;

  static std::vector< QuadFormInt<R,5> > nippToForms(NippEntry entry);
  
  static std::vector<std::vector< QuadFormInt<R,5> > >
  getQuinaryForms(const Integer<R> & disc);

  static void greedy(SquareMatrixInt<R,n>& q, Isometry<R,n>& s, size_t dim = n);

  std::unordered_map< QuadFormInt<R,n>, Isometry<R,n> > generateOrbit() const;

  inline bool isReduced() const { return this->_is_reduced; }

  size_t numAutomorphisms() const;

  inline void setNumAut(size_t num_aut)
  { this->_num_aut = num_aut; this->_num_aut_init = true; return;}

  inline void setReduced()
  { this->_is_reduced = true; return; }
  
protected:

  // we save these for quick access
  
  bool _is_reduced;
  size_t _num_aut;
  bool _num_aut_init;
  
  // helper functions
  
  // reduce the form to a Minkowski reduced form
  // This is non-constant because we update the members
  // updates also the automorphism group of the lattice
  
  static size_t iReduce(SquareMatrixInt<R,n> & qf,
			Isometry<R,n> & isom,
			std::set< Isometry<R,n> > & auts,
			bool calc_aut = true);

  static bool permutationReduction(SquareMatrixInt<R,n> & qf,
				   Isometry<R,n> & isom,
				   std::set< Isometry<R,n> > & auts,
				   bool calc_aut = true);
  
  static bool signNormalization(SquareMatrixInt<R,n> & qf,
				Isometry<R,n> & isom,
				std::set< Isometry<R,n> > & auts,
				bool calc_aut = true);
  
  static bool normEchelon(SquareMatrixInt<R,n> & qf, Isometry<R,n> & isom);
  
  static bool neighborReduction(SquareMatrixInt<R,n> & qf,
				Isometry<R,n> & isom,
				std::set< Isometry<R,n> > & auts,
				bool calc_aut = true);
  
  static size_t generateAuts(std::set< Isometry<R,n> > & auts);

  static VectorInt<R,n-1> voronoiBounds(size_t dim = n);
  
  // static helper functions

  static std::vector< std::vector<size_t> > allPerms(size_t m);
  
  static int hasse(const VectorInt<R,n>& , const R & );

  // update in-place q and iso according to the closest vector
  // to the space spanned by the n-1 first ones
  static void closestLatticeVector(SquareMatrixInt<R,n> &q,
				   Isometry<R,n> & iso,
				   size_t dim = n);

  static bool signNormalizationSlow(SquareMatrixInt<R,n> & qf,
				    Isometry<R,n> & isom,
				    std::set< Isometry<R,n> > & auts);

  static bool signNormalizationFast(SquareMatrixInt<R,n> & qf,
				    Isometry<R,n> & isom);

  static std::vector<uint8_t> bitTranspose(const std::vector< uint8_t > & mat);
  static uint8_t bitEchelonForm(std::vector< uint8_t > & mat,
				std::vector< uint8_t > & trans);
  
  static std::vector<uint8_t> kernel(const std::vector< uint8_t > & mat);

  std::unordered_map< QuadFormInt<R,n>, Isometry<R,n> >
  permutationOrbit() const;
  
  std::unordered_map<  QuadFormInt<R,n>, Isometry<R,n> >
  signOrbit() const;
};

// Here we find that we must instantiate the following classes due to
// partial template specilization of hashValue.
  
template<size_t n>
class Z_QuadForm : public QuadFormInt<Z,n>
{
public:
  Z_QuadForm() : QuadFormInt<Z,n>() {}

  // a more general constructor
  // We adhere to magma convention - giving the rows
  // up to the diagonal
  Z_QuadForm(const typename QuadFormInt<Z,n>::SymVec& coeffs)
    : QuadFormInt<Z,n>(coeffs) {}

  Z_QuadForm(const SquareMatrixInt<Z,n> & B)
    : QuadFormInt<Z,n>(B) {}

  using QuadFormInt<Z,n>::operator==;
    
  using QuadFormInt<Z,n>::discriminant;
  
  W64 hashValue(void) const;

  using QuadFormInt<Z,n>::evaluate;
  using QuadFormInt<Z,n>::reduce;
  
  static Z_QuadForm<3> getQuadForm(const std::vector<Z_PrimeSymbol>& input);

  using QuadFormInt<Z,n>::generateOrbit;

};

template<size_t n>
class Z64_QuadForm : public QuadFormInt<Z64,n>
{
public:
  Z64_QuadForm() : QuadFormInt<Z64,n>() {}

  // a more general constructor
  // We adhere to magma convention - giving the rows
  // up to the diagonal
  Z64_QuadForm(const typename QuadFormInt<Z64,n>::SymVec& coeffs)
    : QuadFormInt<Z64,n>(coeffs) {}

  Z64_QuadForm(const SquareMatrixInt<Z64,n> & B)
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
class Z128_QuadForm : public QuadFormInt<Z128,n>
{
public:
  Z128_QuadForm() : QuadFormInt<Z128,n>() {}

  // a more general constructor
  // We adhere to magma convention - giving the rows
  // up to the diagonal
  Z128_QuadForm(const typename QuadFormInt<Z128,n>::SymVec& coeffs)
    : QuadFormInt<Z128,n>(coeffs) {}

  Z128_QuadForm(const SquareMatrixInt<Z128,n> & B)
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
