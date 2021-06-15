#ifndef __POLYNOMIAL_H_
#define __POLYNOMIAL_H_

#include <map>
#include <set>

#include "birch.h"

// !! TODO - this works for arbitrary euclidean domains,
// can write the general version and specialize to Fp

template<typename R, typename S>
class PolynomialFp
{
public:

  // create the zero polynomial
  PolynomialFp(std::shared_ptr<const Fp<R,S>> GF);
  // create the constant polynomial
  PolynomialFp(const FpElement<R,S> & a);
  PolynomialFp(std::shared_ptr<const Fp<R,S>> GF, const R & a);
  
  // create the polynomial x_i
  static PolynomialFp<R,S> x(std::shared_ptr<const Fp<R,S>> GF, size_t i);

  // create a polynomial from a bilinear form
  template<size_t n>
  PolynomialFp(const SquareMatrixFp<R,S,n> & );
  
  // access

  // get methods
  FpElement<R,S> constCoefficient(void) const;
  // coefficient of x_i
  FpElement<R,S> coefficient(size_t i) const;
  // coefficient of x_i x_j
  FpElement<R,S> coefficient(size_t i, size_t j) const;

  FpElement<R, S>
  coefficient(const std::multiset<size_t> & mon) const;

  inline const std::map< std::multiset<size_t>, FpElement<R,S> > & monomials(void) const
  {return _mons; }

  PolynomialFp<R,S> quadraticPart(void) const;
  std::vector< FpElement<R,S> > linearPart(size_t rank) const;

  int degree(size_t i) const;

  // conversion, assignment operator
  PolynomialFp<R,S> & operator=(const PolynomialFp<R,S> & );
  PolynomialFp<R,S> & operator=(const FpElement<R,S> & );
  
  // arithmetic
  PolynomialFp<R,S> operator-() const;
  PolynomialFp<R,S> operator+(const PolynomialFp<R,S> & ) const;
  PolynomialFp<R,S> operator-(const PolynomialFp<R,S> & ) const;
  PolynomialFp<R,S> operator*(const PolynomialFp<R,S> & ) const;
  PolynomialFp<R,S> operator*(const FpElement<R,S> & ) const;

  PolynomialFp<R,S> & operator+=(const PolynomialFp<R,S> & );
  PolynomialFp<R,S> & operator-=(const PolynomialFp<R,S> & );
  PolynomialFp<R,S> & operator*=(const PolynomialFp<R,S> & );
  PolynomialFp<R,S> & operator*=(const FpElement<R,S> & );

  FpElement<R,S> evaluate(const std::vector< FpElement<R,S> > & ) const;
  PolynomialFp<R,S> evaluate(const std::vector<PolynomialFp<R,S> > & vec) const;

  // booleans
  bool isZero() const;
  bool operator==(const PolynomialFp<R,S> & ) const;
  bool operator!=(const PolynomialFp<R,S> & ) const;
  bool operator==(const FpElement<R,S> & ) const;
  bool operator!=(const FpElement<R,S> & ) const;

  
protected:
  std::shared_ptr<const Fp<R,S>> _GF;
  
  std::map< std::multiset<size_t>, FpElement<R,S> > _mons;
  
};

template<typename R, typename S>
inline PolynomialFp<R,S> operator*(const FpElement<R,S> & a,
				   const PolynomialFp<R,S>  & poly)
{ return poly*a; }

template<typename R, typename S>
std::ostream& operator<<(std::ostream&, const PolynomialFp<R,S>&);

#include "Polynomial.inl"

#endif // __POLYNOMIAL_H_
