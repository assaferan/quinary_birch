#ifndef __POLYNOMIAL_H_
#define __POLYNOMIAL_H_

#include <map>
#include <set>

#include "birch.h"
#include "PolynomialRing.h"
#include "RingElement.h"

// !! TODO - this works for arbitrary euclidean domains,
// can write the general version and specialize to Fp

template<typename R, typename S>
class PolynomialFp : public virtual RingElement< PolynomialFp<R,S>, PolynomialRingFp<R,S> >
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
  PolynomialFp<R,S> & operator=(const PolynomialFp<R,S> & ) override;
  PolynomialFp<R,S> & operator=(const FpElement<R,S> & );
  
  // arithmetic
  PolynomialFp<R,S> operator-() const override;
  PolynomialFp<R,S> operator+(const PolynomialFp<R,S> & ) const override;
  PolynomialFp<R,S> operator-(const PolynomialFp<R,S> & ) const override;
  PolynomialFp<R,S> operator*(const PolynomialFp<R,S> & ) const override;
  PolynomialFp<R,S> operator*(const FpElement<R,S> & ) const;

  PolynomialFp<R,S> & operator+=(const PolynomialFp<R,S> & ) override;
  PolynomialFp<R,S> & operator-=(const PolynomialFp<R,S> & ) override;
  PolynomialFp<R,S> & operator*=(const PolynomialFp<R,S> & ) override;
  PolynomialFp<R,S> & operator*=(const FpElement<R,S> & );

  FpElement<R,S> evaluate(const std::vector< FpElement<R,S> > & ) const;
  PolynomialFp<R,S> evaluate(const std::vector<PolynomialFp<R,S> > & vec) const;

  // booleans
  bool isZero(void) const override;
  bool isOne(void) const override;
  bool operator==(const PolynomialFp<R,S> & ) const override;
  bool operator!=(const PolynomialFp<R,S> & ) const override;
  bool operator==(const FpElement<R,S> & ) const;
  bool operator!=(const FpElement<R,S> & ) const;

  inline PolynomialFp<R,S> * getPtr(void) override {return this;}
  inline const PolynomialFp<R,S> * getPtr(void) const override {return this;}

  inline PolynomialFp<R,S> & makeZero(void) override
  { this->_mons.clear(); return (*this);}

  inline PolynomialFp<R,S> & makeOne(void) override
  { std::multiset<size_t> empty_set; this->_mons[empty_set] = this->_GF->one(); return (*this); }

  inline static PolynomialFp<R,S> zero(std::shared_ptr< const Fp<R,S> > GF)
  { PolynomialFp<R,S> z(GF); return z.makeZero(); }

  inline static PolynomialFp<R,S> one(std::shared_ptr< const Fp<R,S> > GF)
  { PolynomialFp<R,S> z(GF); return z.makeOne(); }

  void print(std::ostream&) const override;

  inline std::shared_ptr<const PolynomialRingFp<R,S> > parent(void) const override
  { return std::make_shared< const PolynomialRingFp<R,S> >(this->_GF);}
    
protected:
  std::shared_ptr<const Fp<R,S>> _GF;
  
  std::map< std::multiset<size_t>, FpElement<R,S> > _mons;
  
};

template<typename R, typename S>
inline PolynomialFp<R,S> operator*(const FpElement<R,S> & a,
				   const PolynomialFp<R,S>  & poly)
{ return poly*a; }

#include "Polynomial.inl"

#endif // __POLYNOMIAL_H_
