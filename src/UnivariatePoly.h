#Ifndef __UNIVARIATE_POLY_H_
#define __UNIVARIATE_POLY_H_

#include "birch.h"

// !! TODO - make UnivariatePoly inherit from RingElement
// This will make some of the methods here redundant, and avoid code duplication.

template<class R, class Parent>
class UnivariatePoly
{
  // We would want R to be a Euclidean domain, to be able to perform things such as gcd and reduction
  static_assert(std::is_base_of<EuclideanDomainElement<R,Parent>,R>::value);
  static_assert(std::is_base_of<EuclideanDomain<Parent,R>,Parent>::value);
public:
  // create the zero polynomial
  UnivariatePoly(std::shared_ptr<const Parent> base_ring) : _base(base_ring) {}
  // create the constant polynomial
  UnivariatePoly(const R &);

  // create polynomial from coefficients
  UnivariatePoly(const std::vector<R> &);

  // create the polynomial x^i
  static UnivariatePoly<R,Parent> x(std::shared_ptr<const Parent> base_ring, size_t i = 1);
  
  // access
  // get methods
  R constCoefficient(void) const {return this->coefficient(0); }
  
  // coefficient of x^i
  R coefficient(size_t i) const;

  // leading coefficient
  R lead(void) const { return this->coefficient(this->degree()); }

  const std::vector<R> & coefficients(void) const
  {return this->_coeffs; }

  // if poly == 0, returns -1
  int degree(void) const {return this->_coeffs.size()-1; }

  inline std::shared_ptr<const Parent> baseRing(void) const
  { return this->_base; }
  
  R content(void) const;

  // conversion, assignment operator
  template<class T, class TParent>
  UnivariatePoly<R,Parent> & operator=(const UnivariatePoly<T,TParent> & );
  UnivariatePoly<R,Parent> & operator=(const R & );
  
  // arithmetic
  UnivariatePoly<R,Parent> operator-() const;
  UnivariatePoly<R,Parent> operator+(const UnivariatePoly<R,Parent> & ) const;
  UnivariatePoly<R,Parent> operator-(const UnivariatePoly<R,Parent> & ) const;
  UnivariatePoly<R,Parent> operator*(const UnivariatePoly<R,Parent> & ) const;
  UnivariatePoly<R,Parent> operator/(const UnivariatePoly<R,Parent> & ) const;
  UnivariatePoly<R,Parent> operator%(const UnivariatePoly<R,Parent> & ) const;
  UnivariatePoly<R,Parent> operator*(const R & ) const;
  UnivariatePoly<R,Parent> operator/(const R & ) const;
  UnivariatePoly<R,Parent> operator%(const R & ) const;

  UnivariatePoly<R,Parent> & operator+=(const UnivariatePoly<R,Parent> & );
  UnivariatePoly<R,Parent> & operator-=(const UnivariatePoly<R,Parent> & );
  UnivariatePoly<R,Parent> & operator*=(const UnivariatePoly<R,Parent> & );
  UnivariatePoly<R,Parent> & operator/=(const UnivariatePoly<R,Parent> & );
  UnivariatePoly<R,Parent> & operator%=(const UnivariatePoly<R,Parent> & );
  UnivariatePoly<R,Parent>& operator*=(const R & );
  UnivariatePoly<R,Parent>& operator/=(const R & );
  UnivariatePoly<R,Parent>& operator%=(const R & );
  
  UnivariatePoly<R,Parent> evaluate(const UnivariatePoly<R,Parent> &) const;
  R evaluate(const R &) const;
  
  template<class S, class SParent>
  Matrix<S,SParent> evaluate(const Matrix<S,SParent> &) const;

  // booleans
  bool isZero() const {return this->_coeffs.empty();}
  bool operator==(const UnivariatePoly<R,Parent> & ) const;
  bool operator!=(const UnivariatePoly<R,Parent> & ) const;
  bool operator==(const R & ) const;
  bool operator!=(const R & ) const;
  
  // algorithms
  UnivariatePoly<R,Parent> derivative(void) const;

  static void divRem(const UnivariatePoly<R,Parent> & f,
		     const UnivariatePoly<R,Parent> & g,
		     UnivariatePoly<R,Parent> & q,
		     UnivariatePoly<R,Parent> & r);
  
  static UnivariatePoly<R,Parent> gcd(const UnivariatePoly<R,Parent> &,
				      const UnivariatePoly<R,Parent> &);

  static UnivariatePoly<R,Parent> xgcd(const UnivariatePoly<R,Parent> & f,
				       const UnivariatePoly<R,Parent> & g,
				       UnivariatePoly<R,Parent> & s,
				       UnivariatePoly<R,Parent> & t);
  
protected:
  std::vector<R> _coeffs;
  std::shared_ptr< const Parent > _base;

  void _eliminateDeg(void);
 
};

template<class R, class Parent>
UnivariatePoly<R,Parent> operator*(const R & a,
				   const UnivariatePoly<R,Parent>  & poly)
{ return poly*a; }

template<class R, class Parent>
std::ostream& operator<<(std::ostream&, const UnivariatePoly<R, Parent> &);

namespace std
{
  template<typename R>
  struct hash< UnivariatePoly<R,Parent> >
  {
    Z64 operator()(const UnivariatePoly<R,Parent>& p) const
    {
      return p.hashValue();
    }
  };
}

#include "UnivariatePoly.inl"

#endif // __UNIVARIATE_POLY_H
