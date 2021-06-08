#ifndef __POLYNOMIAL_H_
#define __POLYNOMIAL_H_

#include "Ring.h"

class PolynomialRing : public Ring
{
public:
  PolynomialRing(const Ring & R) : base_(R) {}

  Polynomial zero(void) const;
  Polynomial one(void) const;

protected:
  const Ring & base_;
};

template<typename R>
class Polynomial : public RingElement
{
  static_assert(std::is_base_of<RingElement, R>::value, "R must inherit from RingElement");
public:
  
  // create the zero polynomial
  Polynomial(const Ring & R) : RingElement(PolynomialRing(R)) {}
  // create the constant polynomial
  Polynomial(const R &);

  // create polynomial from coefficients
  Polynomial(const std::vector<R> &);

  // create the polynomial x^i
  static Polynomial<R> x(size_t i = 1);
  
  // access
  // get methods
  R const_coefficient() const {return this->coefficient(0); }
  
  // coefficient of x^i
  R coefficient(size_t i) const;

  // leading coefficient
  R lead() const { return this->coefficient(this->degree()); }

  const std::vector<R> & coefficients() const
  {return this->coeffs; }

  // if poly == 0, returns -1
  int degree() const {return this->coeffs.size()-1; }

  R content() const;

  // conversion, assignment operator
  template<typename T>
  Polynomial<R> & operator=(const Polynomial<T> & );
  Polynomial<R> & operator=(const R & );
  
  // arithmetic
  Polynomial<R> operator-() const;
  Polynomial<R> operator+(const Polynomial<R> & ) const;
  Polynomial<R> operator-(const Polynomial<R> & ) const;
  Polynomial<R> operator*(const Polynomial<R> & ) const;
  Polynomial<R> operator/(const Polynomial<R> & ) const;
  Polynomial<R> operator%(const Polynomial<R> & ) const;
  Polynomial<R> operator*(const R & ) const;
  Polynomial<R> operator/(const R & ) const;
  Polynomial<R> operator%(const R & ) const;

  Polynomial<R> & operator+=(const Polynomial<R> & );
  Polynomial<R> & operator-=(const Polynomial<R> & );
  Polynomial<R> & operator*=(const Polynomial<R> & );
  Polynomial<R> & operator/=(const Polynomial<R> & );
  Polynomial<R> & operator%=(const Polynomial<R> & );
  Polynomial<R>& operator*=(const R & );
  Polynomial<R>& operator/=(const R & );
  Polynomial<R>& operator%=(const R & );

  template<typename S, typename T>
  PolynomialFp<S, T> mod(std::shared_ptr< const Fp<S, T> >) const;
  
  Polynomial<R> evaluate(const Polynomial<R> &) const;
  R evaluate(const R &) const;
  template<typename S>
  Matrix<S> evaluate(const Matrix<S> &) const;

  // booleans
  bool is_zero() const {return this->coeffs.empty();}
  bool operator==(const Polynomial<R> & ) const;
  bool operator!=(const Polynomial<R> & ) const;
  bool operator==(const R & ) const;
  bool operator!=(const R & ) const;

  // hash
  W64 hash_value(void) const;
  
  // algorithms
  Polynomial<R> derivative() const;

  static void div_rem(const Polynomial<R> & f,
		      const Polynomial<R> & g,
		      Polynomial<R> & q,
		      Polynomial<R> & r);
  
  static Polynomial<R> gcd(const Polynomial<R> &,
			       const Polynomial<R> &);

  static Polynomial<R> xgcd(const Polynomial<R> & f,
				const Polynomial<R> & g,
				Polynomial<R> & s,
				Polynomial<R> & t);
  
  std::unordered_map< Polynomial<R>, size_t > factor() const;
  
protected:
  std::vector<R> coeffs;

  void eliminate_deg();
  
  // these helper methods are needed for factorization
  
  template<typename S, typename T>
  void hensel_step(std::vector<Polynomial<R> > & u,
		   std::vector<Polynomial<R> > & v,
		   std::shared_ptr< const Fp<S,T> > GF,
		   size_t i) const;

  template<typename S, typename T>
  std::vector< Polynomial<R> >
  hensel_lift(const std::vector<PolynomialFp<S, T> > & g,
	      size_t a) const;

  std::vector< Polynomial<R> >
  trial_factor(const std::vector<Polynomial< R > > & u,
	       const R & N) const;

  void
  find_trial_factor(const std::vector< Polynomial<R> > & u,
		    const R & N,
		    size_t & j,
		    std::set<size_t> & C,
		    size_t & s,
		    std::vector< Polynomial<R> > & gs);

  R landau_mignotte() const;

  std::vector< Polynomial<R> > squarefree_factor() const;

  static std::set< std::set<size_t> >
  subsets(const std::set<size_t> & S, size_t k);
};

template<typename R>
Polynomial<R> operator*(const R & a,
			    const Polynomial<R>  & poly)
{ return poly*a; }

template<typename R>
std::ostream& operator<<(std::ostream&, const Polynomial<R> &);

namespace std
{
  template<typename R>
  struct hash< Polynomial<R> >
  {
    Z64 operator()(const Polynomial<R>& p) const
    {
      return p.hash_value();
    }
  };
}

template<typename R, typename S>
class PolynomialFp : public Polynomial< FpElement<R,S> >
{
public:
  PolynomialFp(std::shared_ptr< const Fp<R,S>> GF)
  { this->GF_ = GF; }

  // create the constant polynomial
  PolynomialFp(const FpElement<R, S> & a);
  
  // create polynomial from coefficients
  PolynomialFp(const std::vector< FpElement<R,S> > & v)
    : Polynomial< FpElement<R,S> >(v)
  {this->GF_ = v[0].field(); }

  PolynomialFp(const Polynomial< FpElement<R, S> > & other);
  
  // create the polynomial x^i
  static PolynomialFp<R,S> x(std::shared_ptr< const Fp<R,S>> GF,
				 size_t i = 1);
  
  // access
  const std::shared_ptr< const Fp<R,S> > & field() const
  {return this->GF_;}

  FpElement<R,S> coefficient(size_t i) const;
  FpElement<R,S> content() const;

  // arithmetic
  PolynomialFp<R, S> operator*(const PolynomialFp<R,S> & ) const;
  PolynomialFp<R, S> operator/(const PolynomialFp<R,S> & ) const;
  PolynomialFp<R, S> operator%(const PolynomialFp<R,S> & ) const;

  PolynomialFp<R,S> & operator*=(const PolynomialFp<R,S> & );
  PolynomialFp<R,S> & operator/=(const PolynomialFp<R,S> & );
  PolynomialFp<R,S> & operator%=(const PolynomialFp<R,S> & );

  PolynomialFp<R,S> derivative() const;
  
  std::vector< PolynomialFp<R,S> > sqf_factor() const;

  Polynomial<R> lift() const;
  
  PolynomialFp<R,S>
  pow_mod(size_t, const PolynomialFp<R,S> & ) const;

  // assignment and conversion
  using Polynomial< FpElement<R,S> >::operator=;
  PolynomialFp<R,S> & operator=(const Polynomial< FpElement<R,S> > &);
  
  
  // couldn't make it work with inheritance
  static PolynomialFp<R,S> gcd(const PolynomialFp<R,S> & f,
				   const PolynomialFp<R,S> & g);
				    
  static PolynomialFp<R,S> xgcd(const PolynomialFp<R,S> & f,
				    const PolynomialFp<R,S> & g,
				    PolynomialFp<R,S> & s,
				    PolynomialFp<R,S> & t);

  static void div_rem(const PolynomialFp<R,S> & f,
		      const PolynomialFp<R,S> & g,
		      PolynomialFp<R,S> & q,
		      PolynomialFp<R,S> & r);
  
  /*
  using Polynomial< FpElement<R,S> >::div_rem;
  using Polynomial< FpElement<R,S> >::gcd;
  using Polynomial< FpElement<R,S> >::xgcd;
  */
  
protected:
  std::shared_ptr< const Fp<R,S>> GF_;
  
  std::vector< PolynomialFp<R,S> >
  cz_eq_deg_partial_factor(size_t r) const;

  std::vector< PolynomialFp<R,S> > cz_eq_deg_factor(size_t r) const;

  std::vector< PolynomialFp<R,S> > cz_distinct_deg_factor() const;
  
};

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
  // PolynomialFp(std::shared_ptr<const Fp<R,S>> GF, size_t i);
  static PolynomialFp<R,S> x(std::shared_ptr<const Fp<R,S>> GF, size_t i);

  // create a polynomial from a bilinear form
  template<size_t n>
  PolynomialFp(const SquareMatrixFp<R, S, n> & );
  
  // access

  // get methods
  FpElement<R, S> const_coefficient() const;
  // coefficient of x_i
  FpElement<R, S> coefficient(size_t i) const;
  // coefficient of x_i x_j
  FpElement<R, S> coefficient(size_t i, size_t j) const;

  FpElement<R, S>
  coefficient(const std::multiset<size_t> & mon) const;

  const std::map< std::multiset<size_t>, FpElement<R,S> > & monomials() const
  {return mons; }

  PolynomialFp<R,S> quadratic_part() const;
  std::vector< FpElement<R,S> > linear_part(size_t rank) const;

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
  bool is_zero() const;
  bool operator==(const PolynomialFp<R, S> & ) const;
  bool operator!=(const PolynomialFp<R, S> & ) const;
  bool operator==(const FpElement<R, S> & ) const;
  bool operator!=(const FpElement<R, S> & ) const;

  
protected:
  std::shared_ptr<const Fp<R,S>> GF;
  
  std::map< std::multiset<size_t>, FpElement<R,S> > mons;
  
};

template<typename R, typename S>
PolynomialFp<R,S> operator*(const FpElement<R,S> & a,
			      const PolynomialFp<R,S>  & poly)
{ return poly*a; }

template<typename R, typename S>
std::ostream& operator<<(std::ostream&, const PolynomialFp<R,S>&);

#include "Polynomial.inl"

#endif // __POLYNOMIAL_H_
