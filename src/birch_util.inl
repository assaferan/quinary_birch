#include <memory>

namespace birch_util
{
  template<typename From, typename To>
  inline PrimeSymbol<To> convertPrimeSymbol(const PrimeSymbol<From>& symbol)
  {
    PrimeSymbol<To> temp;
    temp.p = convertInteger<From,To>(symbol.p);
    temp.power = symbol.power;
    temp.ramified = symbol.ramified;
    return temp;
  }

  template<typename From, typename To, size_t n>
  inline QuadFormZZ<To,n> convertQuadForm(const QuadFormZZ<From,n>& q)
  {
    SquareMatrixInt<To,n> mat;
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < n; j++)
	mat(i,j) = convertInteger<From,To>(q.bilinearForm()(i,j).num());
    QuadFormZZ<To,n> qq(mat);
    if (q.isReduced()) {
      qq.setReduced();
      qq.setNumAut(qq.numAutomorphisms());
    }
    
    return qq;
  }

  template<typename From, typename To, size_t n>
  inline Isometry<To,n> convertIsometry(const Isometry<From,n>& s)
  {
    SquareMatrixInt<To,n> mat;
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < n; j++)
	mat(i,j) = convertInteger<From,To>(s.a(i,j).num());
    To scale = convertInteger<From,To>(s.getScale().num());
    Isometry<To,n> ss(mat,scale);
    return ss;
  }

  template<typename From, typename To, size_t n>
  inline GenusRep<To,n> convertGenusRep(const GenusRep<From,n>& from)
  {
    GenusRep<To,n> to;
    to.q = birch_util::convertQuadForm<From,To,n>(from.q);
    to.s = birch_util::convertIsometry<From,To,n>(from.s);
    to.sinv = birch_util::convertIsometry<From,To,n>(from.sinv);
    to.parent = from.parent;
    to.p = convertInteger<From,To>(from.p);
    for (auto& pair : from.es)
      {
	to.es[convertInteger<From,To>(pair.first)] = pair.second;
      }
    return to;
  }


  template<typename R, size_t n>
  inline MatrixInt<R> convertMatrix(const SquareMatrixInt<R,n> & from)
  {
    std::shared_ptr<const IntegerRing<R> > ZZ = std::make_shared<const IntegerRing<R> >();
    MatrixInt<R> mat(ZZ, n, n);
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < n; j++)
	mat(i,j) = from(i,j);
    return mat;
  }
  
  template<typename R>
  inline R myPow(const std::map<R,int>& pairs)
  {
    R x = 1;
    for (const auto& pair : pairs)
      {
	for (int k=0; k<pair.second; k++)
	  {
	    x *= pair.first;
	  }
      }
    return x;
  }

  template<typename R>
  inline R gcd(const R & x, const R & y)
  {
    Z x_Z = convertInteger<R,Z>(x);
    Z y_Z = convertInteger<R,Z>(y);
    Z d_Z;
    mpz_gcd(d_Z.get_mpz_t(), x_Z.get_mpz_t(), y_Z.get_mpz_t());
    return convertInteger<Z,R>(d_Z);
  }

  template<typename R>
  inline R lcm(const R & x, const R & y)
  {
    Z x_Z = convertInteger<R,Z>(x);
    Z y_Z = convertInteger<R,Z>(y);
    Z d_Z;
    mpz_lcm(d_Z.get_mpz_t(), x_Z.get_mpz_t(), y_Z.get_mpz_t());
    return convertInteger<Z,R>(d_Z);
  }

  template<typename R>
  inline R xgcd(const R & x, const R & y, R & s, R& t)
  {
    Z x_Z = convertInteger<R,Z>(x);
    Z y_Z = convertInteger<R,Z>(y);
    Z d_Z, s_Z, t_Z;
    mpz_gcdext(d_Z.get_mpz_t(), s_Z.get_mpz_t(), t_Z.get_mpz_t(), x_Z.get_mpz_t(), y_Z.get_mpz_t());
    s = convertInteger<Z,R>(s_Z);
    t = convertInteger<Z,R>(t_Z);
    return convertInteger<Z,R>(d_Z);
  }

  extern int char_vals[256];
  
  inline int charVal(W64 x)
  {
    if (x <= 0xff) return char_vals[x];
    int value = 1;
    while (x)
      {
	value *= char_vals[x&0xff];
	x >>= 8;
      }
    return value;
  }
  
}

// or for std::vector
template<typename R>
inline std::ostream& operator<<(std::ostream& os, const std::vector<R>& v)
{
  if (v.size() >= 1) {
    for (size_t i = 0; i < v.size() - 1; i++)
      os << v[i] << ',';
  
    os << v[v.size()-1];
  }
  return os;
}


