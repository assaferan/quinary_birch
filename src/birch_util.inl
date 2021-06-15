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
