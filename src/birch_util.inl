namespace birch_util
{
  template<typename From, typename To>
  inline PrimeSymbol<To> convert_PrimeSymbol(const PrimeSymbol<From>& symbol)
  {
    PrimeSymbol<To> temp;
    temp.p = convert_Integer<From,To>(symbol.p);
    temp.power = symbol.power;
    temp.ramified = symbol.ramified;
    return temp;
  }

  template<typename From, typename To, size_t n>
  inline QuadForm<To,n> convert_QuadForm(const QuadForm<From,n>& q)
  {
    SquareMatrix<To, n> mat;
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < n; j++)
	mat(i,j) = convert_Integer<From,To>(q.bilinear_form()(i,j));
    QuadForm<To,n> qq(mat);
    if (q.is_reduced()) {
      qq.set_reduced();
      qq.set_num_aut(qq.num_automorphisms());
    }
    
    return qq;
  }

  template<typename From, typename To, size_t n>
  inline Isometry<To,n> convert_Isometry(const Isometry<From,n>& s)
  {
    SquareMatrix<To, n> mat;
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < n; j++)
	mat(i,j) = convert_Integer<From,To>(s.a(i,j));
    To scale = convert_Integer<From, To>(s.get_scale());
    Isometry<To,n> ss(mat, scale);
    return ss;
  }

  template<typename From, typename To, size_t n>
  inline GenusRep<To,n> convert_GenusRep(const GenusRep<From,n>& from)
  {
    GenusRep<To,n> to;
    to.q = birch_util::convert_QuadForm<From,To,n>(from.q);
    to.s = birch_util::convert_Isometry<From,To,n>(from.s);
    to.sinv = birch_util::convert_Isometry<From,To,n>(from.sinv);
    to.parent = from.parent;
    to.p = convert_Integer<From,To>(from.p);
    for (auto& pair : from.es)
      {
	to.es[convert_Integer<From,To>(pair.first)] = pair.second;
      }
    return to;
  }

  template<typename R>
  inline R my_pow(const std::map<R,int>& pairs)
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

  inline int char_val(W64 x)
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
