#include <limits>
#include <map>
#include <random>

#include "birch_util.h"
#include "Fp.h"
#include "FpElement.h"
#include "ParseNipp.h"
#include "Matrix.h"
#include "SquareMatrix.h"

template<typename R, size_t n>
inline QuadFormInt<R,n>&
QuadFormInt<R,n>::operator=(const QuadFormInt<R,n> & other)
{
  R_QuadForm<R,n>::operator=(other);
  if (this != &other) {
    this->_is_reduced = other._is_reduced;
    this->_num_aut = other._num_aut;
    this->_num_aut_init = other._num_aut_init;
  }
  return *this;
}

template<typename R, size_t n>
template<typename S, typename T>
inline std::shared_ptr< QuadFormFp<S,T,n> >
QuadFormInt<R,n>::mod(std::shared_ptr< Fp<S,T> > GF) const
{
  SquareMatrixFp<S,T,n> q_mod(GF);
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      q_mod(i,j) = GF->mod(this->_B(i,j));
  R p = GF->prime();
  if (p == 2) {
    for (size_t i = 0; i < n; i++) {
      R value = this->_B(i,i) / 2;
      q_mod(i,i) = GF->mod(value);
    }
  }
  std::shared_ptr< QuadFormFp<S,T,n> > q =
    std::make_shared< QuadFormFp<S,T,n> >(q_mod);
  return q;
}

template<typename R, size_t n>
inline std::vector< QuadFormInt<R,5> >
QuadFormInt<R,n>::nippToForms(NippEntry entry)
{
  std::vector< QuadFormInt<R,5> > forms;
  size_t triangular[5];
  for (size_t j = 0; j < 5; j++)
    triangular[j] = j*(j-1)/2;
  typename QuadFormInt<R,5>::SymVec form;
  for (LatticeRecord lat : entry.lattices)
    {
      size_t form_idx = 0;
      for (size_t col = 0; col < 5; col++)
	{
	  for (size_t row = 0; row < col; row++)
	    {
	      form[form_idx++] = lat.form[5+triangular[col]+row]; 
	    }
	  form[form_idx++] = 2*lat.form[col];
	}
      forms.push_back(form);
    }
  return forms;
}

template<typename R, size_t n>
inline std::vector<std::vector< QuadFormInt<R,5> > >
QuadFormInt<R,n>::getQuinaryForms(const Integer<R> & disc)
{
  std::vector< std::vector< QuadFormInt<R,5> > > all_forms;

  std::vector<R> nipp_maxs = {0,256,270,300,322,345,400,440,480,500,513};
  size_t table_idx = 0;
  while (nipp_maxs[table_idx+1] < disc) table_idx++;
  std::ostringstream nipp_fname;
  nipp_fname << "lattice_db/nipp" << nipp_maxs[table_idx]+1 << "-";
  nipp_fname << nipp_maxs[table_idx+1] << ".txt";
  
#ifdef DEBUG_LEVEL_FULL
  std::cerr << "nipp_fname = " << nipp_fname.str() << std::endl;
#endif
  
  std::vector<NippEntry> nipps =
    ParseNipp::parseDisc(nipp_fname.str(),
			 birch_util::convert_Integer<R,Z>(disc));
  
  for (NippEntry nipp : nipps)
    {
#ifdef DEBUG_LEVEL_FULL
      std::cerr << "disc = " << nipp.disc << std::endl;
      std::cerr << "genus = " << nipp.genus << std::endl;
      std::cerr << "mass = " << nipp.mass[0] << "/" << nipp.mass[1] << std::endl;
      std::cerr << "Hasse symbols = ";
      for (short int symb : nipp.HasseSymb)
	std::cerr << symb << " ";
      std::cerr << std::endl;
      std::cerr << "lattices = " << std::endl;
      for (LatticeRecord lat : nipp.lattices)
	{
	  for (Z num : lat.form)
	    std::cerr << num << " ";
	  std::cerr << ";\t" << lat.numAut << std::endl; 
	}
#endif
      all_forms.push_back(QuadFormInt<R, 5>::nippToForms(nipp));
    }
  
  return all_forms;
}

// This is somewhat of a duplicate for cholesky,
// but this one keeps everything integral.
// Do we want to have a version for arbitrary euclidean domains?
// Maybe replace cholesky in SquareMatrix with this version.
template<typename R, size_t n>
inline VectorInt<R,n> QuadFormInt<R,n>::orthogonalizeGram() const
{
  std::shared_ptr<const IntegerRing<R> > ring = this->baseRing();
  VectorInt<R,n> D(ring);
  SquareMatrixInt<R,n> L(ring);
  Integer<R> prod_diag = ring->one();
  Integer<R> d = ring->zero();
  Integer<R> inner_sum = ring->zero();
  // This works but inefficiently - for some reason we get O(n^4) operations.
  // !! TODO - check it out later
  // Oh I see - we should do the L update in two passes...
  for (size_t i = 0; i < n; i++)
    {
      L(i,i) = prod_diag;
      d = prod_diag;
      for (size_t j = 0; j < i; j++)
	{
	  L(i,j).makeZero();
	  for (size_t k = j; k < i; k++)
	    {
	      inner_sum.makeZero();
	      for (size_t r = 0; r <= k; r++)
		inner_sum += L(k, r)*(this->_B(i,r))*L(k,j);
	      inner_sum *= -L(i, i) / D[k];
	      L(i,j) += inner_sum;
	    }
	  d = d.gcd(L(i, j));
	}
      for (size_t j = 0; j <= i; j++)
	L(i,j) /= d;
      D[i].makeZero();
      for (size_t j = 0; j <= i; j++)
	for (size_t k = 0; k <= i; k++)
	  D[i] += L(i, j)*(this->_B(j,k))*L(i, k);
      prod_diag = prod_diag.lcm(D[i]);
      for (size_t j = i+1; j < n; j++)
	L(i,j).makeZero();
    }

  // Recall that this is an even lattice, so all entries in D
  // are even, and we are more interested in their half values,
  // which corresponds to the quadratic form.
  for (size_t i = 0; i < n; i++)
    D[i] /= 2;
  
#ifdef DEBUG_LEVEL_FULL
  std::cerr << "L=" << std::endl << L << std::endl;
#endif
  return D;
}

template<typename R, size_t n>
inline int QuadFormInt<R,n>::hasse(const VectorInt<R,n> & D, const R & p)
{
  std::shared_ptr<const IntegerRing<R> > ring = D.baseRing();
  int hasse = 1;
  Integer<R> prod = ring->one();
  for (size_t i = 0; i < n; i++)
    prod *= D[i];
  for (size_t i = 0; i < n-1; i++)
    {
      prod /= D[i];
      hasse *= D[i].hilbertSymbol(prod, p);
    }
  return hasse;
}

// !! - TODO - merge the code duplication among the following two functions

template<typename R, size_t n>
inline Integer<R>
QuadFormInt<R,n>::invariants(std::set<Integer<R>> & F, size_t& I) const
{
  VectorInt<R,n> D = this->orthogonalizeGram();
  std::set<Integer<R> > P;
  F.clear();
  I = 0;
  
  P.insert(2);
  for (size_t i = 0; i < n; i++)
    {
      if (D[i] < 0) I++;
      typename Integer<R>::FactorData facs = D[i].factorization();
      for (std::pair<Integer<R>, size_t> fa : facs)
	  if (fa.second % 2 == 1)
	    P.insert(fa.first);
    }
  for (Integer<R> p : P)
     if (hasse(D,p) == -1) F.insert(p);

  Integer<R> prod = 1;
  for (size_t i = 0; i < n; i++)
    prod *= D[i];
  
  return prod;
}

template<typename R, size_t n>
inline Integer<R>
QuadFormInt<R,n>::invariants(typename QuadFormInt<R,n>::QFInv &F,
			     size_t& I) const
{
  VectorInt<R, n> D = this->orthogonalizeGram();
  std::set<Integer<R> > P;
  F.clear();
  I = 0;
  
  P.insert(2);
  for (size_t i = 0; i < n; i++)
    {
      if (D[i] < 0) I++;
      typename Integer<R>::FactorData facs = D[i].factorization();
      for (std::pair<Integer<R>, size_t> fa : facs)
	  if (fa.second % 2 == 1)
	    P.insert(fa.first);
    }
  for (Integer<R> p : P)
    F.insert(std::make_pair(p, hasse(D,p)));

  Integer<R> prod = 1;
  for (size_t i = 0; i < n; i++)
    prod *= D[i];
  
  return prod;
}

template<typename R, size_t n>
inline typename QuadFormInt<R,n>::jordan_data
QuadFormInt<R,n>::jordanDecomposition(const Integer<R> & p) const
{
  bool even = (p == 2);
  std::shared_ptr< const RationalField<R> > QQ = RationalField<R>::getInstance().getPtr();
  SquareMatrixRat<R,n> S = SquareMatrixRat<R,n>::identity(QQ);
  SquareMatrixRat<R,n> G(QQ);
  MatrixRat<R> F(QQ,n,n);
  
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      F(i,j) = this->_B(i,j);
  
  size_t k = 0;
  // virtually infinity
  size_t old_val = std::numeric_limits<size_t>::max();
  std::vector<size_t> blocks;
  jordan_data jordan;
  
  while (k < n)
    {
#ifdef DEBUG_LEVEL_FULL
      std::cerr << "k = " << k << std::endl;
#endif
      // G = SFS^t
      // !! TODO - can we write simply G = S*(this->_B)*S.transpose() ?
     for (size_t i = 0; i < n; i++)
       for (size_t j = 0; j < n; j++)
	 G(i, j) =
	   SquareMatrixInt<R,n>::innerProduct(this->_B, S, i, j);
#ifdef DEBUG_LEVEL_FULL
     std::cerr << "G = " << std::endl;
     G.prettyPrint(std::cerr);
#endif
     size_t ii = k;
     // infty
     size_t m = std::numeric_limits<size_t>::max();

     Integer<R> zero = 0;
     
     for (size_t i = k; i < n; i++)
       {
	 if (G(i,i) != zero) {
	   size_t val = G(i,i).valuation(p);
	   if (val < m)
	     {
	       m = val;
	       ii = i;
	     }
	 }
       }
     std::pair<size_t, size_t> i_pair = std::make_pair(ii, ii);
     for (size_t i = k; i < n; i++)
       for (size_t j = i+1; j < n; j++)
	 {
	   if (G(i,j) != zero) {
	     size_t tmp = G(i,j).valuation(p);
	     if (tmp < m)
	       {
		 m = tmp;
		 i_pair.first = i;
		 i_pair.second = j;
	       }
	   }
	 }
     
#ifdef DEBUG_LEVEL_FULL
     std::cerr << "i_pair = (" << i_pair.first << "," << i_pair.second << ")";
     std::cerr << std::endl << "m = " << m << std::endl;
#endif
     
     if (m != old_val)
       {
	 blocks.push_back(k);
	 old_val = m;
	 jordan.exponents.push_back(m);
       }
     
#ifdef DEBUG_LEVEL_FULL
     std::cerr << "blocks = " << blocks << std::endl;
     std::cerr << "jordan.exponents = " << jordan.exponents << std::endl;
#endif
     
     if ((even) && (i_pair.first != i_pair.second))
       {
	 S.swap_rows(i_pair.first, k);
	 S.swap_rows(i_pair.second, k+1);
	 
	 // T12 = S[k]*F*S[k+1]^t
	 Rational<R> T12 =
	   SquareMatrixRat<R,n>::inner_product(this->_B, S, k, k+1);

	 // multiply S[k] by p^val(T12,p)/T12
	 // Check whether we have to change to rational here
	 for (size_t i = 0; i < n; i++) {
	   Integer<R> val = (1 << T12.valuation(p));
	   S(k,i) *= val / T12;
	 }
	 Rational<R> T11 =
	   SquareMatrixInt<R,n>::innerProduct(this->_B, S, k, k);
	 Rational<R> T22 =
	   SquareMatrixInt<R,n>::innerProduct(this->_B, S, k+1, k+1);
	 T12 = SquareMatrixInt<R,n>::innerProduct(this->_B, S, k, k+1);
	 Rational<R> d = T11*T22-T12*T12;
	 for (size_t l = k+2; l < n; l++)
	   {
	     Rational<R> tl =
	       T12*SquareMatrixInt<R,n>::innerProduct(this->_B,S,k+1,l) -
	       T22*SquareMatrixInt<R,n>::innerProduct(this->_B,S,k,l);
	     Rational<R> ul =
	       T12*SquareMatrixInt<R,n>::innerProduct(this->_B,S,k,l) -
	       T11*SquareMatrixInt<R,n>::innerProduct(this->_B,S,k+1,l);
	     for (size_t i = 0; i < n; i++)
	       S(l,i) += (tl/d)*S(k,i) + (ul/d)*S(k+1,i);
	   }
	 k += 2;
       }
     else
       {
	 if (i_pair.first == i_pair.second) {
	   
#ifdef DEBUG_LEVEL_FULL
	   std::cerr << "swapping rows" << std::endl;
#endif
	   S.swapRows(i_pair.first, k);
	   
#ifdef DEBUG_LEVEL_FULL
	   std::cerr << "S = " << std::endl;
	   S.prettyPrint(std::cerr);
#endif
	 }
	 else
	   {
	     Integer<R> one = 1;
	     // std::cerr << "adding rows" << std::endl;
	     S.addRow(i_pair.first, i_pair.second, one);
	     
#ifdef DEBUG_LEVEL_FULL
	     std::cerr << "S = " << std::endl;
	     S.prettyPrint(std::cerr);
	     std::cerr << "swapping rows" << std::endl;
#endif
	     S.swapRows(i_pair.first, k);
#ifdef DEBUG_LEVEL_FULL
	     std::cerr << "S = " << std::endl;
	     S.prettyPrint(std::cerr);
#endif
	   }
	 Rational<R> nrm = SquareMatrixInt<R,n>::innerProduct(this->_B, S, k, k);

#ifdef DEBUG_LEVEL_FULL
	 std::cerr << "nrm = " << nrm << std::endl;
#endif
	 
	 Rational<R> X[n];
	 for (size_t i = 0; i < n; i++)
	   X[i] = SquareMatrixInt<R,n>::innerProduct(this->_B, S, k, i);
	 
#ifdef DEBUG_LEVEL_FULL
	 std::cerr << "X = " << X << std::endl;;
#endif
	 for (size_t l = k+1; l < n; l++)
	     for (size_t i = 0; i < n; i++)
	       S(l,i) -= X[l]/nrm * S(k, i);
	 
#ifdef DEBUG_LEVEL_FULL
         std::cerr << "S = " << std::endl;
	 S.prettyPrint(std::cerr);
#endif
	 k += 1;
       }
    }
  blocks.push_back(n);
  
#ifdef DEBUG_LEVEL_FULL
  std::cerr << "blocks = " << blocks << std::endl;
#endif
  
  for (size_t i = 0; i < blocks.size()-1; i++) {
    size_t nrows = blocks[i+1]-blocks[i];
    std::vector< Rational<R> > data(nrows*n);
    size_t idx = 0;
    for (size_t row = 0; row < nrows; row++)
      for (size_t col = 0; col < n; col++)
	data[idx++] = S(blocks[i]+row, col);
    MatrixRat<R> mat(data, nrows, n);
    jordan.matrices.push_back(mat);
  }
  
  for (MatrixRat<R> m  : jordan.matrices) {
#ifdef DEBUG_LEVEL_FULL
    std::cerr << "m = " << m << std::endl;
    std::cerr << "F = " << F << std::endl;
    std::cerr << "m^t = " << m.transpose() << std::endl;
    Rational<R> tmp_rat = m(0,0)*F(0,0);
    
    std::cerr << "tmp_rat = " << tmp_rat << std::endl;
    MatrixRat<R> tmp = m*F;
    
    std::cerr << "m*F = " << tmp << std::endl;
    MatrixRat<R> tmp2 = m.transpose();
    MatrixRat<R> tmp3 = tmp*tmp2;

    std::cerr << "m*F*m^t = " << tmp3 << std::endl;
#endif
    jordan.grams.push_back(m*F*m.transpose());
  }
#ifdef DEBUG_LEVEL_FULL
  std::cerr << "jordan.matrices = " << std::endl;
  for (size_t i = 0; i < jordan.matrices.size(); i++)
    std::cerr << std::endl << jordan.matrices[i] << std::endl;

  std::cerr << "jordan.grams = " << std::endl;
  for (size_t i = 0; i < jordan.grams.size(); i++)
    std::cerr << std::endl << jordan.grams[i] << std::endl;

  std::cerr << "jordan.exponents = ";
  for (size_t i = 0; i < jordan.exponents.size(); i++)
    std::cerr << jordan.exponents[i] << " ";
  std::cerr << std::endl;
#endif
  return jordan;
}

template<typename R, size_t n>
inline VectorInt<R,n-1> QuadFormInt<R,n>::voronoiBounds(size_t dim)
{
  // !! TODO - check what the real bounds are !!
  VectorInt<R,n-1> bounds;
  for (size_t i = 0; i < dim; i++)
    bounds[i] = 1;
  return bounds;
}

template<typename R, size_t n>
inline void QuadFormInt<R,n>::closestLatticeVector(SquareMatrixInt<R,n> &q,
					    Isometry<R,n> & iso,
					    size_t dim)
{
  Isometry<R,n> g, min_g;
  SquareMatrixInt<R,n> x_gram;
  SquareMatrixInt<R,n-1> H_int;
  VectorInt<R,n-1> v_int;
  std::shared_ptr< const IntegerRing<R> > ZZ = IntegerRing<R>::getInstance().getPtr();

#ifdef DEBUG_LEVEL_FULL
  std::cerr << "finding closest_lattice_vector with gram:" << std::endl;
  q.prettyPrint(std::cerr, dim);
#endif
  
  for (size_t i = 0; i < dim-1; i++) {
    for (size_t j = 0; j < dim-1; j++) {
      H_int(i,j) = q(i,j);
    }
    v_int[i] = q(i,dim-1);
  }

  H_int = H_int.adjugate(dim-1);
  
#ifdef DEBUG_LEVEL_FULL
  std::cerr << "H_int = " << std::endl;
  H_int.prettyPrint(std::cerr, dim-1);

  std::cerr << "v_int = " << std::endl;
  v_int.prettyPrint(std::cerr, dim-1);
#endif

  VectorInt<R,n-1> y_int = v_int*H_int.transpose();

#ifdef DEBUG_LEVEL_FULL
  std::cerr << "y_int = " << std::endl;
  y_int.prettyPrint(std::cerr, dim-1);
#endif
  
  VectorInt<R,n-1> voronoi = voronoiBounds(dim-1);
  VectorInt<R,n-1> x, x_min, x_max, x_num;
  VectorInt<R,n-1> x_closest;

  // This can be calculated more efficiently
  Integer<R> det = ZZ->zero();
  for (size_t i = 0; i < dim-1; i++)
    det += H_int(0,i)*q(i,0);
  det = (det > 0) ?  det : -det;
  for (size_t i = 0; i < dim-1; i++) {
    Integer<R> tmp =  y_int[i] - det*voronoi[i];
    x_min[i] = ((tmp >= 0) ? tmp+det-1 : tmp)/det;
  }
  for (size_t i = 0; i < dim-1; i++) {
    Integer<R> tmp =  y_int[i] + det*voronoi[i];
    x_max[i] = ((tmp >= 0) ? tmp : tmp-det+1)/det;
  }
  
  for (size_t i = 0; i < dim-1; i++)
    x_num[i] = x_max[i] - x_min[i] + 1;
  Integer<R> num_xs = 1;
  for (size_t i = 0; i < dim-1; i++)
    num_xs *= x_num[i];
  // This should be infinity
  Integer<R> min_dist = std::numeric_limits<R>::max();
  for (Integer<R> x_idx = 0; x_idx < num_xs; x_idx++) {
    Integer<R> tmp = x_idx;
    for (size_t i = 0; i < dim-1; i++) {
      size_t j = dim-2-i;
      x[j] = x_min[j] + (tmp % x_num[j]);
      tmp /= x_num[j];
    }
    for (size_t i = 0; i < dim-1; i++)
      g(i,dim-1) = -x[i];
    x_gram = g.transform(q);
    if (x_gram(dim-1,dim-1) < min_dist) {
      min_dist = x_gram(dim-1,dim-1);
      min_g = g;
      x_closest = x;
    }
  }
  
#ifdef DEBUG_LEVEL_FULL
  std::cerr << "x_closest = " << std::endl;
  x_closest.prettyPrint(std::cerr, dim-1);
#endif
  
  iso = iso*min_g;
  q = min_g.transform(q);

#ifdef DEBUG_LEVEL_FULL
  std::cerr << "returning isometry: " << std::endl;
  iso.a.prettyPrint(std::cerr, dim);
  std::cerr << "transformed gram to: " << std::endl;
  q.prettyPrint(std::cerr, dim);
#endif
  return;
}

// to avoid recursive template instantiation,
// we supply a parameter defining the level of recursion
// and use only this part of the matrices
// All containers will have size n, but we will only use dim entries
template<typename R, size_t n>
inline void QuadFormInt<R,n>::greedy(SquareMatrixInt<R,n>& gram,
			      Isometry<R,n>& s,
			      size_t dim)
{

#ifdef DEBUG_LEVEL_FULL
  Isometry<R,n> s0 = s;
  SquareMatrixInt<R,n> q0 = gram;
#endif
  
  if (dim == 1) return;

  // temp isometry
  Isometry<R,n> temp;

  std::pair<R,size_t> perm_pair[n];
  VectorInt<size_t,n> perm;
  do {
    for (size_t i = 0; i < dim; i++)
      perm_pair[i] = std::make_pair(gram(i,i), i);
    std::sort(perm_pair, perm_pair+dim);
    
    for (size_t i = 0; i < dim; i++)
      perm[i] = perm_pair[i].second;

    // this is to make sure we do not touch these rows
    for (size_t i = dim; i < n; i++)
      perm[i] = i;

    temp.setIdentity();
    temp.updatePerm(perm);

    // update isometry
    // s.update_perm(perm);
    s = s*temp;
    
    // update gram
    gram = temp.transform(gram);
    
#ifdef DEBUG_LEVEL_FULL
    assert((s0.inverse()*s).transform(q0) == gram);
#endif

    // !! - TODO - do we really need iso here
    // or could we simply pass s?
    Isometry<R,n> iso;
	
    greedy(gram, iso, dim-1);

    s = s*iso;
    // !! TODO - one can use subgram to save computations
    // This transformation already happens inside greedy(dim-1)
    //     gram = iso.transform(gram);

#ifdef DEBUG_LEVEL_FULL
    assert((s0.inverse()*s).transform(q0) == gram);
#endif
    
    closestLatticeVector(gram, s, dim);

#ifdef DEBUG_LEVEL_FULL
    assert((s0.inverse()*s).transform(q0) == gram);
#endif
    
  } while (gram(dim-1,dim-1) < gram(dim-2,dim-2));
  return;

}

template<typename R, size_t n>
inline std::vector< std::vector<size_t> >
QuadFormInt<R,n>::allPerms(size_t m)
{
  std::vector< std::vector<size_t> > perms;
  if (m == 1) {
    std::vector<size_t> id(1);
    id[0] = 0;
    perms.push_back(id);
    return perms;
  }
  std::vector< std::vector<size_t> > rec_perms = allPerms(m-1);
  for (std::vector<size_t> perm : rec_perms) {
    perm.push_back(m-1);
    perms.push_back(perm);
    for (size_t i = 0; i < m-1; i++) {
      std::vector<size_t> perm_new(m);
      for (size_t j = 0; j < m-1; j++)
	perm_new[j] = perm[j];
      perm_new[m-1] = perm[i];
      perm_new[i] = m-1;
      perms.push_back(perm_new);
    }
    
  }
  return perms;
}

template<typename R, size_t n>
inline bool QuadFormInt<R,n>::permutationReduction(SquareMatrixInt<R,n> & qf,
					     Isometry<R,n> & isom,
					     std::set< Isometry<R, n> > & auts,
					     bool calc_aut)
{
  bool is_reduced = true;
  std::map<Integer<R>, std::vector<size_t> > stable_sets;
  Isometry<R, n> s_final;
  SquareMatrixInt<R,n> q0, q1;
  q0 = qf;
  
  for (size_t i = 0; i < n; i++) {
    Integer<R> val = qf(i,i);
    auto search = stable_sets.find(val);
    if (search == stable_sets.end()) {
      std::vector<size_t> empty_vec;
      stable_sets[val] = empty_vec;
    }
    stable_sets[val].push_back(i);
  }
  // !! TODO - Here I go one by one, but in magma
  // Kohel tries all possible permutations (all products)
  // Could it really matter?
  typename std::map<Integer<R>, std::vector<size_t> >::const_iterator iter;
  for (iter = stable_sets.begin(); iter != stable_sets.end(); iter++) {
    //    R key = iter->first;
    std::vector<size_t> value = iter->second;
    std::vector< std::vector<size_t> > val_perms = allPerms(value.size());
    for (std::vector<size_t> perm : val_perms) {
      VectorInt<size_t, n> large_perm;
      for (size_t k = 0; k < n; k++)
	large_perm[k] = k;
      for (size_t idx = 0; idx < value.size(); idx++) {
	large_perm[value[idx]] = value[perm[idx]];
      }
      Isometry<R,n> s;
      s.updatePerm(large_perm);
      q1 = s.transform(qf);
      if (q1 < q0) {
	q0 = q1;
	s_final = s;
	is_reduced = false;
      }
      else if ((calc_aut) && (q1 == q0)) {
	auts.insert(isom*s*s_final.inverse()*isom.inverse());
      }
    }
  }
  qf = q0;
  isom = isom*s_final;
  return is_reduced;
}

template<typename R, size_t n>
inline bool QuadFormInt<R,n>::signNormalizationSlow(SquareMatrixInt<R,n> & qf,
					     Isometry<R,n> & isom,
					     std::set< Isometry<R,n> > & auts)
{
  bool is_reduced = true;
  W16 prime = 2;
  std::random_device rd;
  W64 seed = rd();
  std::shared_ptr<W16_F2> GF2 = std::make_shared<W16_F2>(prime,seed);
  // !! - TODO - we don't really need set here, can use vector
  std::set< W16_VectorFp<n> > boundary_basis;
  std::set< std::pair<size_t, size_t> > priority_set;
  
  size_t count = 0;
  for (size_t j = 1; j < n; j++)
    for (size_t k = 0; k < n-j; k++) {
      W16_MatrixFp w_F2(GF2, boundary_basis.size()+1, n);
      typename std::set< W16_VectorFp<n> >::const_iterator bb_ptr;
      bb_ptr = boundary_basis.begin();
      for (size_t row = 0; row < boundary_basis.size(); row++) {
	for (size_t col = 0; col < n; col++)
	  w_F2(row, col) = (*bb_ptr)[col];
	bb_ptr++;
      }
      for (size_t col = 0; col < n; col++)
	w_F2(boundary_basis.size(), col) = GF2->mod(0);
      w_F2(boundary_basis.size(), k) = GF2->mod(1);
      w_F2(boundary_basis.size(), k+j) = GF2->mod(1);
      if ((w_F2.rank() > count) && (qf(k,k+j) != 0)) {
	priority_set.insert(std::make_pair(k,k+j));
	W16_VectorFp<n> last_row(GF2);
	for (size_t col = 0; col < n; col++)
	  last_row[col] = GF2->mod(0);
	last_row[k] = GF2->mod(1);
	last_row[k+j] = GF2->mod(1);
	boundary_basis.insert(last_row);
	count++;
      }
    }

  W16_MatrixFp w_F2(GF2, priority_set.size(), n+1);
  std::set< std::pair<size_t, size_t> >::const_iterator ps_ptr;
  ps_ptr = priority_set.begin();
  for (size_t row = 0; row < priority_set.size(); row++) {
    for (size_t col = 0; col <= n; col++)
	w_F2(row, col) = GF2->mod(0);
    
    w_F2(row, ps_ptr->first) = GF2->mod(1);
    w_F2(row, ps_ptr->second) = GF2->mod(1);

    // the affine coordinate
    if (qf(ps_ptr->first, ps_ptr->second) < 0)
      w_F2(row, n) = GF2->mod(1);
    
    ps_ptr++;
  }

  W16_MatrixFp ker = w_F2.kernel();
  // The last row of ker should now be a solution to the affine equation
  // The rows above are the kernel
#ifdef DEBUG
  for (size_t row = 0; row + 1 < ker.nrows(); row++)
    assert(ker(row, n) == GF2->mod(0));
  if (ker.nrows() >= 1)
    assert(ker(ker.nrows()-1, n) == GF2->mod(1));
#endif
  //  W16_MatrixFp tmp(GF2, ker.nrows(), ker.nrows());
  // W16_MatrixFp::row_echelon(ker, tmp);
  Isometry<R,n> s;
  is_reduced = true;
  for (size_t row = 0; row + 1 < ker.nrows(); row++) {
    is_reduced = false;
    for (size_t i = 0; i < n; i++)
      s(i,i) = (ker(row, i) + ker(ker.nrows()-1, i) == 1) ? -1 : 1;
    if (s.transform(qf) == qf) {
      auts.insert(isom*s*isom.inverse());
      is_reduced = true;
      // to be compatible with magma implementation for debugging
      for (size_t i = 0; i < n; i++) s(i,i) = 1;
    }
  }
  qf = s.transform(qf);
  isom = isom*s;
  return is_reduced;
}

// Returns the matrix obtained by a basis permutation 
// such that QF[i,i] le QF[j,j] for all i le j. }
// !! TODO - maybe use this one also to update the automorphism group
// (by permutation). Though it seems this is covered by permutation_reduction
template<typename R, size_t n>
inline bool QuadFormInt<R,n>::normEchelon(SquareMatrixInt<R,n> & qf,
				   Isometry<R,n> & isom)
{
#ifdef DEBUG
  SquareMatrixInt<R,n> qf_orig = qf;
#endif
  bool is_reduced = true;
  Isometry<R,n> s, u0;
  for (size_t i = 0; i < n-1; i++) {
    if (qf(i+1,i+1) < qf(i,i)) {
      s.setIdentity();
      s(i+1, i+1) = 0;
      s(i, i) = 0;
      s(i,i+1) = 1;
      s(i+1, i) = 1;
      qf = s.transform(qf);
      u0 = u0*s;
#ifdef DEBUG
      assert(u0.transform(qf_orig) == qf);
#endif
      is_reduced = false;
    }
  }
  if (u0.a != SquareMatrixInt<R,n>::identity())
    is_reduced = (is_reduced) && norm_echelon(qf, isom);
  isom = isom*u0;
  return is_reduced;
}

template<typename R, size_t n>
inline bool QuadFormInt<R,n>::neighborReduction(SquareMatrixInt<R,n> & qf,
					 Isometry<R,n> & isom,
					 std::set< Isometry<R,n> > & auts,
					 bool calc_aut)
{
#ifdef DEBUG_LEVEL_FULL
  SquareMatrixInt<R,n> qf_orig = qf;
  Isometry<R,n> isom_orig = isom;
#endif
  bool is_reduced = true;
  std::vector< std::set< Vector<Integer<R>, n> > > local_neighbors(1);
  Isometry<R,n> b0;
  VectorInt<R,n> vec;
  vec[0] = 1;
  for (size_t i = 1; i < n; i++)
    vec[i] = 0;
  local_neighbors[0].insert(vec);
  size_t num_free = 1;
  for (size_t i = 1; i < n; i++) {
    num_free *= 3;
    std::set< Vector<R, n> > free_hood;
    for (size_t x_idx = 0; x_idx < num_free; x_idx++) {
      size_t tmp = x_idx;
      VectorInt<R,n> x;
      for (size_t j = 0; j < i; j++) {
	// we separate because tmp is unsigned, which might lead to overflow
	x[j] = (tmp % 3);
	x[j]--;
	tmp /= 3;
      }
      x[i] = 1;
      for (size_t j = i+1; j < n; j++)
	x[j] = 0;
      Integer<R> norm = VectorInt<R,n>::innerProduct(x*qf, x);
      if (norm < qf(i,i)) {
	// !! TODO - this doesn't really happen
	// when we start with something which is
	// already greedy reduced
	b0.setIdentity();
	for (size_t j = 0; j < n; j++)
	  b0(j,i) = x[j];
	qf = b0.transform(qf);
	isom = isom*b0;
#ifdef DEBUG_LEVEL_FULL
	assert((isom_orig.inverse() * isom).transform(qf_orig) == qf);
#endif
	norm_echelon(qf, isom);
#ifdef DEBUG_LEVEL_FULL
	assert((isom_orig.inverse() * isom).transform(qf_orig) == qf);
#endif
	return false;
      }
      else if (norm == qf(i,i)) {
	free_hood.insert(x);
      }
    }
    local_neighbors.push_back(free_hood);
  }

  std::map< Integer<R>, std::vector<size_t> > norms;
  for (size_t i = 0; i < n; i++) {
    Integer<R> val = qf(i,i);
    auto search = norms.find(val);
    if (search == norms.end()) {
      std::vector<size_t> empty_vec(0);
      norms[val] = empty_vec;
    }
    norms[val].push_back(i);
  }
  typename std::map< Integer<R>, std::vector<size_t> >::const_iterator iter;
  // !! TODO - here there's duplication.
  // we can simply call local_neighbors[norms[i]],
  // changin local_neighbors to depend on the norm
  for (iter = norms.begin(); iter != norms.end(); iter++) {
    std::vector<size_t> inds = iter->second;
    std::set< VectorInt<R,n> > X;
    for (size_t i : inds) {
      X.insert(local_neighbors[i].begin(), local_neighbors[i].end());
    }
    for (size_t i  : inds)
      local_neighbors[i] = X;
  }

  size_t nbs_size = 1;
  for (size_t i = 0; i < local_neighbors.size(); i++)
    nbs_size *= local_neighbors[i].size();
  
#ifdef DEBUG_LEVEL_FULL
  std::cerr << "Original NeighborSpace size: " << nbs_size << std::endl;
#endif
  
  std::vector< std::vector< VectorInt<R,n> > > neighbor_space;
  for (VectorInt<R,n> x : local_neighbors[0]) {
    std::vector< VectorInt<R,n> > singleton;
    singleton.push_back(x);
    neighbor_space.push_back(singleton);
  }
  for (size_t i = 1; i < n; i++) {
    std::vector< std::vector< VectorInt<R,n> > > ns1;
    R norm = qf(i,i);
    std::vector<size_t> inds;
    for (size_t j = 0; j < i; j++)
      if (qf(j,j) == norm) inds.push_back(j);
    for (VectorInt<R,n> y : local_neighbors[i]) {
      for (std::vector< VectorInt<R,n> > c : neighbor_space) {
	bool include = true;
	for (size_t j : inds)
	  if (c[j] == y) {
	    include = false;
	    break;
	  }
	if ((include) &&
	    (abs(VectorInt<R,n>::innerProduct(c[i-1]*qf, y)) >= abs(qf(i,i-1)))) {
	  c.push_back(y);
	  ns1.push_back(c);
	}
	else {
	  for (size_t j = 1; j < i; j++) {
	    if (abs(VectorInt<R,n>::innerProduct(c[j-1]*qf, c[j])) >
		abs(qf(j,j-1))) {
	      c.push_back(y);
	      ns1.push_back(c);
	      break;
	    }
	  }
	}
      }
    }
    neighbor_space = ns1;
  }
  
#ifdef DEBUG_LEVEL_FULL
  std::cerr << "Reduced to NeighborSpace of size:";
  std::cerr << neighbor_space.size() << std::endl;
#endif
  
  // !! - TODO - we can from the beginning store c as a matrix
  for (std::vector< VectorInt<R,n> > c : neighbor_space) {
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < n; j++)
	// note that we transpose
	b0(i,j) = c[j][i];
    if (abs(b0.a.determinant()) == 1) {
      SquareMatrixInt<R,n> q0 = b0.transform(qf);
      Isometry<R,n> u;
      std::set< Isometry<R,n> > tmp_auts;
      signNormalization(q0, u, tmp_auts, calc_aut);
      if (q0 < qf) {
	qf = q0;
        isom = isom*b0*u;
	is_reduced = false;

#ifdef DEBUG_LEVEL_FULL
	assert((isom_orig.inverse() * isom).transform(qf_orig) == qf);
#endif
	return false;
      }
      else if ((calc_aut) && (q0 == qf)) {
	auts.insert(isom*b0*u*isom.inverse());
      }
    }
  }
  return is_reduced;
}

// This generates the entire automorphism group.
// We don't really need to do this.
// implement a small version of Todd-Coxeter Schreier-Sims ?!
template<typename R, size_t n>
inline size_t QuadFormInt<R,n>::generateAuts(std::set< Isometry<R,n> > & auts)
{
  size_t num_aut;
  do {
    num_aut = auts.size();
    typename std::set< Isometry<R,n> >::const_iterator iter1, iter2;
    for (iter1 = auts.begin(); iter1 != auts.end(); iter1++) {
      for (iter2 = auts.begin(); iter2 != auts.end(); iter2++) {
	Isometry<R,n> prod = (*iter1)*(*iter2);
	if (auts.find(prod) == auts.end())
	  auts.insert(prod);
      }
    }
    // if this condition is fullfilled we are closed under taking
    // products. Since this is a finite group, we are done.
  } while (num_aut != auts.size());
  return num_aut;
}

template<typename R, size_t n>
inline size_t QuadFormInt<R,n>::numAutomorphisms() const
{
  if (this->_num_aut_init) return this->_num_aut;
  SquareMatrixInt<R,n> qf = this->_B;
  Isometry<R,n> isom;
  std::set< Isometry<R,n> > auts;
  return iReduce(qf, isom, auts, true);
}

// !! - TODO - think whether we want to save this as a member.
// Right now it seems to me that most of the time we don't need it,
// so there is no need to carry it around.
template<typename R, size_t n>
inline std::set<Isometry<R,n>> QuadFormInt<R,n>::properAutomorphisms() const
{
  SquareMatrixInt<R,n> qf = this->_B;
  Isometry<R,n> isom;
  std::set< Isometry<R,n> > auts;
  iReduce(qf, isom, auts, true);
  return auts;
}

template<typename R, size_t n>
inline QuadFormInt<R,n> QuadFormInt<R,n>::reduce(const QuadFormInt<R,n> & q,
					  Isometry<R,n> & isom,
					  bool calc_aut)
{
  assert(q.bilinearForm().isPositiveDefinite());

  std::set< Isometry<R,n> > auts;
  SquareMatrixInt<R,n> qf = q.bilinearForm();
  size_t num_aut = iReduce(qf, isom, auts, calc_aut);
  QuadFormInt<R,n> q_red(qf);
  if (calc_aut) {
    q_red._num_aut = num_aut;
    q_red._num_aut_init = true;
  }
  q_red._is_reduced = true;
  return q_red;
}

template<typename R, size_t n>
inline size_t QuadFormInt<R,n>::iReduce(SquareMatrixInt<R,n> & qf,
					Isometry<R,n> & isom,
					std::set< Isometry<R,n> > & auts,
					bool calc_aut)
{
#ifdef DEBUG_LEVEL_FULL
  SquareMatrixInt<R,n> q0 = qf;
  Isometry<R,n> s0 = isom;
#endif
  greedy(qf, isom);
#ifdef DEBUG_LEVEL_FULL
  assert((s0.inverse()*isom).transform(q0) == qf);
#endif
  
  bool is_reduced;
  do {
    is_reduced = true;
    is_reduced = (permutationReduction(qf, isom, auts, calc_aut)) && (is_reduced);
#ifdef DEBUG_LEVEL_FULL
    assert((s0.inverse()*isom).transform(q0) == qf);
    for (Isometry<R, n> s : auts) {
      assert((s0.inverse()*s*s0).transform(q0) == q0);
    }
#endif    
    is_reduced = (signNormalization(qf, isom, auts, calc_aut)) && (is_reduced);
#ifdef DEBUG_LEVEL_FULL
    assert((s0.inverse()*isom).transform(q0) == qf);
    for (Isometry<R,n> s : auts) {
      assert((s0.inverse()*s*s0).transform(q0) == q0);
    }
#endif
    
    is_reduced = (neighborReduction(qf, isom, auts, calc_aut)) && (is_reduced);
    
#ifdef DEBUG_LEVEL_FULL
    assert((s0.inverse()*isom).transform(q0) == qf);
    for (Isometry<R,n> s : auts) {
      assert((s0.inverse()*s*s0).transform(q0) == q0);
    }
#endif
  } while (!is_reduced);
  if (!calc_aut)
    return auts.size();
  return generateAuts(auts);
}

template<typename R, size_t n>
inline bool QuadFormInt<R,n>::signNormalization(SquareMatrixInt<R,n> & qf,
						Isometry<R,n> & isom,
						std::set< Isometry<R,n> > & auts,
						bool calc_aut)
{
  if (!calc_aut)
    return signNormalizationFast(qf, isom);
  return signNormalizationSlow(qf, isom, auts);
}

template<typename R, size_t n>
inline std::vector<uint8_t>
QuadFormInt<R,n>::bitTranspose(const std::vector< uint8_t > & mat)
{
  assert(n+1 <= 8);

  std::vector<uint8_t> trans(n+1);

  for (uint8_t row = 0; row <= n; row++) {
    trans[row] = 0;
    for (uint8_t col = 0; col < mat.size(); col++) {
      trans[row] |= (((mat[col] >> row) & 1) << col);
    }
  }

  return trans;
}

// returns the transformation and the rank
// performs the echelonization in place
template<typename R, size_t n>
inline uint8_t QuadFormInt<R,n>::bitEchelonForm(std::vector< uint8_t > & mat,
						std::vector< uint8_t > & trans)
{
  assert(mat.size() <= 8);
  
  trans.resize(mat.size());
  
  for (size_t row = 0; row < mat.size(); row++)
    trans[row] = (1 << row);
  
  uint8_t pivot_row;
  pivot_row = 0;
  uint8_t pivot_col;
  pivot_col = 0;
 
  uint8_t row;  
  uint8_t val;
  
  while ((pivot_row < mat.size()) && (pivot_col <= n)) {
    val = 0;
    for (row = pivot_row ; (!val) && (row < mat.size()); row++) {
      val = (mat[row] >> pivot_col) & 1;
    }
    if (!val) {
      pivot_col++;
    }
    else {
      row--;
      if (row != pivot_row) {
	// swapping rows
	mat[pivot_row] ^= mat[row];
	mat[row] ^= mat[pivot_row];
	mat[pivot_row] ^= mat[row];
      
	trans[pivot_row] ^= trans[row];
	trans[row] ^= trans[pivot_row];
	trans[pivot_row] ^= trans[row];
      }
      
      for (row = pivot_row+1; row < mat.size(); row++) {
	val = (mat[row] >> pivot_col) & 1;
	if (val) {
	  mat[row] ^= mat[pivot_row];
	  trans[row] ^= trans[pivot_row];
	}
      }
      
      pivot_row++;
      pivot_col++;
    }
  }
  return pivot_row;
}

template<typename R, size_t n>
inline std::vector<uint8_t>
QuadFormInt<R,n>::kernel(const std::vector< uint8_t > & mat)
{
  std::vector<uint8_t> ker;

  std::vector<uint8_t> mat_t = bitTranspose(mat);
  
  uint8_t rank = bitEchelonForm(mat_t, ker);
  // getting the zero rows

  ker.erase(ker.begin(), ker.begin() + rank);
  
  return ker;
}

// !! TODO - use bit slicing to make this faster
// Also does not need to compute the rank every time
template<typename R, size_t n>
bool QuadFormInt<R,n>::signNormalizationFast(SquareMatrixInt<R,n> & qf,
					     Isometry<R,n> & isom)
{
  bool is_reduced = true;

  // This assumes n < 8
  assert(n < 8);

  std::vector< uint8_t > bb_vecs;
  uint8_t vec;
  // vec after reduction - for echelon form
  uint8_t ech_vec;
  
  // There should be a more efficient way of doing this,
  // but it will help me keep track of things for now
  // save the pivots of each row, this is always sorted
  std::vector<uint8_t> pivots;
  
  // for each k (from 0 to n-1) save the row number
  // in which we will want it to be placed
  // If it is already a pivot, this will be -1
  int8_t place_pivots[n] = {0};

  // the position of the row where k is s pivot.
  // If it is not, then it is -1.
  int8_t inv_pivots[n];
  for (uint8_t i = 0; i < n-1; i++) inv_pivots[i] = -1;
  
  for (size_t j = 1; j < n; j++) {
    // vec will always have the only the bits k and k+j on 
    vec = 1 | (1 << j);
    for (size_t k = 0; k < n-j; k++) {
      if (qf(k,k+j) != 0) {
	int lead = k;
	ech_vec = vec;

	if (qf(k,k+j) < 0)
	  ech_vec |= (1 << n);

	// while we already have this as pivot, we 
	while ((lead >= 0) && (inv_pivots[lead] >= 0) && (lead != n)) {
	  ech_vec ^= bb_vecs[inv_pivots[lead]];
	  lead = ffs(ech_vec)-1;
	}

	// If it is not a pivot, we put it in its proper place
	// and update the arrays tracking the pivots
	if ((lead >= 0) && (lead != n)) {
	  // echelonize the rows already in bb_vecs above the pivot
	  for (uint8_t row = 0; row < place_pivots[lead]; row++) {
	    uint8_t bit = (bb_vecs[row] >> lead) & 1;
	    if (bit) bb_vecs[row]^= ech_vec;
	  }
	  bb_vecs.insert(bb_vecs.begin()+place_pivots[lead], ech_vec);
	  inv_pivots[lead] = place_pivots[lead];
	  place_pivots[lead] = -1;
	  // next time we will get any of these, we will want
	  // them after this one
	  for (size_t r = lead+1; r < n; r++) {
	    place_pivots[r]++;
	    if (inv_pivots[r] >= 0)
	      inv_pivots[r]++;
	  } 
	}
      }
      vec <<= 1;
    }
  }

  // The last row of ker should now be a solution to the affine equation
  // The rows above are the kernel
  std::vector<uint8_t> ker_bit = kernel(bb_vecs);

#ifdef DEBUG_LEVEL_FULL
  Isometry<R,n> s;
  for (size_t i = 0; i < n; i++)
    s(i,i) = ((ker_bit[ker_bit.size()-1] >> i) & 1) ? -1 : 1;
  if (s.transform(qf) == qf) {
    is_reduced = true;
    // to be compatible with magma implementation for debugging
    for (size_t i = 0; i < n; i++) s(i,i) = 1;
  }
  SquareMatrixInt<R,n> qf_orig = qf;
  Isometry<R,n> isom_orig = isom;
#endif 
  is_reduced = true;
  
  if (!ker_bit.empty()) {
    for (size_t row = 0; row+1 < n; row++) {
      uint8_t bit_row = (ker_bit[ker_bit.size()-1] >> row) & 1;
      for (size_t col = row+1; col < n; col++) {
	uint8_t bit_col = (ker_bit[ker_bit.size()-1] >> col) & 1;
	if (bit_col^bit_row) {
	  qf(row,col) = -qf(row,col);
	  qf(col,row) = -qf(col,row);
	  is_reduced = false;
	}
      }
    }
#ifdef DEBUG_LEVEL_FULL
    assert(s.transform(qf_orig) == qf);
#endif
    if (!is_reduced)
      for (size_t col = 0; col < n; col++)
	if ((ker_bit[ker_bit.size()-1] >> col) & 1) {
	  for (size_t row = 0; row < n; row++)
	    isom(row,col) = -isom(row,col);
	}
#ifdef DEBUG_LEVEL_FULL
    assert(isom_orig*s == isom);
#endif
  }
  
  return is_reduced;
}

template<typename R, size_t n>
inline std::unordered_map<QuadFormInt<R,n>, Isometry<R,n> >
QuadFormInt<R,n>::generateOrbit() const
{
  Isometry<R,n> s;
  QuadFormInt<R,n> qf(this->bilinearForm());
  size_t num = 0;
  std::unordered_map< QuadFormInt<R,n>, Isometry<R,n> > orbit;
  typename std::unordered_map< QuadFormInt<R,n>,
			       Isometry<R,n> >::const_iterator i, j;
  orbit.insert(std::make_pair(qf, s));
  while (num < orbit.size()) {
    num = orbit.size();
    for (i = orbit.begin(); i != orbit.end(); i++) {
      std::unordered_map< QuadFormInt<R,n>, Isometry<R,n> >
	perms = (i->first).permutationOrbit();
      for (j = perms.begin(); j != perms.end(); j++) {
	assert((i->second*j->second).transform(this->bilinearForm()) ==
	       j->first.bilinearForm());
        orbit[j->first] = i->second*j->second;
      }
    }
    for (i = orbit.begin(); i != orbit.end(); i++) {
      std::unordered_map< QuadFormInt<R,n>, Isometry<R,n> >
	signs = (i->first).signOrbit();
      for (j = signs.begin(); j != signs.end(); j++) {
	assert((i->second*j->second).transform(this->bilinearForm()) ==
	       j->first.bilinearForm());
        orbit[j->first] = i->second*j->second;
      }
    }
  }
  return orbit;
}


template<typename R, size_t n>
inline std::unordered_map< QuadFormInt<R,n>, Isometry<R,n> >
QuadFormInt<R,n>::permutationOrbit() const
{
  std::unordered_map< QuadFormInt<R,n>, Isometry<R,n> > orbit; 
  std::map<R, std::vector<size_t> > stable_sets;
  
  SquareMatrixInt<R,n> q1;
  
  for (size_t i = 0; i < n; i++) {
    Integer<R> val = this->bilinearForm()(i,i);
    auto search = stable_sets.find(val);
    if (search == stable_sets.end()) {
      std::vector<size_t> empty_vec;
      stable_sets[val] = empty_vec;
    }
    stable_sets[val].push_back(i);
  }
  
  typename std::map<Integer<R>, std::vector<size_t> >::const_iterator iter;
  for (iter = stable_sets.begin(); iter != stable_sets.end(); iter++) {
    std::vector<size_t> value = iter->second;
    std::vector< std::vector<size_t> > val_perms =
      QuadFormInt<R,n>::allPerms(value.size());
    for (std::vector<size_t> perm : val_perms) {
      VectorInt<size_t,n> large_perm;
      for (size_t k = 0; k < n; k++)
	large_perm[k] = k;
      for (size_t idx = 0; idx < value.size(); idx++) {
	large_perm[value[idx]] = value[perm[idx]];
      }
      Isometry<R,n> s;
      s.updatePerm(large_perm);
      q1 = s.transform(this->bilinearForm());
      greedy(q1, s);
      QuadFormInt<R,n> q(q1);

      assert(s.transform(this->bilinearForm()) == q.bilinearForm());

      orbit[q] = s;
    }
  }
  return orbit;
}

template<typename R, size_t n>
inline std::unordered_map< QuadFormInt<R,n>, Isometry<R,n> >
QuadFormInt<R,n>::signOrbit() const
{
  std::unordered_map< QuadFormInt<R,n>, Isometry<R,n> > orbit;
  Isometry<R,n> s;
  SquareMatrixInt<R,n> q;
  
  for (size_t signs = 0; signs < (1 << n); signs++) {
    size_t tmp = signs;
    for (size_t j = 0; j < n; j++) {
      s(j,j) = ((tmp & 1) ? -s(j,j) : s(j,j));
      tmp >>= 1;
    }
    q = s.transform(this->bilinearForm());
    greedy(q,s);
    QuadFormInt<R,n> qq(q);
    assert(s.transform(this->bilinearForm()) == qq.bilinearForm());
    orbit[qq] = s;
  }
  return orbit;
}

// general hashValue

template<size_t n>
inline W64 Z_QuadForm<n>::hashValue(void) const
{
  W64 fnv = FNV_OFFSET;
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j <= i; j++)
      fnv = (fnv ^ mpz_get_si((this->_B(i,j)).get_mpz_t())) * FNV_PRIME;

  return fnv;
}

template<size_t n>
inline W64 Z64_QuadForm<n>::hashValue(void) const
{
  W64 fnv = FNV_OFFSET;
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j <= i; j++)
      fnv = (fnv ^ this->_B(i,j)) * FNV_PRIME;
  return fnv;
}

template<size_t n>
inline W64 Z128_QuadForm<n>::hashValue(void) const
{
  W64 fnv = FNV_OFFSET;
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j <= i; j++)
      fnv = (fnv ^ this->_B(i,j)) * FNV_PRIME;
  return fnv;
}
