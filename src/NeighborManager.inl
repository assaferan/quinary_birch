#include "Polynomial.h"
#include "PolynomialRing.h"
#include "QuadFormFp.h"
#include "SquareMatrixInt.h"
#include "VectorInt.h"

template<typename R, typename S, typename T, size_t n>
inline NeighborManager<R,S,T,n>::NeighborManager(const QuadFormZZ<T,n>& q,
						 std::shared_ptr<Fp<R,S>> GF,
						 size_t k)
  : _vec(GF), _b(GF), _quot_gram()
{
  T p = GF->prime();
  
  this->_q = q;
  this->_disc = q.discriminant();

  std::shared_ptr<QuadFormFp<R,S,n> > qp = q.mod(GF);

  this->_b = qp->bilinearForm();
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      this->_quot_gram(i,j) = (this->_q(i,j)) % (p*p);
	
  this->_GF = GF;
  assert(qp->isotropicVector(this->_vec));

#ifdef DEBUG
  R prime = GF->prime();
  if (prime != 2) assert( qp->evaluate(this->_vec) == 0 );
#endif

  this->_p_std_gram = std::make_shared<SquareMatrixFp<R,S,n> >(GF);
  this->_p_basis = std::make_shared<SquareMatrixFp<R,S,n> >(GF);

#ifdef DEBUG
  qp->decompose(*this->_p_std_gram, *this->_p_basis, true);
#else
  qp->decompose(*this->_p_std_gram, *this->_p_basis);
#endif
  
  this->_p_q_std = std::make_shared<PolynomialFp<R,S> >(*this->_p_std_gram);

#ifdef DEBUG_LEVEL_FULL
  std::cerr << "Performed Witt Decomposition on" << std::endl;
  qp->bilinearForm().prettyPrint(std::cerr);
  std::cerr << "Resulting gram matrix is " << std::endl;
  this->_p_std_gram->prettyPrint(std::cerr);
  std::cerr << "Resulting basis is " << std::endl;
  this->_p_basis->prettyPrint(std::cerr);
#endif

  // Count the rows at the end of the matrix which are exactly zero.
  size_t idx = n;
  while ((idx >= 1) && (*this->_p_std_gram)[idx-1].isZero()) {
    idx--;
  }

  // The dimension of the radical.
  this->_rad_dim = n - idx;

  // Determine the dimension of the totally hyperbolic subspace.
  idx = 1;
  while ((idx <= n-_rad_dim) && ((*_p_std_gram)(idx-1, idx-1) == 0) )
    idx++;

  // Dimension of the anistotropic subspace.
  this->_aniso_dim = n - _rad_dim - idx + 1;

  // The number of hyperbolic planes in the Witt decomposition.
  this->_witt_index = (idx - 1) / 2;

  this->_pivots = __pivots(n-_rad_dim, _aniso_dim, k);
  this->_pivot_ptr = 0;
  this->_k = k;
  this->_skew_dim = k*(k-1)/2;
  this->_p_skew = std::make_shared< MatrixFp<R,S> >(this->_GF, k, k);

  this->nextIsotropicSubspace();
  this->_liftSubspace();
}

//!! TODO - make gram work only modulo p^2

template<typename R, typename S, typename T, size_t n>
inline SquareMatrixInt<T,n>
NeighborManager<R,S,T,n>::__gram(const SquareMatrixInt<T,n> & B, bool quot) const
{
  T p = this->_GF->prime();
  std::shared_ptr<const IntegerRing<T> > ZZ = std::make_shared<const IntegerRing<T> >();
  SquareMatrixInt<T,n> gram(ZZ);
  
  if (p == 2)
    gram = B * this->_q.bilinearForm() * B.transpose();
  else
    gram =  B * (this->_quot_gram) * B.transpose();
  
  // !! TODO - this isn't necessary only for debugging versus magma
  if ((quot) || (p != 2))
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < n; j++)
	gram(i,j) = gram(i,j) % (p*p);
  
  return gram;
}

template<typename R, typename S, typename T, size_t n>
inline void NeighborManager<R,S,T,n>::_liftSubspace(void)
{
  if ((this->_iso_subspace).empty()) return;
  Integer<T> p = this->_GF->prime();

  assert(this->_pivot_ptr >= 1);
  
  // Get the pivots for the bases of the isotropic subspaces.
  std::vector<size_t> pivots = this->_pivots[this->_pivot_ptr-1];

#ifdef DEBUG_LEVEL_FULL
  std::cerr << "before lifting, p_basis is " << std::endl << (*this->_p_basis);
  std::cerr << std::endl;
  std::cerr << "iso_subspace is " << this->_iso_subspace << std::endl;
  std::cerr << "pivots = ";
  for (size_t i = 0; i < pivots.size(); i++)
    std::cerr << pivots[i];
  std::cerr << std::endl;
#endif

  SquareMatrixFp<R,S,n> basis = *this->_p_basis;
  // Set up the correct basis vectors.
  for (size_t i = 0; i < this->_k; i++)
    for (size_t j = pivots[i]+1; j < n; j++)
      basis.addCol(pivots[i], j, this->_iso_subspace[i][j]);

#ifdef DEBUG_LEVEL_FULL
  std::cerr << "the correct basis vectors are" << std::endl << basis;
  std::cerr << std::endl;
#endif
  
  // Extract our target isotropic subspace modulo p
  std::vector< VectorFp<R,S,n> > x,z,u;
  for (size_t i = 0; i < this->_k; i++) {
    x.push_back(basis.transpose()[pivots[i]]);
  }

#ifdef DEBUG_LEVEL_FULL
  std::cerr << "x = " << x << std::endl;
#endif
  
  // Extract the hyperbolic complement modulo pR.
  std::vector<size_t> paired(this->_k);
  size_t h_dim = 2 * this->_witt_index; 
  for (size_t i = 0; i < this->_k; i++)
    paired[i] = h_dim - 1 - pivots[this->_k-1-i];
  for (size_t i = 0; i < this->_k; i++) {
    z.push_back(basis.transpose()[paired[i]]);
  }

#ifdef DEBUG_LEVEL_FULL
  std::cerr << "z = " << z << std::endl;
#endif
  
  // Extract the remaining basis vectors.
  std::set<size_t> exclude;
  exclude.insert(pivots.begin(), pivots.end());
  exclude.insert(paired.begin(), paired.end());
  for (size_t i = 0; i < n; i++) {
    std::set<size_t>::const_iterator iter = exclude.find(i);
    if (iter == exclude.end())
      u.push_back(basis.transpose()[i]);
  }

#ifdef DEBUG_LEVEL_FULL
  std::cerr << "u = " << u << std::endl;
#endif

  std::shared_ptr<const IntegerRing<T> > ZZ = std::make_shared<const IntegerRing<T> >();
  
  // Convert to coordinates modulo p^2.
  this->_X.resize(this->_k);
  this->_Z.resize(this->_k);
  this->_U.resize(n - 2*this->_k);
  
  // Build the coordinate matrix.
  // !! TODO - the mod p is not necessary, good for debugging  
  SquareMatrixInt<T,n> B;
  for (size_t i = 0; i < this->_k; i++)
    for (size_t j = 0; j < n; j++)
      B(i,j) = (this->_X[i][j] = (x[i][j].lift()) % p.num());

  for (size_t i = 0; i < this->_k; i++)
    for (size_t j = 0; j < n; j++)
      B(this->_k+i,j) = (this->_Z[i][j] = (z[i][j].lift()) % p.num());

  for (size_t i = 0; i < n - 2*this->_k; i++)
    for (size_t j = 0; j < n; j++)
      B(2*this->_k+i,j) = (this->_U[i][j] = (u[i][j].lift()) % p.num());

#ifdef DEBUG_LEVEL_FULL
  std::cerr << "X = " << this->_X << std::endl;
  std::cerr << "Z = " << this->_Z << std::endl;
  std::cerr << "U = " << this->_U << std::endl;
#endif
  
  // Compute the Gram matrix of the subspace with respect to the spaces
  //  we will perform the following computations upon.

  SquareMatrixInt<T,n> gram = __gram(B);

  // Lift Z so that it is in a hyperbolic pair with X modulo p^2.
  std::vector< VectorInt<T,n> > Z_new(this->_k);
  for (size_t i = 0; i < this->_k; i++) {
    Z_new[i] = this->_Z[i];
    for (size_t j = 0; j < this->_k; j++) {
      Integer<T> delta = (i == j) ? Integer<T>::one() : Integer<T>::zero();
      // we split the operation due to signed type issues
      Integer<T> a = gram(this->_k-j-1, i + this->_k);
      // a nonnegative value with the same residue mod p*p
      // !!! TODO - we might be able to get rid of that
      // since T is always signed - check!
      a = (a / (p*p) + Integer<T>::one())*p*p-a+delta;
      if (a >= (p*p))
	a -= p*p;
      Z_new[i] += a.num() * this->_Z[j];
    }
  }
  this->_Z = Z_new;
  
#ifdef DEBUG
  for (size_t i = 0; i < this->_k; i++)
    for (size_t j = 0; j < n; j++)
      this->_Z[i][j] = this->_Z[i][j] % (p*p).num();
  
#ifdef DEBUG_LEVEL_FULL
  std::cerr << "after setting <X,Z> = 1" << std::endl;
  std::cerr << "Z = " << this->_Z << std::endl;
#endif // DEBUG_LEVEL_FULL
  
  // Verify that X and Z form a hyperbolic pair.
  // Compute the Gram matrix thusfar.
  for (size_t i = 0; i < this->_k; i++)
    for (size_t j = 0; j < n; j++)
      B(this->_k+i,j) = this->_Z[i][j];
  
  SquareMatrixInt<T,n> temp = __gram(B);
  for (size_t i = 0; i < this->_k; i++)
    for (size_t j = 0; j < this->_k; j++)
      // This is beacuse negative % is negative
      assert((temp(i, this->_k+j) - ((i+j == this->_k-1) ? 1 : 0)) % (p*p).num() == 0);	
#endif // DEBUG
  
  if (p == 2) {
    for (size_t i = 0; i < this->_k; i++)
      for (size_t j = 0; j < n; j++)
	B(this->_k+i,j) = this->_Z[i][j];
    gram = __gram(B, false);
  }
  // Lift X so that it is isotropic modulo p^2.
  std::vector< VectorInt<T,n> > X_new(this->_k);
  Integer<T> two = T(2);
  Integer<T> half = (p*p+Integer<T>::one())/two;
  for (size_t i = 0; i < this->_k; i++) {
    X_new[i] = this->_X[i];
    Integer<T> gram2 = gram(i,i)/2 + ((gram(i,i) % 2 == 0) ? 0 : half.num());
    for (size_t j = this->_k-1-i; j < this->_k; j++) {
      Integer<T> scalar = (i+j == this->_k-1) ? gram2 : gram(i, this->_k-1-j);
      scalar = (scalar / (p*p) + Integer<T>::one())*p*p-scalar;
      if (scalar >= p*p)
	scalar -= p*p;
      X_new[i] +=  scalar.num()  * this->_Z[j];
    }
  }
  this->_X = X_new;

#ifdef DEBUG
  for (size_t i = 0; i < this->_k; i++)
    for (size_t j = 0; j < n; j++)
      this->_X[i][j] = this->_X[i][j] % (p*p).num();

#ifdef DEBUG_LEVEL_FULL
  std::cerr << "after setting <X,X> = 0" << std::endl;
  std::cerr << "X = " << this->_X << std::endl;
#endif // DEBUG_LEVEL_FULL
  
  // Verify that X is isotropic modulo p^2.
  for (size_t i = 0; i < this->_k; i++)
    for (size_t j = 0; j < n; j++)
      B(i,j) = this->_X[i][j];

  // The Gram matrix on this basis.
  temp = __gram(B);

  // Verify all is well.
  for (size_t i = 0; i < this->_k; i++)
    for (size_t j = 0; j < this->_k; j++)
      assert(temp(i,j) % (p*p).num() == 0);
  
#endif // DEBUG

  // Lift Z so that it is isotropic modulo p^2.
  for (size_t i = 0; i < this->_k; i++) {
    Z_new[i] = this->_Z[i];
    for (size_t j = this->_k-1-i; j < this->_k; j++) {
      Integer<T> scalar = gram(this->_k+i, 2*this->_k-1-j);
      if (i+j == this->_k-1)
	scalar = (scalar/two) + ((scalar % two).isZero() ? Integer<T>::zero() : half);
      scalar = (scalar / (p*p) + Integer<T>::one())*p*p-scalar;
      if (scalar >= p*p)
	scalar -= p*p;
      Z_new[i] += scalar.num() * this->_X[j];
    }
  }
  this->_Z = Z_new;

#ifdef DEBUG
  for (size_t i = 0; i < this->_k; i++)
    for (size_t j = 0; j < n; j++)
      this->_Z[i][j] = this->_Z[i][j] % (p*p).num();

#ifdef DEBUG_LEVEL_FULL
  std::cerr << "after setting <Z,Z> = 0" << std::endl;
  std::cerr << "Z = " << this->_Z << std::endl;
#endif // DEBUG_LEVEL_FULL
  
  // Verify that Z is isotropic modulo p^2.
  for (size_t i = 0; i < this->_k; i++)
    for (size_t j = 0; j < n; j++)
      B(this->_k+i,j) = this->_Z[i][j];

  // The Gram matrix on this basis.
  temp = __gram(B);

  // Verify all is well.
  for (size_t i = 0; i < this->_k; i++)
    for (size_t j = 0; j < this->_k; j++)
      assert(temp(this->_k+i,this->_k+j) % (p*p).num() == 0);
  
#endif //DEBUG
  
  // The Gram matrix thusfar.
  for (size_t i = 0; i < this->_k; i++)
    for (size_t j = 0; j < n; j++)
      B(i,j) = this->_X[i][j];

  for (size_t i = 0; i < this->_k; i++)
    for (size_t j = 0; j < n; j++)
      B(this->_k+i,j) = this->_Z[i][j];

  gram = __gram(B);

  // Make U orthogonal to X+Z.
  for (size_t i = 0; i < this->_k; i++)
    for (size_t j = 0; j < n - 2*this->_k; j++) {
      // Clear components corresponding to X.
      Integer<T> scalar = gram(2*this->_k-1-i, 2*this->_k+j);
      scalar = (scalar / (p*p) + Integer<T>::one())*p*p-scalar;
      if (scalar >= p*p)
	scalar -= p*p;
      this->_U[j] += scalar.num() * this->_X[i];
      
      // Clear components corresponding to Z.
      scalar = gram(this->_k-1-i, 2*this->_k+j);
      scalar = (scalar / (p*p) + Integer<T>::one())*p*p-scalar;
      if (scalar >= p*p)
	scalar -= p*p;
      this->_U[j] += scalar.num() * this->_Z[i];
    }

#ifdef DEBUG
  for (size_t i = 0; i < n-2*this->_k; i++)
    for (size_t j = 0; j < n; j++)
      this->_U[i][j] = this->_U[i][j] % (p*p).num();

#ifdef DEBUG_LEVEL_FULL
  std::cerr << "after setting <U,X+Z> = 0" << std::endl;
  std::cerr << "U = " << this->_U << std::endl;
#endif // DEBUG_LEVEL_FULL
  
  // Verify that U is now orthogonal to X+Z.
  for (size_t i = 0; i < n-2*this->_k; i++)
    for (size_t j = 0; j < n; j++)
      B(2*this->_k+i,j) = this->_U[i][j];

  // The Gram matrix on this basis.
  temp = __gram(B);

  // Verify all is well.
  for (size_t i = 0; i < 2*this->_k; i++)
    for (size_t j = 2*this->_k; j < n; j++)
      assert(temp(i,j) % (p*p).num() == 0);

  // make sure that all the entries of U are between 0 and p^2
  for (size_t i = 0; i < n-2*this->_k; i++)
    for (size_t j = 0; j < n; j++)
      this->_U[i][j] = this->_U[i][j] % (p*p).num();
  
#endif // DEBUG

  return;
}

template<typename R, typename S, typename T, size_t n>
inline void NeighborManager<R,S,T,n>::nextIsotropicSubspace(void)
{
  if (this->_params.empty()) {
    // Move to the next pivot.
    this->_pivot_ptr++;
    
    // If we've exceeded the list of pivots, we're done.
    if (this->_pivot_ptr > this->_pivots.size()) {
      // Reset the pivot pointer so that we can loop through
      //  the isotropic subspaces again if needed.
      this->_pivot_ptr = 0;
      this->_iso_subspace.clear();
      return;
    }
    
    // Initialize the new pivot.
    this->__initializePivots();
  }

  // The list of evaluation values.
  FpElement<R,S> zero(this->_GF, 0);
  std::vector<FpElement<R,S> > eval_list(n*this->_k, zero);

  // Produce the isotropic subspace corresponding to the current
  //  parameters.
  for (size_t i = 0; i < this->_params.size(); i++)
    eval_list[this->_free_vars[i]] = this->_params[i];

  // The basis for the current isotropic subspace.
  this->_iso_subspace.clear();
  for (size_t i = 0; i < this->_k; i++) {
    VectorFp<R,S,n> vec_fp(this->_GF);
    for (size_t j = 0; j < n; j++) {
      vec_fp[j] = (*(this->_p_isotropic_param))(i,j).evaluate(eval_list);
    }
    this->_iso_subspace.push_back(vec_fp);
  }

  if (this->_free_vars.size() != 0) {
    // The current position in the parameterization.
    size_t pos = 0;
    FpElement<R,S> zero(this->_GF, 0);
    // Terminate loop once we found the next new subspace, or we
    //  hit the end of the list.
    do {
      // Increment position.
      pos++;
      // Manually move to the next element.
      assert(pos <= this->_params.size());
      this->_params[pos-1]++;
    } while ((pos != this->_free_vars.size()) && (this->_params[pos-1] == zero));
  }

  // If we've hit the end of the list, indicate we need to move on to the
  //  next pivot.

  bool all_zero = true;
  for (size_t i = 0; i < this->_params.size(); i++)
    if (this->_params[i] != 0) {
      all_zero = false;
      break;
    }

  if (all_zero) {
    this->_params.clear();
  }
  
  return;
}

// !! TODO - Is it enough to use greedy here as well?
template<typename R, typename S, typename T, size_t n>
inline GenusRep<T,n> NeighborManager<R,S,T,n>::getReducedNeighborRep(void)
{
  GenusRep<T,n> rep;

  rep.q = this->buildNeighbor(rep.s);
  // rep.q = QuadFormZZ<T,n>::reduce(rep.q, rep.s);
  SquareMatrixInt<T,n> qf = rep.q.bilinearForm();
  QuadFormZZ<T,n>::greedy(qf, rep.s);
  rep.q = qf;
  return rep;
}

// to representative of the line
// !! TODO  - do this only with actions of FpElement
template<typename R, typename S, typename T, size_t n>
inline VectorInt<R,n>
NeighborManager<R,S,T,n>::transformVector(const GenusRep<T,n>& dst,
					  VectorInt<R,n> src)
{
  VectorInt<T,n> temp;
  VectorFp<R,S,n> temp_mod = *mod(src, this->_GF);
  for (size_t i = 0; i < n; i++)
    temp[i] = temp_mod[i].lift();

  // should that be inverse mod p?
  Isometry<T,n> sinv = dst.s.inverse();
  temp = sinv * temp;

  for (size_t i = 0; i < n; i++)
    assert( temp[i] % sinv.getScale() == 0 );

  for (size_t i = 0; i < n; i++)
    temp[i] /= sinv.getScale();

  VectorInt<R,n> vec;
  for (size_t i = 0; i < n; i++)
    vec[i] = this->_GF->mod(temp[i]).lift();

  for (size_t i = n; i > 0; i--) {
    if (vec[i-1] != 0)
      {
	R inv = this->_GF->inverse(vec[i-1]);
	for (size_t j = 0; j < i-1; j++)
	  vec[j] = this->_GF->mod(this->_GF->mul(vec[j], inv)).lift();
	vec[i-1] = 1;
	break;
      }
  }

  return vec;
}

template<typename R, typename S, typename T, size_t n>
inline void NeighborManager<R,S,T,n>::_updateSkewMatrix(size_t & row, size_t & col)
{
  bool done;
  do {
    // Flag for determining whether we are done updating
    //  the skew matrix.
    done = true;
     
    // Increment value of the (row,col) position.
    (*(this->_p_skew))(row, col)++;
      
    // Update the coefficient of the skew matrix reflected
    //  across the anti-diagonal.
    (*(this->_p_skew))(this->_k-1-col, this->_k-1-row) = -(*(this->_p_skew))(row,col);
      
    // If we've rolled over, move on to the next position.
    if ((*(this->_p_skew))(row,col) == 0) {
      // The next column of our skew matrix.
      col++;
      // Are we at the end of the column?
      if (row+col == this->_k-1) {
	// Yes. Move to the next row.
	row++;
	// And reset the column.
	col = 0;
      }
      // Indicate we should repeat another iteration.
      done = false;
    }
  } while ((!done) && (row+col != this->_k-1));
  return;
}

template<typename R, typename S, typename T, size_t n>
inline void NeighborManager<R,S,T,n>::_updateSkewSpace(void)
{
  assert(this->_X.size() == this->_k);
  assert(this->_Z.size() == this->_k);
  
  Integer<T> p = this->_GF->prime();
  // Update the skew space.
  for (size_t i = 0; i < this->_k ; i++) {
    for (size_t j = 0; j < this->_k; j++){
      // !! TODO - I got rid here of X_skew,
      // check that it didn't destroy anything
      Integer<T> val = (*(this->_p_skew))(i,j).lift();
      this->_X[i] += p.num() * (val.num() * this->_Z[j]);
    }
  }
}

template<typename R, typename S, typename T, size_t n>
inline void NeighborManager<R,S,T,n>::getNextNeighbor(void)
{
  size_t row,col;
  // The starting position of the skew vector to update.
  row = 0;
  col = 0;
  // Update the skew matrix (only for k >= 2).
  if (this->_skew_dim != 0) {
    this->_updateSkewMatrix(row, col);
  }
  
  // If we haven't rolled over, update the skew space and return...
  if (row+col < this->_k-1) {
    this->_updateSkewSpace();
    return;
  }

  // ...otherwise, get the next isotropic subspace modulo p.
  this->nextIsotropicSubspace();

  // Lift the subspace if we haven't reached the end of the list.
  if (!(this->_iso_subspace.empty())) {
    this->_liftSubspace();
  }
  else {
    this->_X.clear();
    this->_Z.clear();
    this->_U.clear();
  }
    
  return;
}

template<typename R, typename S, typename T, size_t n>
inline QuadFormZZ<T,n> NeighborManager<R,S,T,n>::buildNeighbor(Isometry<T,n>& s) const
{
  T p = this->_GF->prime();
  T p2 = p*p;
  T p3 = p2*p;
  std::shared_ptr<const IntegerRing<T> > ZZ = std::make_shared< const IntegerRing<T> >(); 
  SquareMatrixInt<T,n> qq;

  // fill the isomtery by the X,Z,U
  // if lift_subspace was successful,
  // <X,X>, <X,U>,<Z,Z>,<Z,U> in p^2 and <X,Z> = 1 mod p^2

  // For now, we follow thw magma implementation
  // start by scaling the basis
  
  for (size_t i = 0; i < this->_k; i++)
    for (size_t j = 0; j < n; j++)
      s(i,j) = this->_X[i][j];

  for (size_t i = 0; i < this->_k; i++)
    for (size_t j = 0; j < n; j++)
      s(this->_k+i,j) = p2*this->_Z[i][j];

  for (size_t i = 0; i < n-2*this->_k; i++)
    for (size_t j = 0; j < n; j++)
      s(2*this->_k+i,j) = p*this->_U[i][j];

  SquareMatrixInt<T,n> hermite = s.hermiteForm(p3);
  
  s = hermite.transpose();
  s.setScale(p);
	 
  //	s.swap_cols(0, pivot);
  // need to adjust determinant for s to be in SO
  // This transforms using the isometry (and rescales already)
  qq = s.transform(this->_q.bilinearForm());

  QuadFormZZ<T,n> retval(qq);
	
  if (std::is_same<T,Z64>::value)
    {
      // If we're not using arbitrary precision, throw an exception if
      // the discriminant of the p-neighbor isn't correct.
      if (retval.discriminant() != this->_disc)
	{
	  throw std::overflow_error(
				    "An overflow has occurred. The p-neighbor's discriminant "
				    "does not match the original.");
	}
    }
  return retval;
}

template<typename R, typename S, typename T, size_t n>
inline VectorInt<R,n> NeighborManager<R,S,T,n>::_isotropicVector_p2(R t) const
{
  VectorInt<R,n> res;
  res[n-1] = 1;
    
  // Stub
  // !! TODO - do something appropriate here
       
  return res;
}

// A helper function for computing valid pivots.
template<typename R, typename S, typename T, size_t n>
inline std::vector< std::vector<size_t> >
NeighborManager<R,S,T,n>::__pivots(size_t dim, size_t aniso, size_t k)
{
  std::vector< std::vector<size_t> > pivs;
  // Base case.
  if (k == 1) {
    for (size_t i = 0; i < dim - aniso; i++) {
      std::vector<size_t> singleton(1);
      singleton[0] = i;
      pivs.push_back(singleton);
    }
    return pivs;
  }

  // Retrieve lower-dimensional maximal pivots.
  pivs = NeighborManager<R,S,T,n>::__pivots(dim-2, aniso, k-1);
  for (size_t i = 0; i < pivs.size(); i++)
    for (size_t j = 0; j < pivs[i].size(); j++)
      pivs[i][j]++;

  size_t num = pivs.size();
  // Determine the first set of pivots.
  pivs.insert(pivs.end(), pivs.begin(), pivs.end());
  for (size_t i = 0; i < num; i++){
    pivs[i].insert(pivs[i].begin(), 0);
    pivs[i+num].push_back(dim-aniso-1);
  } 

  // Add additional pivots when we're not in the maximal case.
  if (2*k <= dim - aniso) {
    std::vector< std::vector<size_t> > pivs2 = NeighborManager<R,S,T,n>::__pivots(dim-2, aniso, k);
    for (size_t i = 0; i < pivs2.size(); i++)
      for (size_t j = 0; j < pivs2[i].size(); j++)
	pivs2[i][j]++;
    pivs.insert(pivs.begin(), pivs2.begin(), pivs2.end());
  }
  return pivs;
}

template<typename R, typename S, typename T, size_t n>
inline void NeighborManager<R,S,T,n>::__initializePivots(void)
{
  assert(this->_pivot_ptr >= 1);

  std::vector<size_t> pivot = this->_pivots[this->_pivot_ptr-1];
  size_t rank = (this->_k)*n;
  
  // Keep a list of non-free variables from which we will determine the
  //  free variables when we are done.
  std::vector<size_t> remove;

  // Initialize matrix that will determine parameterization.
  std::vector< PolynomialFp<R, S> > data;
  for (size_t i = 0; i < rank; i++) {
    PolynomialFp<R,S> x_i = PolynomialFp<R,S>::x(this->_GF, i);
    data.push_back(x_i);
  }
  this->_p_isotropic_param =
    std::make_shared< Matrix< PolynomialFp<R, S>, PolynomialRingFp<R,S> > >(data, this->_k, n);

  FpElement<R,S> zero(this->_GF, 0);
  FpElement<R,S> one(this->_GF, 1);
  // Setup the columns corresponding to the pivots.
  for (size_t row = 0; row < this->_k; row++)
    for (size_t col = 0; col < this->_k; col++) {
      (*this->_p_isotropic_param)(row, pivot[col]) = (row == col) ? one : zero;
      remove.push_back(row*n + pivot[col]);
    }

  // Clear the rows prior to the pivot positions (but not the radical).
  for (size_t row = 0; row < this->_k; row++)
    for (size_t col = 0; col < pivot[row]; col++) {
      (*this->_p_isotropic_param)(row, col) = zero;
      remove.push_back(row*n + col);
    }

  // Check if one or more of the anisotropic coordinates need to be zero.
  for (size_t row = 0; row < this->_k; row++) {
    if (pivot[row] >= this->_witt_index) {
      for (size_t col = 0; col < this->_aniso_dim; col++) {
	(*this->_p_isotropic_param)(row, n-1-this->_rad_dim-col) = zero;
	remove.push_back((row+1)*n-1-this->_rad_dim-col);
      }
    }
  }

  // Determine the number of rows of the matrix that we'll echelonize.
  size_t rows = this->_k*(this->_k+1)/2;

  // Here we will save the quadratic polynomials
  data.clear();
  
  
  // The matrix that we're going to echelonize.
  MatrixFp<R,S> mat(this->_GF, rows, rank);
  
  // The current row to fill in in the matrix.
  size_t row = 0;
  for (size_t i = 0; i < this->_k; i++)
    for (size_t j = i; j < this->_k; j++) {
      // The appropriate vector that we want to be isotropic.
      std::vector<PolynomialFp<R, S> > vec;
      for (size_t r = 0; r < n; r++)
	vec.push_back((i == j) ? (*this->_p_isotropic_param)(i,r) :
		      (*this->_p_isotropic_param)(i,r) + (*this->_p_isotropic_param)(j,r));
   
      PolynomialFp<R,S> f = this->_p_q_std->evaluate(vec);
      // Degree 2 terms are inhomogenous.
      data.push_back(-f.quadraticPart());
      // The other terms are linear
      // so fill in the matrix accordingly.
      std::vector< FpElement<R,S> > l = f.linearPart(rank);
      for (size_t i = 0; i < rank; i++)
	mat(row, i) = l[i];
      // Move on to the next row.
      row++;
    }

#ifdef DEBUG //_LEVEL_FULL
  std::cerr << "The matrix before echelon is mat = " << mat << std::endl;
  std::cerr << "The last entry is the quadratic data = " << data << std::endl;
#endif
  
  // Compute the Echelon form of the matrix.
  MatrixFp<R,S> trans(this->_GF, rows, rows);
  MatrixFp<R,S>::rowEchelon(mat, trans);

  // The evaluation list for replacing variables with their dependence
  //  relations.
  std::vector< PolynomialFp<R,S> > eval_list;
  for (size_t i = 0; i < rank; i++) {
    PolynomialFp<R,S> x_i = PolynomialFp<R,S>::x(this->_GF, i);
    eval_list.push_back(x_i);
  }

  for (size_t i = 0; i < rows; i++) {
    // Find the pivot in the i-th row.
    size_t c = 0;
    while ((c < rank) && (mat(i,c) != 1)) c++;
    // Add this pivot to the list of non-free variables.
    remove.push_back(c);

    // If the row is entirely zero, skip it.
    if (c >= rank) continue;

    // Build the multinomial for which x_c is dependent.
    //    PolynomialFp<R,S> f = mat(i, rank);
    PolynomialFp<R,S> f(this->_GF);
    for (size_t j = 0; j < rows; j++) {
      f += trans(i,j) * data[j];
    }
    for (size_t j = 0; j < rank; j++) {
      if (j != c) {
	PolynomialFp<R,S> x_j = PolynomialFp<R,S>::x(this->_GF, j);
	f -= mat(i,j) * x_j;
      }
    }
    eval_list[c] = f;
  }

  // The matrix whose rows parameterize all isotropic subspaces.
  for (size_t row = 0; row < this->_p_isotropic_param->nrows(); row++)
    for (size_t col = 0; col < this->_p_isotropic_param->ncols(); col++)
      (*this->_p_isotropic_param)(row,col) =
	(*this->_p_isotropic_param)(row,col).evaluate(eval_list);

#ifdef DEBUG
  // Verify that we didn't screw up somewhere along the line.
#ifdef DEBUG // _LEVEL_FULL
  std::cerr << "testing parametrization" << std::endl;
#endif // DEBUG_LEVEL_FULL
  for (size_t i = 0; i < this->_k; i++)
    for (size_t j = 0; j < this->_k; j++) {
      std::vector<PolynomialFp<R,S> > vec;
      for (size_t r = 0; r < n; r++)
	vec.push_back((i == j) ? (*this->_p_isotropic_param)(i,r) :
		      (*this->_p_isotropic_param)(i,r) + (*this->_p_isotropic_param)(j,r));
#ifdef DEBUG // _LEVEL_FULL
      std::cerr << "Substituting vec = " << vec << std::endl;
      std::cerr << " in q_std = " << (*this->_p_q_std) << std::endl;
#endif // DEBUG_LEVEL_FULL
      PolynomialFp<R,S> f = this->_p_q_std->evaluate(vec);
#ifdef DEBUG //_LEVEL_FULL
      std::cerr << " yields f = " << f << std::endl;
#endif // DEBUG_LEVEL_FULL
      assert(f == zero);
    }
#endif // DEBUG

  // Determine the free variables.

  for (size_t i = 0; i < rank; i++) {
    bool appears = false;
    for (size_t row = 0; row < this->_p_isotropic_param->nrows(); row++)
      for (size_t col = 0; col < this->_p_isotropic_param->ncols(); col++)
	if ((*this->_p_isotropic_param)(row,col).degree(i) > 0) {
	  appears = true;
	  break;
	}
    if (!appears)
      remove.push_back(i);
  }

  this->_free_vars.clear();
  for (size_t i = 0; i < rank; i++) {
    std::vector<size_t>::const_iterator it;
    it = std::find(remove.begin(), remove.end(), i);
    if (it == remove.end()) {
      this->_free_vars.push_back(i);
      this->_params.push_back(zero);
    }
  }
  return;
}

template<typename R, typename S, typename T, size_t n>
std::vector< VectorFp<R,S,n> > NeighborManager<R,S,T,n>::radical(void) const
{
  SquareMatrixFp<R,S,n> basis = (*this->_p_basis).transpose();
  std::vector< VectorFp<R,S,n> > rad;
  for (size_t i = 0; i < this->_rad_dim; i++)
    rad.push_back(basis[n-1-i]);
  return rad;
}
