#include <cassert>
#include <random>
#include <unordered_set>

// implementation file for header Genus.h

template<typename R, size_t dim>
inline std::set< Integer<R> >
Genus<R,dim>::_wittToHasse(const Integer<R>& det,
			   const std::set<std::pair<Integer<R>,int> > & finite)
{
  std::set< Integer<R> > hasse;
  int c_table[8] = {2, 1, 1, -2, -2, -1, -1, 2};
  int c_mask = c_table[dim % 8];
  Integer<R> c = R((c_mask / 2)*det + c_mask % 2);
  Integer<R> minus_one = -Integer<R>::one();

  for (std::pair<Integer<R>,int> x : finite)
    if (x.second != minus_one.hilbertSymbol(c, x.first))
      hasse.insert(x.first);
  
  return hasse;
}

template<typename R, size_t n>
inline Rational<Z> Genus<R,n>::_localFactor(const MatrixRat<R> & g,
					    const Integer<R> & p)
{
  size_t m = g.ncols();
  Rational<Z> one = Rational<Z>::one();
  Integer<Z> p_sqr = birch_util::convert_Integer<R,Z>((p*p).num());
  Rational<Z> f = one;
  Rational<Z> p_i = one/p_sqr;
  for (size_t i = 2; i+2 <= m; i+= 2)
    {
      f *= (one - p_i);
      p_i /= p_sqr;
    }
  if (m % 2 == 1) {
    if (m != 1) f *= (one - p_i);
    return f;
  }
  size_t r = m / 2;
  R sign = (r % 2 == 0) ? 1 : -1;
  Rational<R> d = g.determinant() * sign;
  
  if (((d.valuation(p)) % 2) == 0) {
    p_i = one;
    for (size_t i = 0; i < r; i++) p_i /= birch_util::convert_Integer<R,Z>(p.num());
    if (d.isLocalSquare(p))
      f *= one - p_i;
    else
      f *= one + p_i;
  }
  return f;
}

template<typename R, size_t n>
inline Rational<Z> Genus<R,n>::_combine(const QuadFormZZ<R,n>& q,
					const Integer<R> & p)
{
  assert(p.num() != 2);
  typename QuadFormZZ<R,n>::jordan_data jordan = q.jordanDecomposition(p);
  Integer<Z> one = Integer<Z>::one();
  Rational<Z> f = one;
  Rational<Z64> e = 0;
  std::vector<size_t> ms;
  size_t m = 0;
  for (MatrixRat<R> g : jordan.grams) {
    ms.push_back(g.ncols());
    m += g.ncols();
  }
  for (size_t i = 0; i < ms.size(); i++) {
    Z64 t = (i == 0) ? 0 : jordan.exponents[i-1];
    Rational<Z64> tmp1((jordan.exponents[i]-t)*(m+1)*m,2);
    Rational<Z64> tmp2(jordan.exponents[i]*(n+1)*ms[i],2);
    e += tmp1-tmp2;
    f *= Genus<R,n>::_localFactor(jordan.grams[i], p);
    m -= ms[i];
  }
  // !! We might run into trouble at 2 here
  // check if we need disc or half-disc
  // size_t v = Math<R>::valuation(q.discriminant(), p);
  MatrixInt<R> q_mat = q.bilinearForm();
  size_t v = q_mat.determinant().valuation(p);
  if ((n % 2 == 0) && (v % 2 == 1)) {
    Rational<Z64> n_rat = (Z64)n;
    e += (n_rat-1)/2;
  }
  assert(e.isIntegral());
  Integer<Z> p_Z = birch_util::convert_Integer<R, Z>(p.num());
  Rational<Z> p_e = p_Z^(e.floor());

  Z pow2 = 1 << (jordan.grams.size()-1);
  Rational<Z> denom = pow2 * f * p_e;
  MatrixRat<R> diag = MatrixRat<R>::diagonalJoin(jordan.grams);
  return Genus<R,n>::_localFactor(diag, p) / denom;
}

template<typename R, size_t n>
Rational<Z> Genus<R,n>::_getMass(const QuadFormZZ<R,n>& q,
				 const std::vector<PrimeSymbol<R>>& symbols)
{
  size_t r = n / 2;

  std::set<std::pair<Integer<R>, int> > hasse;
  // do we need the dummy? We could probably let it go
  size_t dummy;
  Integer<R> det = q.invariants(hasse, dummy);
  std::set<Integer<R> > witt = Genus<R,n>::_wittToHasse(det, hasse);
  std::vector< std::pair<Integer<R>, size_t> > fac = det.factorization();
  std::set<Integer<R> > B;
  // TODO - replace these by set_union, set_difference ?
  for (std::pair<Integer<R>, size_t> fa : fac)
    B.insert(fa.first);
  for (Integer<R> p : witt)
    B.insert(p);
  B.erase(2);
  for (Integer<R> p : B)
    witt.erase(p);
  size_t val2 = det.valuation(R(2));
     
  // mass from infinity and 2
  Rational<Z> mass(1, 1<<r);    
     
  for (size_t i = 1; i < n / 2 + n % 2; i++)
    mass *= -Rational<Z>::bernoulliNumber(2*i)/(2*i);
     
  if (n % 2 == 1)
    {	 
      if (val2 % 2 == 1)
	{
	  mass *= (1 << r) + ((witt.find(2) != witt.end()) ? -1 : 1);
	  mass /= Z(2);
	  witt.erase(2);
	}
      if (witt.find(2) != witt.end())
	{ 
	  mass *= (1 << (n-1)) - 1;
	  mass /= Z(6);
	}
    }
  else
    {
      Integer<R> disc = (r % 2 == 1) ? -det : det;
      if (disc.isSquare())
	mass *= -Rational<Z>::bernoulliNumber(r)/r;
      else
	{
	  Integer<Z> disc_z = birch_util::convert_Integer<R,Z>(disc);
	  mass *= -Rational<Z>::bernoulliNumber(r, disc_z) / r;
	  if (r % 2 == 0)
	    mass *= -Rational<Z>::bernoulliNumber(r) / r;
	  if (val2 % 2 == 1)
	    {
	      mass /= Z(2);
	      witt.erase(2);
	    }
	}
      if (witt.find(2) != witt.end())
	{
	  // checking if disc is a local square at 2
	  int w = 1;
	  if ((val2 % 2 != 1) && ((disc >> val2) % 8 == 1))
	    w = -1;
	  mass *= (1<<(r-1))+w;
	  mass *= (1<<r)+w;
	  mass /= 6;
	}
    }
  // odd places which are not unimodular or have Witt invariant -1.
  for (Integer<R> p : B)
    mass *= Genus<R,n>::_combine(q,p);
     
  return mass.abs();
}

template<typename R, size_t n>
Genus<R,n>::Genus(const QuadFormZZ<R,n>& q,
		  const std::vector<PrimeSymbol<R>>& symbols, W64 seed)
{
  if (seed == 0)
    {
      std::random_device rd;
      seed = rd();
    }
  
  this->disc = q.discriminant();
  this->_seed = seed;
  
  this->_prime_divisors.reserve(symbols.size());
  for (const PrimeSymbol<R>& symb : symbols)
    {
      this->_prime_divisors.push_back(symb.p);
    }
  
  Spinor<R> *spin = new Spinor<R>(this->_prime_divisors);
  this->spinor = std::unique_ptr<Spinor<R>>(spin);
  
  if (symbols.size() > 63)
    {
      throw std::domain_error("Must have 63 or fewer prime divisors.");
    }
  
  size_t num_conductors = 1LL << symbols.size();
  
  this->_conductors.reserve(num_conductors);
  this->_conductors.push_back(1);
  
  size_t bits = 0;
  size_t mask = 1;
  for (size_t c=1; c<num_conductors; c++)
    {
      if (c == 2*mask)
	{
	  ++bits;
	  mask = 1LL << bits;
	}
      R value = this->_prime_divisors[bits] * this->_conductors[c ^ mask];
      this->_conductors.push_back(value);
    }

  // Set the mass. This value is used to determine
  // when the genus has been fully populated.
  this->_mass = this->_getMass(q, symbols);

  // The mass provides a reasonable estimate for the size of the genus.
  Z64 estimated_size = mpz_get_si(this->_mass.ceiling().get_mpz_t());

  // Should this be 1/#aut or 2/#aut? probably depends if this is SO or O
  GenusRep<R,n> rep;
  Isometry<R,n> s;
  rep.q = QuadFormZZ<R,n>::reduce(q,s);
  Z num_aut = rep.q.numAutomorphisms();
  Rational<Z> sum_mass(1, num_aut);
  
  rep.p = 1;
  rep.parent = -1;

  auto *ptr = new HashMap<GenusRep<R,n>>(estimated_size);
  this->_hash = std::unique_ptr<HashMap<GenusRep<R,n>>>(ptr);
  this->_hash->add(rep);

  // A temporary placeholder for the genus representatives before they
  // are fully built.
  GenusRep<R,n> foo;
  foo.p = 1;
  foo.parent = -1;
  
  auto *inv_ptr = new HashMap<GenusRep<R,n>>(estimated_size);
  this->_inv_hash = std::unique_ptr<HashMap<GenusRep<R,n>>>(inv_ptr);

  // add the orbit representatives to the invariants
  std::unordered_map< QuadFormZZ<R,n>, Isometry<R,n> > q_orbit;
  q_orbit = rep.q.generateOrbit();
  typename std::unordered_map<QuadFormZZ<R,n>, Isometry<R,n> >::const_iterator
    iter;

  for (iter = q_orbit.begin(); iter != q_orbit.end(); iter++) {
    foo.q = iter->first;
    foo.s = iter->second;
    this->_inv_hash->add(foo);
    this->_inv_map[this->_inv_hash->indexof(foo)] =
      this->_hash->indexof(rep);
  }
  
  // The spinor primes hash table, used to identify the primes used in
  // constructing the genus representatives.
  auto *ptr2 = new HashMap<W16>();
  this->_spinor_primes = std::unique_ptr<HashMap<W16>>(ptr2);
  
  Z p = 1;
  W16 prime = 1;

  bool done = (sum_mass == this->_mass);
  while (!done)
    {
      // !! TODO - we don't really need to restrict to good primes here,
      // but let's check these first
      // Get the next good prime and build the appropriate finite field.
      do
	{
	  mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
	  prime = mpz_get_ui(p.get_mpz_t());
	}
      while (this->_disc % prime == 0);
      std::shared_ptr<W16_Fp> GF;
      if (prime == 2)
	GF = std::make_shared<W16_F2>(prime, this->_seed);
      else
	GF = std::make_shared<W16_Fp>(prime, this->_seed, true);
      
      size_t current = 0;
      while (!done && current < this->_hash->size())
	{
	  // Get the current quadratic form and build the neighbor manager.
	  const QuadFormZZ<R,n>& mother = this->_hash->get(current)._q;
	  NeighborManager<W16,W32,R,n> manager(mother, GF);

#ifdef DEBUG
	  // Build the affine quadratic form for debugging purposes.
	  std::shared_ptr< W16_QuadForm<n> > qp = mother.mod(GF);
#endif
	  manager.getNextNeighbor();
	  bool prime_done = manager.getIsotropicSubspace().empty();
	  while ((!done) && (!prime_done))
	    {     
	      // Verify that the appropriate vector is isotropic.
	      assert(!manager.getIsotropicSubspace().empty());
	      assert( mother.evaluate(manager.getIsotropicSubspace()[0]) % prime == 0 );
	      
	      // Construct the neighbor, the isometry is stored in s.
	      foo.s.setIdentity();
	      foo.q = manager.buildNeighbor(foo.s);
	      
	      // Verify neighbor discriminant matches.
	      assert( foo.q.discriminant() == mother.discriminant() );
	      // Verify the isometry indeed constructs the neighbor
	      assert( foo.s.transform(mother.bilinearForm()) ==
		      foo.q.bilinearForm() );

	      // Reduce the neighbor to its Eisenstein form and add it to
	      // the hash table.

	      // Here we also want to compute the automorphism group
	      // !!TODO - ?? Do we want this ??
	      // We can compute it only if we need to add it.
	      
	      foo.q = QuadFormZZ<R,n>::reduce(foo.q, foo.s, true);

	      assert( foo.s.transform(mother.bilinearForm()) ==
		      foo.q.bilinearForm() );

	      foo.p = prime;
	      foo.parent = current;

	      bool added = this->_hash->add(foo);
	      if (added)
		{
		  const GenusRep<R,n>& temp = this->_hash->last();
		  Z num_aut = temp.q.numAutomorphisms();
		  Rational<Z> q_mass(1, num_aut);
		  sum_mass += q_mass;
		  done = (sum_mass == this->mass);

		  assert(sum_mass <= this->mass);
		  
		  this->_spinor_primes->add(prime);

		  // add the orbit representatives to the invariants
		  q_orbit = temp.q.generateOrbit();		  

		  assert(q_orbit.find(temp.q) != q_orbit.end());

		  for (iter = q_orbit.begin(); iter != q_orbit.end(); iter++) {

		    assert(iter->second.isIsometry(temp.q, iter->first));

		    foo.q = iter->first;
		    foo.s = temp.s * iter->second;
		    foo.parent = temp.parent;

		    assert(foo.s.isIsometry(this->_hash->at(temp.parent).q,
					     foo.q));

		    this->_inv_hash->add(foo);
		    this->_inv_map[this->_inv_hash->indexof(foo)] =
		      this->_hash->indexof(temp);
		  }
		}
	      manager.getNextNeighbor();
	      prime_done = manager.getIsotropicSubspace().empty();
	    }
	  ++current;
	}
    }
  
  // Initialize the dimensions to zero, we will compute these values below.
  this->_dims.resize(num_conductors, 0);
  
  // Create the lookup table values for each genus rep at each conductor.
  size_t genus_size = this->_hash->size();
  this->_lut_positions.resize(num_conductors, std::vector<int>(genus_size, -1));
  this->_num_auts.resize(num_conductors);

  assert(this->_hash->size() > 0);
  
#ifdef DEBUG
  GenusRep<R,n>& mother = this->_hash->at(0);
#endif
  
  // The genus rep isometries were initialized only to contain the
  // isometry between the parent and its child, we now want to update
  // these isometries so that they are rational isometries between the
  // "mother" quadratic form and the genus rep.
  for (size_t idx=0; idx<this->_hash->size(); idx++)
    {
      GenusRep<R,n>& rep = this->_hash->at(idx);
      
      // Only compute composite isometries if we are not considering the
      // mother form.
      if (idx)
	{
	  GenusRep<R,n>& parent = this->_hash->at(rep.parent);
	  assert( rep.s.transform(parent.q.bilinearForm()) ==
		  rep.q.bilinearForm() );
	  assert( parent.s.transform(mother.q.bilinearForm()) ==
		  parent.q.bilinearForm() );  

	  // Construct the isometries to/from the mother quadratic form.
	  rep.sinv = rep.s.inverse();
	  rep.sinv = rep.sinv * parent.sinv;
	  rep.s = parent.s * rep.s;

	  // Copy the numerators, and increment the genus rep prime.
	  rep.es = parent.es;
	  ++rep.es[rep.p];

	  // Verify that s is an isometry from the mother form to the rep,
	  // and that sinv is an isometry from the rep to the mother form.
	  assert( rep.s.transform(mother.q.bilinearForm()) ==
		  rep.q.bilinearForm() );
	  assert( rep.s.isIsometry(mother.q, rep.q) );
	  assert( rep.sinv.isIsometry(rep.q, mother.q) );

	}
      
      // Determine which subspaces this representative contributes.
      std::set<Isometry<R,n>> auts = rep.q.properAutomorphisms();

      std::vector<bool> ignore(this->_conductors.size(), false);
      for (const Isometry<R,n>& s : auts)
	{
	  Z64 vals = this->_spinor->norm(rep.q, s, 1);
	  
	  for (size_t k=0; k<num_conductors; k++)
	    {
	      if (!ignore[k] && (birch_util::popcnt(vals & k) & 1))
		{
		  ignore[k] = true;
		}
	    }
	}
      
      int num = rep.q.numAutomorphisms();
      for (size_t k=0; k<num_conductors; k++)
	{
	  if (!ignore[k])
	    {
	      this->_lut_positions[k][idx] = this->-dims[k];
	      this->_num_auts[k].push_back(num);
	    }
	  this->_dims[k] += (ignore[k] ? 0 : 1);
	}
    }

  // Do the same for the inv_hash
  for (size_t idx=0; idx<this->_inv_hash->size(); idx++)
    {
      GenusRep<R,n>& rep = this->_inv_hash->at(idx);
      
      // Only compute composite isometries if we are not considering the
      // mother form.
      if (rep.parent != -1)
	{
	  GenusRep<R,n>& parent = this->_hash->at(rep.parent);

	  assert( rep.s.transform(parent.q.bilinearForm()) ==
		  rep.q.bilinearForm() );
	  assert( parent.s.transform(mother.q.bilinearForm()) ==
		  parent.q.bilinearForm() );  

	  // Construct the isometries to/from the mother quadratic form.
	  rep.sinv = rep.s.inverse();
	  rep.sinv = rep.sinv * parent.sinv;
	  rep.s = parent.s * rep.s;

	  // Copy the numerators, and increment the genus rep prime.
	  rep.es = parent.es;
	  ++rep.es[rep.p];

	  // Verify that s is an isometry from the mother form to the rep,
	  // and that sinv is an isometry from the rep to the mother form.
	  assert( rep.s.transform(mother.q.bilinearForm()) ==
		  rep.q.bilinearForm() );
	  assert( rep.s.is_isometry(mother.q, rep.q) );
	  assert( rep.sinv.is_isometry(rep.q, mother.q) );

	}
      
    }
}

template<typename R, size_t n>
template<typename T>
Genus<R,n>::Genus(const Genus<T,n>& src)
{
  // Convert the discriminant.
  this->_disc = birch_util::convert_Integer<T,R>(src._disc);

  // Convert the prime divisors.
  for (const T& p : src._prime_divisors)
    {
      this->_prime_divisors.push_back(birch_util::convert_Integer<T,R>(p));
    }

  // Convert the conductors.
  for (const T& cond : src._conductors)
    {
      this->_conductors.push_back(birch_util::convert_Integer<T,R>(cond));
    }

  // Copy dimensions.
  this->_dims = src._dims;

  // Copy automorphisms counts.
  this->_num_auts = src._num_auts;

  // Copy lookup table dimensions.
  this->_lut_positions = src._lut_positions;

  // Copy mass.
  this->_mass = src._mass;

  // Build a copy of the spinor primes hash table.
  this->_spinor_primes = std::unique_ptr<HashMap<W16>>(new HashMap<W16>(src._spinor_primes->size()));
  for (W16 x : src._spinor_primes->keys())
    {
      this->_spinor_primes->add(x);
    }
  
  // Build a copy of the genus representatives hash table.
  this->_hash = std::unique_ptr<HashMap<GenusRep<R,n>>>(new HashMap<GenusRep<R,n>>(src._hash->size()));
  
  for (const GenusRep<T, n>& rep : src._hash->keys())
    {
      this->_hash->add(birch_util::convert_GenusRep<T,R>(rep));
    }

  this->_inv_hash = std::unique_ptr<HashMap<GenusRep<R,n>>>(new HashMap<GenusRep<R,n>>(src._inv_hash->size()));
  
  for (const GenusRep<T, n>& rep : src._inv_hash->keys())
    {
      this->_inv_hash->add(birch_util::convert_GenusRep<T,R>(rep));
    }

  for (std::pair<size_t,size_t> element : src._inv_map) {
    this->_inv_map[element.first] = element.second;
  }

  // Create Spinor class.
  std::vector<R> primes;
  primes.reserve(src._spinor->primes().size());
  for (const T& p : src._spinor->primes())
    {
      primes.push_back(birch_util::convert_Integer<T,R>(p));
    }
  this->_spinor = std::unique_ptr<Spinor<R>>(new Spinor<R>(primes));
  
  // Copy seed.
  this->_seed = src._seed;
}

template<typename R, size_t n>
Eigenvector<R> Genus<R,n>::eigenvector(const std::vector<Z32>& vec,
				       const R& conductor) const
{
  size_t num_conductors = this->_conductors.size();
  bool found = false;

  size_t k;
  for (k=0; k<num_conductors; k++)
    {
      if (this->_conductors[k] == conductor)
	{
	  found = true;
	  break;
	}
    }

  if (!found)
    {
      throw std::invalid_argument("Invalid conductor.");
    }

  size_t dim = this->_dims[k];
  if (dim != vec.size())
    {
      throw std::invalid_argument("Eigenvector has incorrect dimension.");
    }
  
  size_t fulldim = this->size();

  std::vector<Z32> temp(this->size());
  const std::vector<int>& lut = this->_lut_positions[k];
  
  for (size_t idx=0; idx<fulldim; idx++)
    {
      if (lut[idx] != -1)
	{
	  temp[idx] = vec[lut[idx]];
	}
    }

  return Eigenvector<R>(std::move(temp), k);
}

template<typename R, size_t n>
std::vector<Z32>
Genus<R,n>::eigenvalues(EigenvectorManager<R,n>& vector_manager,
			const R& p) const
{
  R bits16 = birch_util::convert_Integer<Z64,R>(1LL << 16);
  R bits32 = birch_util::convert_Integer<Z64,R>(1LL << 32);

  if (p == 2)
    {
      W16 prime = 2;
      std::shared_ptr<W16_F2> GF = std::make_shared<W16_F2>(prime, this->seed());
      return this->_eigenvectors<W16,W32>(vector_manager, GF, p);
    }
  else if (p < bits16)
    {
      W16 prime = birch_util::convert_Integer<R,W16>(p);
      std::shared_ptr<W16_Fp> GF = std::make_shared<W16_Fp>(prime, this->seed(), true);
      return this->_eigenvectors<W16,W32>(vector_manager, GF, p);
    }
  else if (p < bits32)
    {
      W32 prime = birch_util::convert_Integer<R,W32>(p);
      std::shared_ptr<W32_Fp> GF = std::make_shared<W32_Fp>(prime, this->seed(), false);
      return this->_eigenvectors<W32,W64>(vector_manager, GF, p);
    }
  else
    {
      W64 prime = birch_util::convert_Integer<R,W64>(p);
      std::shared_ptr<W64_Fp> GF = std::make_shared<W64_Fp>(prime, this->seed(), false);
      return this->_eigenvectors<W64,W128>(vector_manager, GF, p);
    }
}

template<typename R, size_t n>
template<typename S, typename T>
std::vector<Z32>
Genus<R,n>::_eigenvectors(EigenvectorManager<R,n>& vector_manager,
			  std::shared_ptr<Fp<S,T>> GF, const R& p) const
{
  std::vector<Z32> eigenvalues(vector_manager.size());

  S prime = GF->prime();

  const GenusRep<R,n>& mother = this->_hash->get(0);

  const Z32 *stride_ptr = vector_manager._strided_eigenvectors.data();

  size_t num_indices = vector_manager._indices.size();
  for (size_t index=0; index<num_indices; index++)
    {
      size_t npos = static_cast<size_t>(vector_manager._indices[index]);
      const GenusRep<R,n>& cur = this->_hash->get(npos);
      NeighborManager<S,T,R,n> neighbor_manager(cur.q, GF);
      neighbor_manager.getNextNeighbor();
      bool done = neighbor_manager.getIsotropicSubspace().empty();
      //for (W64 t=0; t<=prime; t++)
      while (!done)
	{
	  GenusRep<R,n> foo = neighbor_manager.getReducedNeighborRep();
	  
	  size_t rpos = this->_hash->indexof(foo);
	  size_t offset = vector_manager._stride * rpos;
	  __builtin_prefetch(stride_ptr + offset, 0, 0);

	  W64 spin_vals;
	  if (unlikely(rpos == npos))
	    {
	      spin_vals = this->_spinor->norm(foo.q, foo.s, p);
	    }
	  else
	    {
	      const GenusRep<R,n>& rep = this->_hash->get(rpos);
	      foo.s = cur.s * foo.s;
	      R scalar = p;
	      
	      foo.s = foo.s * rep.sinv;

	      scalar *= birch_util::my_pow(cur.es);
	      scalar *= birch_util::my_pow(rep.es);
	      
	      spin_vals = this->_spinor->norm(mother.q, foo.s, scalar);
	    }
	  
	  for (Z64 vpos : vector_manager._position_lut[index])
	    {
	      W64 cond = vector_manager._conductors[vpos];
	      Z32 value = birch_util::char_val(spin_vals & cond);
	      Z32 coord = vector_manager._strided_eigenvectors[offset + vpos];
	      if (likely(coord))
		{
		  eigenvalues[vpos] += (value * coord);
		}
	    }
	  neighbor_manager.getNextNeighbor();
	  done = neighbor_manager.getIsotropicSubspace().empty();
	}
      
      // Divide out the coordinate associated to the eigenvector to
      // recover the actual eigenvalue.
      for (Z64 vpos : vector_manager._position_lut[index])
	{
	  size_t offset = vector_manager._stride * npos;
	  Z32 coord = vector_manager._strided_eigenvectors[offset + vpos];
	  assert( eigenvalues[vpos] % coord == 0 );
	  eigenvalues[vpos] /= coord;
	}
    }
  
  return eigenvalues;
}

template<typename R, size_t n>
inline std::map<R,std::vector<std::vector<int>>>
Genus<R,n>::_heckeMatrixSparseInternal(const R& p) const
{
  size_t num_conductors = this->_conductors.size();
  size_t num_primes = this->_prime_divisors.size();

  std::vector<std::vector<int>> data(num_conductors);
  std::vector<std::vector<int>> indptr;
  std::vector<std::vector<int>> indices(num_conductors);

  W16 prime = birch_util::convert_Integer<R,W16>(p);

  std::shared_ptr<W16_Fp> GF;
  if (prime == 2)
    GF = std::make_shared<W16_F2>(2, this->seed());
  else
    GF = std::make_shared<W16_Fp>((W16)prime, this->seed(), true);

  std::vector<W64> all_spin_vals;
  all_spin_vals.reserve(prime+1);

  std::vector<std::vector<int>> rowdata;
  for (int dim : this->_dims)
    {
      rowdata.push_back(std::vector<int>(dim));
      indptr.push_back(std::vector<int>(dim+1, 0));
    }

  const GenusRep<R,n>& mother = this->_hash->keys()[0];
  size_t num_reps = this->size();
  for (size_t idx=0; idx<num_reps; idx++)
    {
      const GenusRep<R,n>& cur = this->_hash->get(idx);
      NeighborManager<W16,W32,R,n> manager(cur.q, GF);

      manager.getNextNeighbor();
      bool done = manager.getIsotropicSubspace().empty();

      while (!done)
	{
	  GenusRep<R,n> foo = manager.getReducedNeighborRep();

	  assert( foo.s.isIsometry(cur.q, foo.q) );

	  size_t r = this->_hash->indexof(foo);

	  assert( r < this->size() );

	  W64 spin_vals;
	  if (r == idx)
	    {
	      spin_vals = this->_spinor->norm(foo.q, foo.s, p);
	    }
	  else
	    {
	      const GenusRep<R,n>& rep = this->_hash->get(r);
	      foo.s = cur.s * foo.s;
	      R scalar = p;

	      assert( foo.s.isIsometry(mother.q, foo.q) );

	      foo.s = foo.s * rep.sinv;

	      assert( foo.s.isIsometry(mother.q, mother.q) );

	      scalar *= birch_util::my_pow(cur.es);
	      scalar *= birch_util::my_pow(rep.es);

	      spin_vals = this->_spinor->norm(mother.q, foo.s, scalar);
	    }

	  all_spin_vals.push_back((r << num_primes) | spin_vals);
	  manager.getNextNeighbor();
	  done = manager.getIsotropicSubspace().empty();
	}

      for (size_t k=0; k<num_conductors; k++)
	{
	  const std::vector<int>& lut = this->_lut_positions[k];
	  int npos = lut[idx];
	  if (npos == -1) continue;

	  // Populate the row data.
	  std::vector<int>& row = rowdata[k];
	  for (W64 x : all_spin_vals)
	    {
	      int r = x >> num_primes;
	      int rpos = lut[r];
	      if (rpos == -1) continue;

	      int value = birch_util::char_val(x & k);
	      row[rpos] += value;
	    }

	  // Update data and indices with the nonzero values.
	  size_t nnz = 0;
	  size_t pos = 0;
	  std::vector<int>& data_k = data[k];
	  std::vector<int>& indices_k = indices[k];
	  for (int x : row)
	    {
	      if (x)
		{
		  data_k.push_back(x);
		  indices_k.push_back(pos);
		  row[pos] = 0; // Clear the nonzero entry.
		  ++nnz;
		}
	      ++pos;
	    }

	  // Update indptr
	  indptr[k][npos+1] = indptr[k][npos] + nnz;
	}

      all_spin_vals.clear();
    }

  std::map<R,std::vector<std::vector<int>>> csr_matrices;
  for (size_t k=0; k<num_conductors; k++)
    {
      const R& cond = this->_conductors[k];
      csr_matrices[cond] = std::vector<std::vector<int>>();
      csr_matrices[cond].push_back(data[k]);
      csr_matrices[cond].push_back(indices[k]);
      csr_matrices[cond].push_back(indptr[k]);
    }
  return csr_matrices;
}

template<typename R, size_t n>
inline std::map<R,std::vector<int>>
Genus<R,n>::_heckeMatrixDenseInternal(const R& p) const
{
  size_t num_conductors = this->_conductors.size();
  size_t num_primes = this->_prime_divisors.size();

  // Allocate memory for the Hecke matrices and create a vector to store
  // pointers to the raw matrix data.
  std::vector<int*> hecke_ptr;
  hecke_ptr.reserve(num_conductors);
  std::vector<std::vector<int>> hecke_matrices;
  for (size_t k=0; k<num_conductors; k++)
    {
      size_t dim = this->_dims[k];
      hecke_matrices.push_back(std::vector<int>(dim * dim));
      hecke_ptr.push_back(hecke_matrices.back().data());
    }

  W16 prime = birch_util::convert_Integer<R,W16>(p);
  std::vector<W64> all_spin_vals;
  all_spin_vals.reserve(prime+1);

  std::shared_ptr<W16_Fp> GF;
  if (prime == 2)
    GF = std::make_shared<W16_F2>(2, this->seed());
  else
    GF = std::make_shared<W16_Fp>((W16)prime, this->seed(), true);

  const GenusRep<R,n>& mother = this->_hash->keys()[0];
  size_t num_reps = this->size();

  // Create hash tables for storing isotropic vectors to be skipped
  // at later iterations.
  std::vector<HashMap<W16_Vector<n> >> vector_hash(num_reps);

  for (size_t idx=0; idx<num_reps; idx++)
    {
      const GenusRep<R,n>& cur = this->_hash->get(idx);
      NeighborManager<W16,W32,R,n> manager(cur.q, GF);

      manager.getNextNeighbor();
      bool done = manager.getIsotropicSubspace().empty();

      while (!done)
	{
	  GenusRep<R,n> foo;
	  
	  W16_Vector<n> vec;
	  for (size_t i = 0; i < n; i++)
	    vec[i] = GF->mod(manager.getIsotropicSubspace()[0][i]).lift();

	  // !! TODO - this stopped working when we started storing orbits - fix this. We skip this part for now
	  /*
	  // If this vector has already been identified, skip it. This
	  // avoids unnecessarily computing the neighbor, reducing it,
	  // and testing for isometry. The Hermitian symmetry property
	  // of the Hecke matrix will account for this once we finish
	  // processing neighbors.
	  
	  if (vector_hash[idx].exists(vec)) {
	    manager.get_next_neighbor();
	    done = manager.get_isotropic_subspace().empty();
	    continue;
	  }
	  */
	  // Build the neighbor and reduce it.
	  foo.q = manager.buildNeighbor(foo.s);
	  SquareMatrix<R, n> qf = foo.q.bilinearForm();
	  QuadFormZZ<R,n>::greedy(qf, foo.s);
	  QuadFormZZ<R,n> qred(qf);
	  foo.q = qred;

	  assert( foo.s.isIsometry(cur.q, foo.q) );

	  size_t r_inv = this->_inv_hash->indexof(foo);
	  size_t r = this->_inv_map.at(r_inv);

	  assert( r < this->size() );

	  W64 spin_vals;
	  if (r > idx)
	    {
	      const GenusRep<R,n>& rep = this->_inv_hash->get(r);
	      const GenusRep<R,n>& rep_inv = this->_inv_hash->get(r_inv);
	      
	      GenusRep<R,n> tmp = rep;

	      tmp.s = foo.s * rep_inv.sinv * rep.s;
	      
	      assert( tmp.s.isIsometry(cur.q, tmp.q) );
	      
	      W16_Vector<n> result = manager.transformVector(foo, vec);
	      
	      vector_hash[r].add(result);
	      
	      foo.s = cur.s * foo.s;
	      R scalar = p;

	      assert( foo.s.isIsometry(mother.q, foo.q) );

	      foo.s = foo.s * rep_inv.sinv;

	      assert( foo.s.isIsometry(mother.q, mother.q) );

	      scalar *= birch_util::my_pow(cur.es);
	      scalar *= birch_util::my_pow(rep_inv.es);

	      spin_vals = this->_spinor->norm(mother.q, foo.s, scalar);
	    }
	  else if (r == idx)
	    {
	      spin_vals = this->_spinor->norm(foo.q, foo.s, p);
	    }
	  else {
	    manager.getNextNeighbor();
	    done = manager.getIsotropicSubspace().empty();
	    continue;
	  }

	  all_spin_vals.push_back((r << num_primes) | spin_vals);
	  manager.getNextNeighbor();
	  done = manager.getIsotropicSubspace().empty();
	}

      for (size_t k=0; k<num_conductors; k++)
	{
	  const std::vector<int>& lut = this->_lut_positions[k];
	  int npos = lut[idx];
	  if (unlikely(npos == -1)) continue;

	  int *row = hecke_ptr[k];

	  for (W64 x : all_spin_vals)
	    {
	      int r = x >> num_primes;
	      int rpos = lut[r];
	      if (unlikely(rpos == -1)) continue;

	      row[rpos] += birch_util::char_val(x & k);
	    }

	  hecke_ptr[k] += this->_dims[k];
	}

      all_spin_vals.clear();
    }

  // Copy the upper diagonal entries to the lower diagonal using the
  // Hermitian symmetry property and then move the matrix into an
  // associatively map before returning.
  std::map<R,std::vector<int>> matrices;
  for (size_t k=0; k<num_conductors; k++)
    {
      std::vector<int>& matrix = hecke_matrices[k];
      size_t dim = this->_dims[k];
      size_t dim2 = dim * dim;
      const std::vector<size_t>& auts = this->_num_auts[k];

      // Copy upper diagonal matrix to the lower diagonal.
      for (size_t start=0, row=0; start<dim2; start+=dim+1, row++)
	{
	  int row_auts = auts[row];
	  for (size_t dst=start+dim, src=start+1, col=row+1; col<dim; src++, col++, dst+=dim)
	    {
	      if (matrix[src])
		{
		  int col_auts = auts[col];
		  if (col_auts == row_auts)
		    {
		      matrix[dst] = matrix[src];
		    }
		  else
		    {
		      assert( (matrix[src] * col_auts) % row_auts == 0 );

		      matrix[dst] = matrix[src] * col_auts / row_auts;
		    }
		}
	    }
	}

      // Move the matrix in the corresponding entry in the map.
      matrices[this->_conductors[k]] = std::move(hecke_matrices[k]);
    }
  return matrices;
}

// !! - TOOD - maybe it's better to return here already eigenvectors
template<typename R, size_t n>
inline std::vector< Matrix<int> >
Genus<R,n>::_decompositionRecurse(const Matrix<int> & V_basis,
				  const Integer<R> & p, size_t k) const
{
  // This will hold the bases of the irreducible spaces
  std::vector< Matrix<int> > decomp;
  
  if (V_basis.nrows() == 0)
    return decomp;

#ifdef DEBUG
  std::cerr << "Decomposing spaces of dimension " << this->_dims;
  std::cerr << "using T_" << p << "." << std::endl;
#endif

  // !! - TODO - check that results are stored and we don't
  // recompute for different values of k
  std::map<R,std::vector<int>> T_p_dense = heckeMatrixDense(p.num());
  
  Matrix<int> T_p(T_p_dense[k],this->_dims[k], this->_dims[k]);
  T_p = T_p.restrict(V_basis);
    
  UnivariatePoly<Z> f = T_p.charPoly();
  std::unordered_map< UnivariatePoly<Z>, size_t > fac = f.factor();

  for( std::pair< UnivariatePoly<Z>, size_t > fa : fac) {
    UnivariatePoly<Z> f = fa.first;
    size_t a = fa.second;
    
#ifdef DEBUG_LEVEL_FULL
    std::cerr << "Cutting out subspace using f(T_" << p;
    std::cerr << "), where f = " << f << "." << std::endl;
#endif

    Matrix<int> fT = f.evaluate(T_p);
    Matrix<int> W_basis = fT.kernel();

    if (a == 1)
      decomp.push_back(W_basis);
    else {
      Integer<R> q;
      if (W_basis.nrows() == V_basis.nrows())
	q = p.nextPrime();
      else
	q = R(2);
      std::vector< Matrix<int> > sub = this->_decompositionRecurse(W_basis, q, k);
      decomp.insert(decomp.end(), sub.begin(), sub.end());
    }
  }
  
  return decomp;
}

// !! TODO - support non-squarefree (when there are oldforms)
template<typename R, size_t n>
inline std::vector< Matrix<int> > Genus<R,n>::_decomposition(size_t k) const
{
  std::vector< Matrix<int> > decomp;
  if (this->_dims[k] == 0)
    return decomp;

  MatrixInt<int> M_basis = MatrixInt<int>::identity(this->_dims[k]);

  Integer<R> p = R(2);

  return this->_decompositionRecurse(M_basis, p, k);
}

template<typename R, size_t n>
inline std::map<R, std::vector< std::vector< NumberFieldElement<Z> > > >
Genus<R,n>::eigenvectors(void)
{
  std::map<R, std::vector< std::vector< NumberFieldElement<Z> > > > evecs;

  for (size_t k = 0; k < this->_conductors.size(); k++){
    std::vector< MatrixInt<int> > decomp = this->_decomposition(k);
    std::vector< std::vector< NumberFieldElement<Z> > > evecs_k;
    for (MatrixInt<int> T : decomp) {
      UnivariatePolyInt<Z> f = T.charPoly();
      std::shared_ptr< NumberField<Z> > K
	= std::make_shared< NumberField<Z> >(f);
      Matrix< NumberFieldElement<Z>, NumberField<Z> > T_K(K, T.nrows(), T.ncols());
      for (size_t row = 0; row < T.nrows(); row++)
	for (size_t col = 0; col < T.ncols(); col++)
	  T_K(row,col) = NumberFieldElement<Z>(K, T(row,col));
      NumberFieldElement<Z> lambda(K, UnivariatePolyRat<Z>::x());
      T_K -= lambda * Matrix< NumberFieldElement<Z> >::identity(T.nrows());
      Matrix< NumberFieldElement<Z>, NumberField<Z> > nullsp = T_K.kernel();
      std::vector< NumberFieldElement<Z> > vec = nullsp[0];
      evecs_k.push_back(vec);
    }
    evecs[this->_conductors[k]] = evecs_k;
  }
  return evecs;
}
