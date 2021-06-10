template<typename R, typename S, size_t n>
inline FpElement<R,S> QuadFormFp<R,S,n>::evaluate(const VectorFp<R,S,n>& v) const {
  R p = this->_B.baseRing()->prime();
  if (p == 2) return this->evaluate_p2(v);
  VectorFp<R,S,n> Bv = (this->bilinearForm()) * v;
  FpElement<R,S> two(GF, 2);
  return VectorFp<R,S,n>::innerProduct(v, Bv)/two;
} 

template<typename R, typename S, size_t n>
inline void QuadFormFp<R,S,n>::splitHyperbolicPlane(const VectorFp<R,S,n>& vec,
						    SquareMatrixFp<R,S,n>& gram,
						    SquareMatrixFp<R,S,n> & basis,
						    size_t start) const
{
  // !! TODO - see if we can get rid of this.
  // SquareMatrixFp<R, S, n> orig_basis = basis;
  // The change of basis which preserving the isometry.
  basis.setIdentity();
  // Make a copy of the Gram matrix.
  gram = this->bilinearForm();
  std::shared_ptr<const Fp<R,S> > GF = gram.baseRing();
  R p = GF->prime();
  // Set the diagonal entries to zero when in characteristic 2.
  // This is because we are decomposing the associated bilinear form
  if (p == 2) {
    for (size_t i = start; i < n; i++) gram(i,i) = 0;
  }
  SquareMatrixFp<R,S,n> original_gram = gram;
  // Find the pivot of the specified vector.
  size_t pivot = start;
  while (vec[pivot] == 0) pivot++;

  assert(pivot < n);
  
  // Perform the necessary basis changes so that vec becomes the first
  //  basis vector.
  basis.multiplyRow(pivot, vec[pivot]);
  gram.multiplyCol(pivot, vec[pivot]);
  gram.multiplyRow(pivot, vec[pivot]);
  for (size_t i = pivot+1; i < n; i++) {
    basis.addRow(pivot, i, vec[i]);
    gram.addCol(pivot, i, vec[i]);
    gram.addRow(pivot, i, vec[i]);
  }
  basis.swapRows(start, pivot);
  gram.swapCols(start, pivot);
  gram.swapRows(start, pivot);

  bool is_in_radical = true;
  
  // Find a basis vector which is not orthogonal to our isotropic vector.
  size_t idx = start;
  for (; idx < n; idx++)
    if (gram(start, idx) != 0) {
      is_in_radical = false;
      break;
    }

  // If the first row is entirely zero, then this vector belongs to the
  //  radical of the form.
  
  if (is_in_radical) {
    if (p == 2) {
      // Recover the quadratic form along the diagonal.
      for (size_t i = start; i < n; i++)
	gram(i, i) = this->evaluate(basis[i]);
    }
    return;
  }
  
  // Swap this basis vector with the second basis.
  basis.swapRows(start+1, idx);
  gram.swapCols(start+1, idx);
  gram.swapRows(start+1, idx);

  // Normalize the second basis vector so the (0,1)-entry is 1.
  FpElement<R,S> scalar = gram(start,start+1).inverse();
  basis.multiplyRow(start+1, scalar);
  gram.multiplyCol(start+1, scalar);
  gram.multiplyRow(start+1, scalar);

  // Determine the appropriate scalar for clearing out the (1,1)-entry.
  if (p == 2)
    scalar = this->evaluate(basis[start+1]);
  else {
    FpElement<R,S> two(GF,2);
    scalar = -gram(start+1,start+1) / two;
  }

  // Clear the (1,1)-entry in the Gram matrix.
  basis.addRow(start+1, start, scalar);
  gram.addCol(start+1, start, scalar);
  gram.addRow(start+1, start, scalar);

  // Clear the remaining entries in the Gram matrix.
  for (size_t i = start+2; i < n; i++) {
    // Clear first row/column.
    scalar = -gram(start, i);
    basis.addRow(i, start+1, scalar);
    gram.addCol(i, start+1, scalar);
    gram.addRow(i, start+1, scalar);

    // Clear second row/column.
    scalar = -gram(start+1, i);
    basis.addRow(i, start, scalar);
    gram.addCol(i, start, scalar);
    gram.addRow(i, start, scalar);
  }

#ifdef DEBUG
  SquareMatrixFp<R,S,n> tmp = basis * original_gram;
  tmp = tmp * basis.transpose();
  assert(tmp == gram);
#endif

  // In characteristic 2, we need to recover the diagonal entries by
  //  evaluating the basis via the quadratic form.
  
  if (p == 2)
    for (size_t i = start; i < n; i++)
      gram(i,i) = this->evaluate(basis[i]);

  // basis = basis * orig_basis;
  
  return;
}

template<typename R, typename S, size_t n>
inline void QuadFormFp<R,S,n>::hyperbolizeForm(SquareMatrixFp<R,S,n> & gram,
					       SquareMatrixFp<R,S,n> & basis,
					       bool deterministic,
					       size_t start) const
{
  std::shared_ptr<const Fp<R,S> > GF = _B.baseRing();
  VectorFp<R,S,n> vec(GF);
  bool found = this->isotropicVector(vec, start, deterministic);
  size_t dim = n - start;
  
  // The space is anisotropic.
  if (!found) {
    SquareMatrixFp<R,S,n> originalGram = gram;
    FpElement<R,S> scalar;
    if (dim == 1) {
      // Check if the (0,0)-entry is a square.
      if (gram(start,start).isSquare()) {
	// If so, make it a 1.
	scalar = gram(start,start).sqrt().inverse();
	basis.multiplyRow(start, scalar);
	gram.multiplyCol(start, scalar);
	gram.multiplyRow(start, scalar);
      }
      return;
    }
    if (GF->prime() == 2) {
      // Make the (0,0)-entry equal to 1.
      assert(gram(start,start).isSquare());

      scalar = gram(start,start).sqrt().inverse();
      basis.multiplyRow(start, scalar);
      gram.multiplyCol(start, scalar);
      gram.multiplyRow(start, scalar);

      // Make the (0,1)-entry equal to 1.
      scalar = gram(start,start+1).inverse();
      basis.multiplyRow(start+1, scalar);
      gram.multiplyCol(start+1, scalar);
      gram.multiplyRow(start+1, scalar);

      return;
    }
    
    // Clear the (0,1)-entry.
    scalar = -gram(start,start+1)/gram(start,start);
    basis.addRow(start+1, start, scalar);
    gram.addCol(start+1, start, scalar);
    gram.addRow(start+1, start, scalar);

    // If the (1,1)-entry is a square, make it the first entry.
    if (gram(start+1,start+1).isSquare()) {
      basis.swapRows(start,start+1);
      gram.swapRows(start,start+1);
      gram.swapCols(start,start+1);
    }

    bool is_square[2];
    for (size_t i = 0; i < 2; i++) {
      // Check if the (i,i)-entry is a square then clear it, if so.
      is_square[i] = gram(start+i,start+i).isSquare();
      if (is_square[i]) {
	scalar = gram(start+i,start+i).sqrt().inverse();
	basis.multiplyRow(start+i, scalar);
	gram.multiplyCol(start+i, scalar);
	gram.multiplyRow(start+i, scalar);
      }
    }

    // If neither are squares, make them -1 (note that this occurs
    //  if and only if -1 is not a square).
    if ((!is_square[0]) && (!is_square[1])) {
      for (size_t i = 0; i < 2; i++) {
	assert((-gram(start+i,start+i)).isSquare());
	scalar = (-gram(start+i,start+i)).sqrt().inverse();
	basis.multiplyRow(start+i, scalar);
	gram.multiplyCol(start+i, scalar);
	gram.multiplyRow(start+i, scalar);
      }
    }

    return;
  }

  assert(this->evaluate(vec) == 0);

  // Attempt to split a hyperbolic plane from the form.
  splitHyperbolicPlane(vec, gram, basis, start);

  // Determine how many dimensions we need to split off.
  size_t lower_dim = gram[0].isZero() ? 1 : 2;

  if (dim > lower_dim) {
    // Split the hyperbolic plane from the form.
    QuadFormFp<R,S,n> q_split(gram);
    // !! TODO - check maybe we have to replace basis here
    SquareMatrixFp<R,S,n> newbasis(GF);
    newbasis.setIdentity();
    q_split.hyperbolizeForm(gram, newbasis, deterministic, start + lower_dim);
    basis = newbasis * basis;
  }

#ifdef DEBUG_LEVEL_FULL
  std::cerr << "After hyperbolize_form with start = " << start << "." << std::endl;
  std::cerr << "Resulting gram matrix is " << std::endl;
  gram.prettyPrint(std::cerr);
  std::cerr << ", ";
  std::cerr << "Resulting basis is " << std::endl;
  basis.prettyPrint(std::cerr);
  std::cerr << std::endl;
#endif 
  
  return;
}

template<typename R, typename S, size_t n>
inline void QuadFormFp<R,S,n>::decompose(SquareMatrixFp<R,S,n> & gram,
					 SquareMatrixFp<R,S,n> & basis,
					 bool deterministic) const
{
  basis.setIdentity();
  hyperbolizeForm(gram, basis, deterministic);

#ifdef DEBUG_LEVEL_FULL
  std::cerr << "After hyperbolize_form." << std::endl;
  std::cerr << "Resulting gram matrix is " << gram << ", ";
  std::cerr << "Resulting basis is " << basis << std::endl;

  std::shared_ptr< const Fp<R,S> > GF = _B.baseRing();
  
  // Verify that everyhing we've done is correct.
  SquareMatrixFp<R,S,n> temp1(GF);
  SquareMatrixFp<R,S,n> temp2(GF);
  if (GF->prime() == 2) {
    // Verify that the basis evaluates correctly on the form.
    for (size_t i = 0; i < n; i++) {
      assert(this->evaluate(basis[i]) == gram(i,i));
    }
    // Zero out the diagonal to verify the bilinear form is correct.
    // Recall that in characteristic 2, the associated bilinear form
    // is b_q(x,y) = q(x+y) - q(x) - q(y), with zeros on the diagonal
    temp1 = gram;
    temp2 = this->bilinearForm();
    for (size_t i = 0; i < n; i++) {
      temp1(i,i) = 0;
      temp2(i,i) = 0;
    }
  }
  else {
    temp1 = gram;
    temp2 = this->bilinearForm();
  }
  // Verify that the bilinear forms are similar.
  assert(basis * temp2 * basis.transpose() == temp1);
#endif

  // Let's bubble the basis vectors which belong to the radical to the
  //  end of the basis list.
  size_t rad = 0;
  size_t pos = n;
  while (pos >= 1) {
    if (gram[pos-1].isZero()) {
      rad++;
      for (size_t i = pos; i < n; i++) {
	basis.swapRows(i-1,i);
	gram.swapRows(i-1,i);
	gram.swapCols(i-1,i);
      }
    }
    pos--;
  }
  // Let's put the hyperbolic planes in our standard antidiagonal form.

  // The upper index of the hyperbolic space.
  size_t upper = n + 1 - rad;
  do {
    upper--;
  } while ((upper >= 1) && (gram(upper-1,upper-1) != 0));

  // Indices of the basis vectors we'll be swapping.
  size_t i = 1;
  size_t j = upper;
  // Keep swapping basis vectors until j is less than or equal to i. Note
  //  that if there are no hyperbolic planes, this does nothing.
  while (i < j) {
    basis.swapRows(i-1,j-2);
    gram.swapCols(i-1,j-2);
    gram.swapRows(i-1,j-2);
    i += 2;
    j -= 2;
  }
  // Since we did everything with row vectors, we need to transpose the
  //  basis, so that the subsequent code that utilizes it doesn't break.
  basis = basis.transpose();
  return;
}

// !! TODO - all F2 operations (including this one) can be made faster
template<typename R, typename S, size_t n>
inline FpElement<R,S> QuadFormFp<R,S,n>::evaluate_p2(const VectorFp<R,S,n>& v)
  const
{
  std::shared_ptr< const Fp<R,S> > GF = this->_B.baseRing();
  FpElement<R,S> val(GF, 0);
  for (size_t i = 0; i < n; i++) {
    val += this->_B(i,i) * v[i];
    for (size_t j = i+1; j < n; j++)
      val += this->_B(i,j) * v[i] * v[j];
  }
  return val;
}

template<typename R, typename S, size_t n>
bool QuadFormFp<R,S,n>::isotropicVector_p2(VectorFp<R,S,n> & vec,
					   size_t start) const
{
  FpElement<R,S> g;
  // If we can find a pair of orthogonal basis vectors,
  //  we can easily construct an isotropic vector.
  for (size_t i = start; i < n-1; i++)
    for (size_t j = i+1; j < n; j++) {
      if (this->bilinearForm()(i,j) == 0) {
	g = this->bilinearForm()(j,j) / this->bilinearForm()(i,i);
	assert(g.isSquare());
	g = g.sqrt();
	vec[i] = g;
	vec[j] = 1;
	return true;
      }
    }

  size_t dim = n - start;
  if (dim == 1) {
    if (this->bilinearForm()(start,start) == 0) {
      vec[start] = 1;
      return true;
    }
    return false;
  }
  
  if (dim == 2) {
    for (size_t i = 0; i < 2; i++)
      if (this->bilinearForm()(start+i,start+i) == 0) {
	vec[start+i] = 1;
	return true;
      }
    FpElement<R,S> a = this->_B(start,start);
    FpElement<R,S> b = this->_B(start,start+1);
    FpElement<R,S> c = this->_B(start+1,start+1);
    if (b == 0) {
      vec[start] = 1;
      vec[start+1] = 1;
      return true;
    }
    // In this case a = b = c = 1, so the form is anisotropic
    return false;
  }
  assert (dim >= 3);
  // Otherwise, while the formulation is a bit more
  //  complicated, we can produce an isotropic vector
  //  by taking a linear combination of the first three
  //  basis vectors as follows:

  // Convenient references.
  FpElement<R,S> a = this->bilinearForm()(start,start);
  FpElement<R,S> b = this->bilinearForm()(start+1,start+1);
  FpElement<R,S> c = this->bilinearForm()(start+2,start+2);
  FpElement<R,S> d = this->bilinearForm()(start+1,start+2);
  FpElement<R,S> e = this->bilinearForm()(start,start+2);
  FpElement<R,S> f = this->bilinearForm()(start,start+1);

  g = (b*e*e/f/f + c + e*d/f)/a;
  assert(g.isSquare());
  vec[start] = g;
  vec[start+1] = e/f;
  vec[start+2] = 1;

  return true;
}


template<typename R, typename S, size_t n>
inline bool QuadFormFp<R,S,n>::isotropicVector(VectorFp<R,S,n> & vec,
					       size_t start,
					       bool deterministic) const
{
  // Check the diagonal
  for (size_t i = start; i < n; i++)
    if (this->_B(i,i) == 0) {
      vec[i] = 1;
      return true;
    }

  size_t dim = n - start;

  // if the one vector was isotropic, we would have found it above.
  if (dim == 1)
    return false;

  std::shared_ptr<const Fp<R,S> > GF = this->_B.baseRing();
  
  if (GF->prime() == 2) return this->isotropicVector_p2(vec, start);

  if (dim == 2) {
    FpElement<R,S> a = this->_B(start,start);
    FpElement<R,S> b = this->_B(start,start+1);
    FpElement<R,S> c = this->_B(start+1,start+1);
    // The form is isotropic if and only if b^2-ac is a square.
    FpElement<R,S> d = b*b-a*c;
    // If not a square, this form is anisotropic.
    if (!d.isSquare()) {
      return false;
    }
    // Since a ne 0 and the form is isotropic, we're done.
    d = d.sqrt();
    vec[start] = -((b+d)/a);
    vec[start+1] = 1;
    return true;
  }

  assert(dim >= 3);

  // isometry on the submatrix of 3 first variables
  SquareMatrixFp<R,S,3> basis = SquareMatrixFp<R,S,3>::identity(GF);
  SquareMatrixFp<R,S,3> subM(GF);
  for (size_t i = 0; i < 3; i++)
    for (size_t j = 0; j < 3; j++)
      subM(i,j) = this->bilinearForm()(start+i,start+j);

  FpElement<R,S> scalar;
  
  // clear the off-diagonal entries
  for (size_t i = 0; i < 2; i++)
    for (size_t j = i+1; j < 3; j++) {
      scalar = -subM(i,j) / subM(i,i);
      subM.addCol(j, i, scalar);
      subM.addRow(j, i, scalar);
      basis.addRow(j, i, scalar);
      if (subM(j,j) == 0) {
	for (size_t k = 0; k < 3; k++)
	  vec[start+k] = basis(j,k);
	return true;
      }
    }

  // Check if the first two variables alone are isotropic.
  FpElement<R,S> d = -subM(0,0)*subM(1,1);
  if (d.isSquare()) {
    d = d.sqrt();
    for (size_t k = 0; k < 3; k++)
      vec[start+k] = basis(0,k) +
	(this->bilinearForm()(start,start)/d) * basis(1,k);
    return true;
  }

  if (deterministic) {
    R p = GF->prime();
    // The quadratic form over three variables.
    QuadFormFp<R,S,3> Q(subM);
    VectorFp<R,S,3> v(GF);
    for (R x = 0; x < p; x++)
      for (R y = 0; y < p; y++) {
	v[0] = x;
	v[1] = y;
	v[2] = 1;
	if (Q.evaluate(v) == 0) {
	  // Found an isotropic vector, return it.
	  for (size_t j = 0; j < 3; j++) {
	    vec[start+j] = v[0]*basis(0,j) + v[1]*basis(1,j) + basis(2,j);
	  }
	  return true;
	}
      }
  }

  // If we're fine with a probabilitistic means of finding
  //  isotropic vectors, we can find them much faster.
  
  FpElement<R,S> a = subM(0,0);
  FpElement<R,S> b = subM(1,1);
  FpElement<R,S> c = subM(2,2);

  VectorFp<R,S,2> v(GF);
  bool nonzero;
  do {
    do {
      do {
	for (size_t i = 0; i < 2; i++)
	  v[i] = GF->random();
      } while ((v[0] == 0) && (v[1] == 0));
      d = -(a*v[0]*v[0] + b*v[1]*v[1])/c;
    } while (!d.isSquare());
    
    d = d.sqrt();
    nonzero = false;
    for (size_t j = 0; j < 3; j++) {
      vec[start+j] = v[0]*basis(0,j) + v[1]*basis(1,j) + d*basis(2,j);
      nonzero = nonzero || (vec[start+j] != 0);
    }
  } while (!nonzero);
  return true;
}
