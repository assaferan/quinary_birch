#ifndef __QUAD_FORM_FP_H_
#define __QUAD_FORM_FP_H_

#include "birch.h"
#include "QuadForm.h"

template<typename R, typename S, size_t n>
class QuadFormFp : public QuadForm< FpElement<R, S> , Fp<R,S>, n>
{
public:
  QuadFormFp(const SquareMatrixFp<R,S,n> & mat) : 
    QuadForm< FpElement<R, S>, Fp<R,S>, n>(mat)
  {}

  QuadFormFp(const typename QuadFormInt<R,n>::SymVec& vec,
	     std::shared_ptr<const Fp<R,S>> GF) :
    QuadForm<FpElement<R, S> ,Fp<R,S>, n>(GF->mod(vec))
  {}

  inline const std::shared_ptr<const Fp<R,S>>& field(void) const
  { return this->_B.baseRing();}

  using QuadForm<FpElement<R,S>,Fp<R,S>,n>::bilinearForm;
  using QuadForm<FpElement<R,S>,Fp<R,S>,n>::discriminant;
  
  inline FpElement<R,S> evaluate(const VectorFp<R,S,n>& v) const;
  
  inline Integer<R> evaluate(const Vector<Integer<R>,IntegerRing<R>,n>& vec) const {
    VectorFp<R,S,n> v = this->_B.baseRing()->mod(vec);
    return (this->evaluate(v)).lift();
  }
  
  bool isotropicVector(VectorFp<R,S,n> &,
		       size_t start = 0,
		       bool deterministic = false) const;

  void decompose(SquareMatrixFp<R,S,n> &,
		 SquareMatrixFp<R,S,n> &,
		 bool deterministic = false) const;

protected:

  // To avoid unnecessary computation, we encode each of the three 2-isotropic
  // vectors as a coordinate of the return vector. Special care must be taken
  // to obtain the actual isotropic vectors when needed.
 
  bool _isotropicVector_p2(VectorFp<R,S,n> &, size_t start = 0) const;

  FpElement<R,S> _evaluate_p2(const VectorFp<R,S,n>& v) const;

  void _hyperbolizeForm(SquareMatrixFp<R,S,n> &,
			SquareMatrixFp<R,S,n> &,
			bool deterministic = false,
			size_t start = 0) const;
  
  void _splitHyperbolicPlane(const VectorFp<R,S,n> &,
			     SquareMatrixFp<R,S,n> &,
			     SquareMatrixFp<R,S,n> &,
			     size_t start = 0) const;
};

#include "QuadFormFp.inl"

#endif // __QUAD_FORM_FP_H_
