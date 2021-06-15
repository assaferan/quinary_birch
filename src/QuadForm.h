#ifndef __QUAD_FORM_H_
#define __QUAD_FORM_H_

#include <ostream>

#include "birch.h"

template<class R, class Parent, size_t n>
std::ostream& operator<<(std::ostream&, const QuadForm<R,Parent,n> &);

template<class R, class Parent, size_t n>
class QuadForm
{
  public:
  typedef R SymVec[n*(n+1)/2];

  // c-tors
  QuadForm(std::shared_ptr<const Parent > ring) : _B(ring) {}
  
  // from a vector of n(n+1)/2 elements
  QuadForm(const SymVec& coeffs);
  QuadForm(const SquareMatrix<R,Parent,n> & B) : _B(B) {}
  
  // assignment
  inline QuadForm<R,Parent,n>& operator=(const QuadForm<R,Parent,n> &);

  // access
  R discriminant(void) const;

  inline bool operator==(const QuadForm<R,Parent,n>& q) const
  { return (this->_B == q._B); }

  inline bool operator!=(const QuadForm<R,Parent,n>& q) const
  {return !((*this)==q);}

  inline R evaluate(const Vector<R,Parent,n>& vec) const
  { R one = baseRing()->one(); R two = one + one; assert(!(two.isZero()));
    return Vector<R,Parent,n>::innerProduct(vec, (this->_B) * vec) / two; }

  inline const SquareMatrix<R,Parent,n> & bilinearForm() const
  { return this->_B; }

  std::shared_ptr<const Parent > baseRing() const {return _B.baseRing();}
  
protected:
  // _B = the matrix representing the
  // bilinear form Q(x+y)-Q(x)-Q(y) (so Q(v) = 1/2 B(v,v))
  SquareMatrix<R,Parent,n> _B;
 
};

#include "QuadForm.inl"

#endif // __QUAD_FORM_H_
