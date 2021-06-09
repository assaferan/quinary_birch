#ifndef __MATRIX_RING_H_
#define __MATRIX_RING_H_

#include "Ring.h"
#include "RingElement.h"

// Here R is derived from Ring and RElt is derived from RingElement

template<class RElt>
class Matrix;

template<class R, class RElt>
class MatrixRing : public virtual Ring< MatrixRing<R, RElt>, Matrix<RElt> >
{
  static_assert(std::is_base_of<Ring, R>::value);
  static_assert(std::is_base_of<RingElement, RElt>::value);
  
public:

  MatrixRing(std::shared_ptr<const R> ring) : _base(ring) {}
  
  inline std::shared_ptr<const MatrixRing<R> > getPtr() const override
  {return std::enable_shared_from_this< const MatrixRing<R> >::shared_from_this();}
  
  // producing the global constants of the ring
  inline Matrix<RElt> zero() const override
  { return Matrix<RElt>::zero(_base); }
    
  inline Matrix<RElt> one() const override
  { return Matrix<RElt>::one(_base); }

protected:

  std::shared_ptr<const R> _base;
};

#endif // __MATRIX_RING_H_
