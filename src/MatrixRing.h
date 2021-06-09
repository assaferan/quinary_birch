#ifndef __MATRIX_RING_H_
#define __MATRIX_RING_H_

#include "Ring.h"
#include "RingElement.h"

// Here R is derived from Ring and RElt is derived from RingElement

template<class RElt, class R>
class Matrix;

template<class RElt, class R>
class MatrixRing : public virtual Ring< MatrixRing<RElt,R>, Matrix<RElt, R> >
{
  static_assert(std::is_base_of<Ring, R>::value);
  static_assert(std::is_base_of<RingElement, RElt>::value);
  
public:

  MatrixRing(std::shared_ptr<const R> ring, size_t n) : _base(ring) {}
  
  inline std::shared_ptr<const MatrixRing<RElt,R> > getPtr() const override
  {return std::enable_shared_from_this< const MatrixRing<RElt,R> >::shared_from_this();}
  
  // producing the global constants of the ring
  inline Matrix<RElt,R> zero() const override
  { return Matrix<RElt,R>::zero(_base); }
    
  inline Matrix<RElt,R> one() const override
  { return Matrix<RElt,R>::one(_base); }

  inline std::shared_ptr<const R> baseRing() const {return _base;}

protected:

  std::shared_ptr<const R> _base;
};

#endif // __MATRIX_RING_H_
