#ifndef __MATRIX_RING_H_
#define __MATRIX_RING_H_

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
  {return std::enable_shared_from_this< const IntegerRing<R> >::shared_from_this();}
  
  // producing the global constants of the ring
  inline Matrix<RElt> zero() const override;
  inline Matrix<RElt> one() const override;

protected:

  std::shared_ptr<const R> _base;
};

#endif // __MATRIX_RING_H_
