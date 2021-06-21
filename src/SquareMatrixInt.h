#ifndef __SQUARE_MATRIX_INT_H_
#define __SQUARE_MATRIX_INT_H_

// We create this one to mkae matrix operations more efficient when the underlying class is integral

#include "SquareMatrix.h"

template <typename R, size_t n>
class SquareMatrixInt : public virtual SquareMatrix< Integer<R>, IntegerRing<R>, n>
{
public:
   // c-tors
  SquareMatrixInt()
    : SquareMatrix< Integer<R>, IntegerRing<R>, n>(std::make_shared<const IntegerRing<R> >())
  {}
  SquareMatrixInt(std::shared_ptr<const IntegerRing<R> > base_ring)
    : SquareMatrix< Integer<R>, IntegerRing<R>, n>(base_ring)
  {}
  
  SquareMatrixInt(const R mat[n][n]);
  SquareMatrixInt(const SquareMatrixInt<R,n> & mat);
  SquareMatrixInt(const SquareMatrix<Integer<R>,IntegerRing<R>,n> & mat)
    : SquareMatrix< Integer<R>, IntegerRing<R>, n>(mat)
  {}

  SquareMatrixInt<R,n> operator*(const SquareMatrixInt &) const;
  SquareMatrixInt<R,n> transpose(void) const;

};

#include "SquareMatrixInt.inl"

#endif // __SQUARE_MATRIX_INT_H_
