#ifndef __SQUARE_MATRIX_INT_H_
#define __SQUARE_MATRIX_INT_H_

// We create this one to mkae matrix operations more efficient when the underlying class is integral

#include "SquareMatrix.h"

template <typename R, size_t n>
class SquareMatrixInt : public virtual SquareMatrix< Integer<R>, IntegerRing<R>, n>
{
public:
   // c-tors
  SquareMatrixInt() : _base(std::make_shared<const IntegerRing<R> >()) {}

  SquareMatrixInt(const R mat[n][n]);
  SquareMatrixInt(const SquareMatrixInt<R,n> & mat);

  SquareMatrixInt operator*(const SquareMatrixInt &) const;
};

#include "SquareMatrixInt.inl"

#endif // __SQUARE_MATRIX_INT_H_
