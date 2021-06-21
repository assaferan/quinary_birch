#include <pmmintrin.h>
#include <immintrin.h>

#include "birch.h"
#include "SquareMatrixInt.h"

// these are specializations we were forces to produce
template<>
SquareMatrixInt<Z64,2> SquareMatrixInt<Z64,2>::operator*(const SquareMatrixInt<Z64,2>& other) const
{
  SquareMatrixInt<Z64,2> prod;

  for (size_t i = 0; i < 2; i++)
    for (size_t j = 0; j < 2; j++) {
      prod._mat[i][j] = 0;
      for (size_t k = 0; k < 2; k++)
	prod._mat[i][j] += this->_mat[i][k]*other._mat[k][j];
    }
  
  return prod;
}

template<>
SquareMatrixInt<Z64,3> SquareMatrixInt<Z64,3>::operator*(const SquareMatrixInt<Z64,3>& other) const
{
  SquareMatrixInt<Z64,3> prod;

  for (size_t i = 0; i < 3; i++)
    for (size_t j = 0; j < 3; j++) {
      prod._mat[i][j] = 0;
      for (size_t k = 0; k < 3; k++)
	prod._mat[i][j] += this->_mat[i][k]*other._mat[k][j];
    }
  
  return prod;
}

template<>
SquareMatrixInt<Z,3> SquareMatrixInt<Z,3>::operator*(const SquareMatrixInt<Z,3>& other) const
{
  SquareMatrixInt<Z,3> prod;

  for (size_t i = 0; i < 3; i++)
    for (size_t j = 0; j < 3; j++) {
      prod._mat[i][j] = 0;
      for (size_t k = 0; k < 3; k++)
	prod._mat[i][j] += this->_mat[i][k]*other._mat[k][j];
    }
  
  return prod;
}

// matrix multiplication is a major bottleneck, hence we attempt to optimize it here
template<>
SquareMatrixInt<Z64,4> SquareMatrixInt<Z64,4>::operator*(const SquareMatrixInt<Z64,4>& B) const
{
  SquareMatrixInt<Z64,4> C;
  const SquareMatrixInt<Z64,4>& A = *this;
  
  const __m256i BCx = _mm256_loadu_si256((const __m256i_u*)&B._mat[0]);
  const __m256i BCy = _mm256_loadu_si256((const __m256i_u*)&B._mat[1]);
  const __m256i BCz = _mm256_loadu_si256((const __m256i_u*)&B._mat[2]);
  const __m256i BCw = _mm256_loadu_si256((const __m256i_u*)&B._mat[3]);

  const Z64* leftRowPointer = &(A._mat[0][0]);
  Z64* resultRowPointer = &(C._mat[0][0]);

  for (size_t i = 0; i < 4; ++i, leftRowPointer += 4, resultRowPointer += 4) {
    __m256i ARx = _mm256_set1_epi64x(leftRowPointer[0]);
    __m256i ARy = _mm256_set1_epi64x(leftRowPointer[1]);
    __m256i ARz = _mm256_set1_epi64x(leftRowPointer[2]);
    __m256i ARw = _mm256_set1_epi64x(leftRowPointer[3]);

    __m256i X = ARx * BCx;
    __m256i Y = ARy * BCy;
    __m256i Z = ARz * BCz;
    __m256i W = ARw * BCw;

    __m256i S = X+Y+Z+W;
    _mm256_store_si256((__m256i_u*)resultRowPointer, S);
  }
  
  return C;
}
