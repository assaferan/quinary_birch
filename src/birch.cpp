#include <iostream>

#include "birch.h"
#include "birch_util.h"
#include "testInteger.h"
#include "testRational.h"
#include "FpElement.h"
#include "Matrix.h"
#include "SquareMatrix.h"
#include "Vector.h"
#include "QuadFormInt.h"

std::ostream & operator<<(std::ostream & os, const Z128 & z)
{
  os << birch_util::convert_Integer<Z128, Z>(z);
  return os;
}

int main()
{
  testInteger<Z64> test1;
  testInteger<Z> test2;
  testInteger<Z128> test3;

  testRational<Z64> test4;
  testRational<Z> test5;
  testRational<Z128> test6;
  
  W64 seed = 1;
  std::shared_ptr< W16_Fp > GF = std::make_shared<W16_Fp>(3, seed);
  W16_FpElement a(GF);

  Matrix< W16_FpElement, W16_Fp > mat = Matrix< W16_FpElement, W16_Fp >::identity(GF,5);
  Vector< W16_FpElement, W16_Fp, 5> vec_fp(GF);
  SquareMatrix< W16_FpElement, W16_Fp, 5> sq_mat(GF);

  /*
  std::vector<Z64_PrimeSymbol> symbols_64;
  Z64_PrimeSymbol p_64;
  std::vector<Z_PrimeSymbol> symbols;
  Z_PrimeSymbol p;
  std::vector<Z128_PrimeSymbol> symbols_128;
  Z128_PrimeSymbol p_128;
  */
  
  Z64_QuadForm<3>::SymVec coeffs_64 = {2,1,2,1,1,2};
  Z_QuadForm<3>::SymVec coeffs = {Z(2),Z(1),Z(2),Z(1),Z(1),Z(2)};

  Z64_QuadForm<3> q0_64(coeffs_64);
    
#ifdef DEBUG
  const Z64_SquareMatrix<3> & B_64 = q0_64.bilinearForm();
  std::cerr << "B_64 = " << B_64 << std::endl;
#endif
    
  Z_QuadForm<3> q0(coeffs);
    
#ifdef DEBUG
  const Z_SquareMatrix<3> & B = q0.bilinearForm();
  std::cerr << "B = " << B << std::endl;
#endif
    
  std::vector<std::vector<Z64_QuadForm<5> > >
    vec_64 = Z64_QuadForm<5>::getQuinaryForms(61);

  std::vector<std::vector<Z128_QuadForm<5> > >
    vec_128 = Z128_QuadForm<5>::getQuinaryForms(61);
    
  std::vector<std::vector<Z_QuadForm<5> > >
    vec = Z_QuadForm<5>::getQuinaryForms(61);

  return 0;
}
