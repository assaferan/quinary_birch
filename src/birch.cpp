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
  /*
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

  
  std::vector<Z64_PrimeSymbol> symbols_64;
  Z64_PrimeSymbol p_64;
  std::vector<Z_PrimeSymbol> symbols;
  Z_PrimeSymbol p;
  std::vector<Z128_PrimeSymbol> symbols_128;
  Z128_PrimeSymbol p_128;
  
  
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

  std::set< Integer<Z64> > F_64;
  std::set<std::pair<Integer<Z64>, int> > F_ext_64;
    
  std::set<Integer<Z> > F;
  std::set<std::pair<Integer<Z>, int> > F_ext;
    
#ifdef DEBUG
  size_t I;

  for (std::vector<Z64_QuadForm<5> > genus : vec_64)
    {
      for (Z64_QuadForm<5> q : genus)
	{
	  std::cerr << q << std::endl;
	  std::cerr << q.discriminant() << std::endl;
	  Integer<Z64> det = q.invariants(F_64,I);
	  std::cerr << "det = " << det << std::endl;
	  std::cerr<< std::endl;
	  for (Integer<Z64> f : F_64)
	    std::cerr << f << " ";
	  std::cerr<< std::endl << I << std::endl << std::endl;
	  det = q.invariants(F_ext_64,I);
	  for (std::pair<Integer<Z64>,int> f : F_ext_64)
	    std::cerr << "Hasse(" << f.first << ") =  " << f.second << " ";
	  std::cerr << std::endl;
	}
    }
    
  for (std::vector<Z_QuadForm<5> > genus : vec)
    {

      for (Z_QuadForm<5> q : genus)
	{
	  std::cerr << q << std::endl;
	  std::cerr << q.discriminant() << std::endl;
	  Integer<Z> det = q.invariants(F,I);
	  std::cerr << "det = " << det << std::endl;
	  std::cerr<< std::endl;
	  for (Integer<Z> f : F)
	    std::cerr << f << " ";
	  std::cerr<< std::endl << I << std::endl << std::endl;
	  det = q.invariants(F_ext,I);
	  for (std::pair<Integer<Z>,int> f : F_ext)
	    std::cerr << "Hasse(" << f.first << ") =  " << f.second << " ";
	  std::cerr << std::endl;
	}

    }
#endif
    
  p_64.p = 61;
  p_64.power = 1;
  p_64.ramified = true;
  symbols_64.push_back(p_64);

  p_128.p = 61;
  p_128.power = 1;
  p_128.ramified = true;
  symbols_128.push_back(p_128);
    
  p.p = 61;
  p.power = 1;
  p.ramified = true;
  symbols.push_back(p);
  */

  Z64_QuadForm<3>::SymVec coeffs = {2,0,2,1,0,6};
  Z64_QuadForm<3> q(coeffs);

  // !! TODO - maybe use these to determine the spinor primes?
  std::set< Integer<Z64> > F;
  size_t I;
  Integer<Z64> det = q.invariants(F,I);

  // For now, using the discriminant to do these
  
  Integer<Z64> disc = q.discriminant();
  Integer<Z64>::FactorData facs = disc.factorization();

  std::vector<Z64_PrimeSymbol> symbols;
  Z64_PrimeSymbol symb;

  for (std::pair<Integer<Z64>, size_t> fa : facs) {
    symb.p = fa.first.num();
    symb.power = fa.second;
    symb.ramified = true;
    symbols.push_back(symb);
  }

  Z64_Genus<3> genus(q, symbols);

  std::map<Z64,size_t> dims = genus.dimensionMap();

  std::cout << "Dimensions are ";
  for (std::pair<Z64, size_t> dim : dims) {
    std::cout << dim.second << " ";
  }
  std::cout << std::endl;
  
  return 0;
}
