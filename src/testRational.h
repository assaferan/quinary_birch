#ifndef __TEST_RATIONAL_H_
#define __TEST_RATIONAL_H_

#include "Rational.h"

template<typename R>
class testRational
{
public:
  testRational();
  
protected:

  static void testConstructor(const R & num, const R & denom);
  static void testAdd(const R & num1, const R & denom1,
		      const R & num2, const R & denom2);
  static void testSub(const R & num1, const R & denom1,
		      const R & num2, const R & denom2);
  static void testMul(const R & num1, const R & denom1,
		      const R & num2, const R & denom2);
  static void testDiv(const R & num1, const R & denom1,
		      const R & num2, const R & denom2);
  static void testPow(const R & num1, const R & denom1,
		      const R & num2, const R & denom2);
  static void testZero(const R & denom);
  static void testOne(const R & val);
  static void testGcd(const R & num1, const R & denom1,
		      const R & num2, const R & denom2);
  static void testEuclid(const R & num1, const R & denom1,
			 const R & num2, const R & denom2);
  
  
};

#include "testRational.inl"

#endif //  __TEST_RATIONAL_H_
