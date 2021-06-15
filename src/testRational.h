#ifndef __TEST_RATIONAL_H_
#define __TEST_RATIONAL_H_

#include "Rational.h"

template<typename R>
class testRational
{
public:
  testRational();
  
protected:

  static void _testConstructor(const R & num, const R & denom);
  static void _testAdd(const R & num1, const R & denom1,
		       const R & num2, const R & denom2);
  static void _testSub(const R & num1, const R & denom1,
		       const R & num2, const R & denom2);
  static void _testMul(const R & num1, const R & denom1,
		       const R & num2, const R & denom2);
  static void _testDiv(const R & num1, const R & denom1,
		       const R & num2, const R & denom2);
  static void _testPow(const R & num1, const R & denom1,
		       const Z64 & e);
  static void _testInverse(const R & num, const R & denom);
  static void _testZero(const R & denom);
  static void _testOne(const R & val);
  static void _testGcd(const R & num1, const R & denom1,
		       const R & num2, const R & denom2);
  static void _testEuclid(const R & num1, const R & denom1,
			  const R & num2, const R & denom2);
  
  
};

#include "testRational.inl"

#endif //  __TEST_RATIONAL_H_
